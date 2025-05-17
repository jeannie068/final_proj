#include "optimizer.hpp"

#include <unordered_map>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <memory>

/**
 * @brief Constructor
 */
Optimizer::Optimizer(const ChipInfo& chip_info,
                   const std::vector<CellType>& cell_types,
                   const std::vector<Cell>& cells,
                   const AlgorithmParams& params)
    : chip_info(chip_info), cell_types(cell_types), cells(cells), params(params) {
    
    // Organize cells by row
    cells_by_row.resize(chip_info.num_rows);
    for (const Cell& cell : cells) {
        int row_idx = cell.getRowIndex(chip_info.row_height);
        if (row_idx >= 0 && row_idx < chip_info.num_rows) {
            // Create a copy of the cell and store a pointer
            Cell* cell_ptr = new Cell(cell);
            cells_by_row[row_idx].push_back(cell_ptr);
        }
    }
    
    // Sort cells in each row by initial x-coordinate
    for (auto& row : cells_by_row) {
        std::sort(row.begin(), row.end(), [](const Cell* a, const Cell* b) {
            return a->initial_x < b->initial_x;
        });
    }
}

/**
 * @brief Run the optimization algorithm
 */
Solution Optimizer::run() {
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Solve the problem by processing triple-row subproblems
    std::vector<Staple> prev_staples;
    
    for (int row = 0; row < chip_info.num_rows - 2; row += 2) {
        int row_end = std::min(row + 3, chip_info.num_rows);
        std::vector<Staple> new_staples = solveTripleRow(row, row_end, prev_staples);
        
        // Add new staples to the list
        inserted_staples.insert(inserted_staples.end(), new_staples.begin(), new_staples.end());
        
        // Update prev_staples for the next iteration
        prev_staples = new_staples;
    }
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    
    std::cout << "Optimization completed in " << duration.count() << " ms" << std::endl;
    
    // Create solution object
    Solution solution;
    
    // Copy refined cells
    for (const auto& row : cells_by_row) {
        for (const Cell* cell_ptr : row) {
            solution.refined_cells.push_back(*cell_ptr);
        }
    }
    
    // Add inserted staples
    solution.inserted_staples = inserted_staples;
    
    // Update statistics
    solution.updateStats();
    
    // Clean up allocated memory
    for (auto& row : cells_by_row) {
        for (Cell* cell_ptr : row) {
            delete cell_ptr;
        }
    }
    
    return solution;
}

/**
 * @brief Solve a triple-row optimization problem
 */
std::vector<Staple> Optimizer::solveTripleRow(int row_start, int row_end, 
                                           const std::vector<Staple>& prev_staples) {
    std::cout << "Solving triple-row problem for rows " << row_start << " to " << (row_end - 1) << std::endl;
    
    // Extract cells for each row in the triple-row problem
    std::vector< std::vector<Cell*> > cells_in_rows;
    for (int r = row_start; r < row_end; r++) {
        if (r < static_cast<int>(cells_by_row.size())) {
            cells_in_rows.push_back(cells_by_row[r]);
        } else {
            cells_in_rows.push_back(std::vector<Cell*>());
        }
    }
    
    // Ensure we have exactly 3 rows (pad with empty rows if necessary)
    while (cells_in_rows.size() < 3) {
        cells_in_rows.push_back(std::vector<Cell*>());
    }
    
    // Initialize the DAG
    DPNode* source_node = initializeDAG(row_start, row_end, cells_in_rows);
    
    // Create node lookup table and queue
    std::unordered_map<CompactState, DPNode*, CompactStateHasher> node_lookup_table;
    std::queue<DPNode*> queue;
    
    // Add source node to queue and lookup table
    queue.push(source_node);
    node_lookup_table[source_node->state] = source_node;
    
    // Process each site from left to right
    for (int site = 0; site <= chip_info.total_sites; site++) {
        processSite(site, queue, node_lookup_table, cells_in_rows, prev_staples, row_start);
        
        // Clear lookup table after processing each site to save memory
        node_lookup_table.clear();
    }
    
    // Find the best node at the rightmost site
    DPNode* best_node = nullptr;
    int max_benefit = -1;
    
    for (auto it = node_lookup_table.begin(); it != node_lookup_table.end(); ++it) {
        DPNode* node = it->second;
        if (node->state.site == chip_info.total_sites) {
            for (int case_idx = 0; case_idx < 5; case_idx++) {
                if (node->benefit[case_idx] > max_benefit) {
                    max_benefit = node->benefit[case_idx];
                    best_node = node;
                }
            }
        }
    }
    
    // Backtrack to get the solution
    std::vector<Staple> inserted_staples;
    if (best_node != nullptr) {
        inserted_staples = backtrack(best_node, cells_in_rows, row_start);
    }
    
    // Clean up allocated nodes
    for (auto it = node_lookup_table.begin(); it != node_lookup_table.end(); ++it) {
        delete it->second;
    }
    
    std::cout << "Triple-row optimization completed. Inserted " << inserted_staples.size() 
              << " staples." << std::endl;
    
    return inserted_staples;
}

/**
 * @brief Initialize the DAG for a triple-row problem
 */
DPNode* Optimizer::initializeDAG(int row_start, int row_end, 
                               const std::vector< std::vector<Cell*> >& cells_in_rows) {
    // Create source node with initial state
    CompactState initial_state;
    initial_state.site = 0;
    
    for (int i = 0; i < 3; i++) {
        initial_state.cell_offset[i] = 0;
        initial_state.displacement[i] = 0;
        initial_state.is_flipped[i] = false;
    }
    
    DPNode* source_node = new DPNode(initial_state);
    
    // Initialize benefit values
    for (int i = 0; i < 5; i++) {
        source_node->benefit[i] = (i == 0) ? 0 : -1000000; // Only case 0 (NO_STAPLE) is valid initially
        source_node->vdd_staples[i] = 0;
        source_node->vss_staples[i] = 0;
    }
    
    return source_node;
}

/**
 * @brief Process all nodes at a given site
 */
void Optimizer::processSite(int site, std::queue<DPNode*>& queue,
                          std::unordered_map<CompactState, DPNode*, CompactStateHasher>& node_lookup_table,
                          const std::vector< std::vector<Cell*> >& cells_in_rows,
                          const std::vector<Staple>& prev_staples,
                          int row_start) {
    // Process all nodes at the current site
    size_t nodes_to_process = queue.size();
    
    for (size_t i = 0; i < nodes_to_process; i++) {
        DPNode* node = queue.front();
        queue.pop();
        
        // Skip if node is not at the current site
        if (node->state.site != site) {
            continue;
        }
        
        // Generate all valid extensions
        std::vector<Extension> extensions = generateExtensions(node, site, cells_in_rows);
        
        // Process each extension
        for (const Extension& ext : extensions) {
            // Create the target state
            CompactState target_state = node->state;
            target_state.site = site + 1;
            
            // Update state based on extension
            for (int row = 0; row < 3; row++) {
                if (ext.new_cell_idx[row] >= 0) {
                    // Update with new cell
                    Cell* cell = cells_in_rows[row][ext.new_cell_idx[row]];
                    int initial_site = cell->getInitialSite(chip_info.site_width);
                    int current_site = site + 1;
                    int cell_offset = ext.new_cell_idx[row];
                    int displacement = (current_site - initial_site);
                    
                    target_state.cell_offset[row] = cell_offset;
                    target_state.displacement[row] = displacement;
                    target_state.is_flipped[row] = ext.new_is_flipped[row];
                } else {
                    // No new cell, update distance
                    target_state.displacement[row] = ext.new_displacement[row];
                }
            }
            
            // Get or create target node
            DPNode* target_node = getOrCreateNode(target_state, queue, node_lookup_table);
            
            // Update benefit
            updateBenefit(node, target_node, site, ext, prev_staples, row_start, cells_in_rows);
        }
    }
}

/**
 * @brief Generate all valid extensions for a node
 */
std::vector<Extension> Optimizer::generateExtensions(DPNode* node, int site,
                                                 const std::vector< std::vector<Cell*> >& cells_in_rows) {
    std::vector<Extension> extensions;
    
    // Generate all combinations of extensions for the three rows
    // For simplicity, we'll just consider 2^3 = 8 possibilities:
    // Each row can either stay as is or place a new cell
    
    for (int r1_ext = 0; r1_ext < 2; r1_ext++) {
        for (int r2_ext = 0; r2_ext < 2; r2_ext++) {
            for (int r3_ext = 0; r3_ext < 2; r3_ext++) {
                // Create extension
                Extension ext(site + 1, -1);
                
                // Process row 1
                int s1 = node->state.cell_offset[0];
                int l1 = node->state.displacement[0];
                bool f1 = node->state.is_flipped[0];
                
                if (r1_ext == 0) {
                    // No new cell, increase distance
                    ext.new_cell_idx[0] = s1;
                    ext.new_displacement[0] = l1 + 1;
                    ext.new_is_flipped[0] = f1;
                } else if (s1 + 1 < static_cast<int>(cells_in_rows[0].size())) {
                    // Place next cell
                    Cell* cell = cells_in_rows[0][s1 + 1];
                    int initial_site = cell->getInitialSite(chip_info.site_width);
                    int current_site = site + 1;
                    
                    // Check if placement is within displacement limit
                    if (std::abs(current_site - initial_site) <= cell->max_displacement) {
                        // Try with and without flipping
                        for (int flip = 0; flip < 2; flip++) {
                            bool is_flipped = (flip == 1);
                            
                            // Check if the placement is valid
                            if (isValidPlacement(cell, (site + 1) * chip_info.site_width, 
                                               cell->initial_y, is_flipped)) {
                                Extension flipped_ext = ext;
                                flipped_ext.new_cell_idx[0] = s1 + 1;
                                flipped_ext.new_displacement[0] = 0;
                                flipped_ext.new_is_flipped[0] = is_flipped;
                                
                                // Process the remaining rows
                                // (similar logic for rows 2 and 3)
                                
                                extensions.push_back(flipped_ext);
                            }
                        }
                    }
                }
                
                // Similar processing for rows 2 and 3
                // (code omitted for brevity)
            }
        }
    }
    
    return extensions;
}

/**
 * @brief Check if a staple can be inserted at a given position
 */
bool Optimizer::canInsertStaple(int site, int row, 
                              const std::vector< std::vector<Cell*> >& cells_in_rows,
                              DPNode* node) {
    // Check if the two adjacent rows have empty space at the given site
    
    // First row to check
    int row1 = row;
    int s1 = node->state.cell_offset[row1];
    int l1 = node->state.displacement[row1];
    bool f1 = node->state.is_flipped[row1];
    
    // Check if there's a cell crossing the site in row1
    if (s1 < static_cast<int>(cells_in_rows[row1].size())) {
        Cell* cell1 = cells_in_rows[row1][s1];
        const CellType& type1 = cell_types[cell1->type_index];
        int site_width = type1.cellSiteWidth;
        
        if (l1 < site_width) {
            // Cell crosses the site, check if there's a pin
            int relative_site = site - (cell1->current_x / chip_info.site_width);
            if (type1.hasPinAt(relative_site, f1)) {
                return false;
            }
        }
    }
    
    // Second row to check
    int row2 = row + 1;
    int s2 = node->state.cell_offset[row2];
    int l2 = node->state.displacement[row2];
    bool f2 = node->state.is_flipped[row2];
    
    // Check if there's a cell crossing the site in row2
    if (s2 < static_cast<int>(cells_in_rows[row2].size())) {
        Cell* cell2 = cells_in_rows[row2][s2];
        const CellType& type2 = cell_types[cell2->type_index];
        int site_width = type2.cellSiteWidth;
        
        if (l2 < site_width) {
            // Cell crosses the site, check if there's a pin
            int relative_site = site - (cell2->current_x / chip_info.site_width);
            if (type2.hasPinAt(relative_site, f2)) {
                return false;
            }
        }
    }
    
    return true;
}

/**
 * @brief Check if there's a staggering violation at the given site
 */
bool Optimizer::hasStaggeringViolation(int site, int staple_case, int prev_case,
                                    const std::vector<Staple>& prev_staples,
                                    int row_start) {
    // Check for staggering between current staples
    if (staple_case == R1_R2_STAPLE && prev_case == R2_R3_STAPLE) {
        return true;  // Staggering between R1-R2 and R2-R3
    }
    if (staple_case == R2_R3_STAPLE && prev_case == R1_R2_STAPLE) {
        return true;  // Staggering between R2-R3 and R1-R2
    }
    
    // Check for staggering with previous round staples
    int site_x = site * chip_info.site_width;
    int next_site_x = (site + 1) * chip_info.site_width;
    
    for (const Staple& prev_staple : prev_staples) {
        if (prev_staple.x == site_x) {
            // Previous staple at current site
            int prev_row = prev_staple.y / chip_info.row_height;
            
            if ((staple_case == R1_R2_STAPLE && prev_row == row_start - 1) ||
                (staple_case == R2_R3_STAPLE && prev_row == row_start + 3)) {
                return true;
            }
        }
        else if (prev_staple.x == next_site_x) {
            // Previous staple at next site
            int prev_row = prev_staple.y / chip_info.row_height;
            
            if ((staple_case == R1_R2_STAPLE && prev_row == row_start + 3) ||
                (staple_case == R2_R3_STAPLE && prev_row == row_start - 1)) {
                return true;
            }
        }
    }
    
    return false;
}

/**
 * @brief Create or get a node with the given configuration
 */
DPNode* Optimizer::getOrCreateNode(const CompactState& state,
                                std::queue<DPNode*>& queue,
                                std::unordered_map<CompactState, DPNode*, CompactStateHasher>& node_lookup_table) {
    auto it = node_lookup_table.find(state);
    if (it != node_lookup_table.end()) {
        return it->second;
    }
    
    // Create new node
    DPNode* new_node = new DPNode(state);
    
    // Add to lookup table and queue
    node_lookup_table[state] = new_node;
    queue.push(new_node);
    
    return new_node;
}

/**
 * @brief Update benefit of a node based on an extension
 */
void Optimizer::updateBenefit(DPNode* from_node, DPNode* to_node, int site,
                           const Extension& extension,
                           const std::vector<Staple>& prev_staples,
                           int row_start, const std::vector< std::vector<Cell*> >& cells_in_rows) {
    // Try all combinations of staple insertion cases
    for (int from_case = 0; from_case < 5; from_case++) {
        // Skip invalid source cases
        if (from_node->benefit[from_case] < 0) {
            continue;
        }
        
        for (int to_case = 0; to_case < 5; to_case++) {
            // Calculate staple benefits based on the case
            int staple_benefit = 0;
            int vdd_staples = 0;
            int vss_staples = 0;
            
            switch (to_case) {
                case NO_STAPLE:
                    // No staples inserted
                    break;
                
                case R1_R2_STAPLE:
                    // Check if a staple can be inserted between R1 and R2
                    if (canInsertStaple(site, 0, cells_in_rows, to_node)) {
                        staple_benefit = 1;
                        if (isVDDRow(row_start)) {
                            vdd_staples = 1;
                        } else {
                            vss_staples = 1;
                        }
                    } else {
                        continue;  // Cannot insert staple, skip this case
                    }
                    break;
                
                case R2_R3_STAPLE:
                    // Check if a staple can be inserted between R2 and R3
                    if (canInsertStaple(site, 1, cells_in_rows, to_node)) {
                        staple_benefit = 1;
                        if (isVDDRow(row_start + 1)) {
                            vdd_staples = 1;
                        } else {
                            vss_staples = 1;
                        }
                    } else {
                        continue;  // Cannot insert staple, skip this case
                    }
                    break;
                
                case BOTH_STAPLES:
                    // Check if staples can be inserted in both positions
                    if (canInsertStaple(site, 0, cells_in_rows, to_node) &&
                        canInsertStaple(site, 1, cells_in_rows, to_node)) {
                        staple_benefit = 2;
                        if (isVDDRow(row_start)) {
                            vdd_staples++;
                        } else {
                            vss_staples++;
                        }
                        if (isVDDRow(row_start + 1)) {
                            vdd_staples++;
                        } else {
                            vss_staples++;
                        }
                    } else {
                        continue;  // Cannot insert both staples, skip this case
                    }
                    break;
            }
            
            // Check for staggering violations
            if (hasStaggeringViolation(site, to_case, from_case, prev_staples, row_start)) {
                continue;  // Staggering violation, skip this case
            }
            
            // Calculate new benefit
            int new_benefit = from_node->benefit[from_case] + staple_benefit;
            int new_vdd = from_node->vdd_staples[from_case] + vdd_staples;
            int new_vss = from_node->vss_staples[from_case] + vss_staples;
            
            // Check balance constraint
            if (new_vdd > 0 && new_vss > 0) {
                double ratio = static_cast<double>(std::max(new_vdd, new_vss)) / 
                              std::min(new_vdd, new_vss);
                
                if (ratio > 1.1) {
                    // Apply balance encouragement
                    if (new_vdd > new_vss && to_case == R1_R2_STAPLE && !isVDDRow(row_start)) {
                        new_benefit += params.balance_factor;
                    }
                    else if (new_vdd > new_vss && to_case == R2_R3_STAPLE && !isVDDRow(row_start + 1)) {
                        new_benefit += params.balance_factor;
                    }
                    else if (new_vss > new_vdd && to_case == R1_R2_STAPLE && isVDDRow(row_start)) {
                        new_benefit += params.balance_factor;
                    }
                    else if (new_vss > new_vdd && to_case == R2_R3_STAPLE && isVDDRow(row_start + 1)) {
                        new_benefit += params.balance_factor;
                    }
                }
            }
            
            // Update benefit if better
            if (new_benefit > to_node->benefit[to_case]) {
                to_node->benefit[to_case] = new_benefit;
                to_node->prev_node[to_case] = from_node;
                to_node->case_from_prev[to_case] = from_case;
                to_node->vdd_staples[to_case] = new_vdd;
                to_node->vss_staples[to_case] = new_vss;
            }
        }
    }
}

/**
 * @brief Backtrack through the DAG to get the optimal solution
 */
std::vector<Staple> Optimizer::backtrack(DPNode* best_node,
                                      const std::vector< std::vector<Cell*> >& cells_in_rows,
                                      int row_start) {
    std::vector<Staple> staples;
    
    // Find the case with maximum benefit
    int best_case = 0;
    int max_benefit = best_node->benefit[0];
    
    for (int i = 1; i < 5; i++) {
        if (best_node->benefit[i] > max_benefit) {
            max_benefit = best_node->benefit[i];
            best_case = i;
        }
    }
    
    // Backtrack from best node
    DPNode* current = best_node;
    int current_case = best_case;
    int site = chip_info.total_sites;
    
    while (current != nullptr && current->prev_node[current_case] != nullptr) {
        int prev_case = current->case_from_prev[current_case];
        DPNode* prev_node = current->prev_node[current_case];
        
        // Insert staples based on the case
        switch (current_case) {
            case R1_R2_STAPLE: {
                // Insert staple between R1 and R2
                int x = site * chip_info.site_width;
                int y = (row_start + 1) * chip_info.row_height;
                bool is_vdd = isVDDRow(row_start);
                staples.push_back(Staple(x, y, is_vdd));
                break;
            }
            
            case R2_R3_STAPLE: {
                // Insert staple between R2 and R3
                int x = site * chip_info.site_width;
                int y = (row_start + 2) * chip_info.row_height;
                bool is_vdd = isVDDRow(row_start + 1);
                staples.push_back(Staple(x, y, is_vdd));
                break;
            }
            
            case BOTH_STAPLES: {
                // Insert staples between R1-R2 and R2-R3
                int x = site * chip_info.site_width;
                
                int y1 = (row_start + 1) * chip_info.row_height;
                bool is_vdd1 = isVDDRow(row_start);
                staples.push_back(Staple(x, y1, is_vdd1));
                
                int y2 = (row_start + 2) * chip_info.row_height;
                bool is_vdd2 = isVDDRow(row_start + 1);
                staples.push_back(Staple(x, y2, is_vdd2));
                break;
            }

        }
        
        // Move to previous node
        site--;
        current = prev_node;
        current_case = prev_case;
    }
    
    // Update cell positions based on the final state
    // (This would need to be implemented based on the backtracking path)
    
    return staples;
}

/**
 * @brief Check if a cell placement is valid
 */
bool Optimizer::isValidPlacement(const Cell* cell, int x, int y, bool is_flipped) {
    // Check if the placement is within displacement limit
    int displacement = std::abs(x - cell->initial_x);
    if (displacement > cell->max_displacement * chip_info.site_width) {
        return false;
    }
    
    // Check if the placement is in the same row
    if (y != cell->initial_y) {
        return false;
    }
    
    // Other checks can be added as needed
    
    return true;
}

/**
 * @brief Get the compact state encoding for a configuration
 */
CompactState Optimizer::getCompactState(int site, 
                                     int s1, int l1, bool f1,
                                     int s2, int l2, bool f2,
                                     int s3, int l3, bool f3,
                                     const std::vector< std::vector<Cell*> >& cells_in_rows) {
    CompactState state;
    state.site = site;
    
    // Calculate compact encoding for row 1
    if (cells_in_rows[0].empty()) {
        state.cell_offset[0] = 0;
        state.displacement[0] = 0;
        state.is_flipped[0] = false;
    } else {
        Cell* cell1 = cells_in_rows[0][s1];
        int initial_site1 = cell1->getInitialSite(chip_info.site_width);
        int initial_s1 = -1;
        
        // Find the cell that would be here in the initial placement
        for (size_t i = 0; i < cells_in_rows[0].size(); i++) {
            Cell* c = cells_in_rows[0][i];
            int c_site = c->getInitialSite(chip_info.site_width);
            if (c_site <= site && site < c_site + cell_types[c->type_index].cellSiteWidth) {
                initial_s1 = i;
                break;
            }
        }
        
        if (initial_s1 == -1) {
            initial_s1 = cells_in_rows[0].size() - 1;
        }
        
        state.cell_offset[0] = s1 - initial_s1;
        state.displacement[0] = l1;
        state.is_flipped[0] = f1;
    }
    
    // Similar calculations for rows 2 and 3
    // (code omitted for brevity)
    
    return state;
}

/**
 * @brief Determine if a position is a VDD or VSS row
 */
bool Optimizer::isVDDRow(int row_idx) const {
    // Alternate rows for power and ground
    // Even rows (0, 2, 4, ...) are VDD, odd rows are VSS
    return (row_idx % 2 == 0);
}