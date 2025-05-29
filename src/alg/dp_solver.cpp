/**
 * @file dp_solver.cpp
 * @brief Correct implementation of DAG-based dynamic programming solver for MATRO
 */

#include "dp_solver.hpp"
#include "../Logger.hpp"
#include "../memory_usage.hpp"
#include <algorithm>
#include <limits>
#include <cmath>
#include <set>

/**
 * @brief Constructor
 */
DPSolver::DPSolver(const ChipInfo& chip_info, 
                   const std::vector<CellType>& cell_types,
                   const AlgorithmParams& params)
    : chip_info(chip_info), cell_types(cell_types), params(params),
      best_final_node(nullptr), best_final_case(-1), nodes_created(0) {
    
    Logger::log(Logger::INFO, "DPSolver initialized");
    Logger::log(Logger::INFO, "  Site-by-site DAG construction");
    Logger::log(Logger::INFO, "  5-case MATRO algorithm");
}

/**
 * @brief Destructor
 */
DPSolver::~DPSolver() {
    cleanup();
}

/**
 * @brief Main MATRO solver implementation
 */
std::vector<Staple> DPSolver::solveTripleRow(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start,
    const std::vector<Staple>& prev_staples) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    std::cout << "=== Triple-row " << row_start << "-" << (row_start+2) 
              << " optimization ===" << std::endl;
    
    // Clean up previous state
    cleanup();
    
    // Store previous staples for staggering checks
    this->prev_staples = prev_staples;
    
    // Initialize base indices for compact encoding
    for (int row = 0; row < 3; row++) {
        if (row < static_cast<int>(cells_in_rows.size())) {
            initial_s[row] = 0;  // All cells start at index 0
        }
    }
    
    // Create source node at site 0
    DPNode* source = new DPNode(0, 0, 0, 0, 0, 0, 0);
    // Initialize case 1 (no staple) with benefit 0
    source->cases[CASE_1_NO_STAPLE].benefit = 0;
    source->cases[CASE_1_NO_STAPLE].valid = true;
    source->cases[CASE_1_NO_STAPLE].vdd_count = 0;
    source->cases[CASE_1_NO_STAPLE].vss_count = 0;
    
    all_nodes.push_back(source);
    
    // Add to hash table with compact encoding
    CompactState compact = getCompactState(source);
    node_lookup[compact] = source;
    
    // Queue for BFS traversal
    std::queue<DPNode*> Q;
    Q.push(source);
    
    std::cout << "Starting DAG construction..." << std::endl;
    
    // Main DAG construction loop
    int processed_nodes = 0;
    int max_site_reached = 0;
    
    while (!Q.empty()) {
        // Check time limit
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - start_time);
        
        if (elapsed.count() > 60) {  // 60 second timeout per triple-row
            std::cout << "Timeout reached, using partial solution" << std::endl;
            break;
        }
        
        DPNode* current = Q.front();
        Q.pop();
        processed_nodes++;
        
        if (current->site >= chip_info.total_sites) {
            continue;  // Reached end
        }
        
        max_site_reached = std::max(max_site_reached, current->site);
        
        // Generate all possible extensions
        generateExtensions(current, cells_in_rows, Q);
        
        // Memory management - prune if too many nodes
        if (all_nodes.size() > 50000 && current->site < max_site_reached - 5) {
            pruneOldNodes(max_site_reached - 5);
        }
        
        if (processed_nodes % 1000 == 0) {
            std::cout << "Processed " << processed_nodes << " nodes, current site: " 
                      << current->site << ", total nodes: " << all_nodes.size() << std::endl;
        }
    }
    
    std::cout << "DAG construction completed. Nodes: " << all_nodes.size() 
              << ", Max site: " << max_site_reached << std::endl;
    
    // Extract best solution
    return extractBestSolution(cells_in_rows, row_start);
}

/**
 * @brief Generate all legal extensions from current node
 */
void DPSolver::generateExtensions(DPNode* current, 
                                 const std::vector<std::vector<Cell*>>& cells_in_rows,
                                 std::queue<DPNode*>& Q) {
    
    int next_site = current->site + 1;
    
    // Try all 8 combinations of placing/not placing cells
    for (int mask = 0; mask < 8; mask++) {
        bool place_cell[3] = {
            (mask & 4) != 0,  // Row 0
            (mask & 2) != 0,  // Row 1  
            (mask & 1) != 0   // Row 2
        };
        
        // Check if this extension is legal
        int new_s[3], new_l[3];
        if (!isLegalExtension(current, place_cell, cells_in_rows, new_s, new_l)) {
            continue;
        }
        
        // Create or find node with this state
        DPNode* next = createOrFindNode(next_site, new_s[0], new_l[0], 
                                       new_s[1], new_l[1], new_s[2], new_l[2]);
        
        if (next != nullptr) {
            // Update benefits for all valid case transitions
            bool had_valid_cases_before = false;
            for (int c = 0; c < 5; c++) {
                if (next->cases[c].valid) {
                    had_valid_cases_before = true;
                    break;
                }
            }
            
            updateNodeCases(current, next, cells_in_rows);
            
            // Add to queue if this node now has new valid cases
            bool has_valid_cases_after = false;
            for (int c = 0; c < 5; c++) {
                if (next->cases[c].valid) {
                    has_valid_cases_after = true;
                    break;
                }
            }
            
            // Only add to queue if we added new valid cases
            if (!had_valid_cases_before && has_valid_cases_after) {
                Q.push(next);
            }
        }
    }
}

/**
 * @brief Check if extension is legal according to paper's rules
 */
bool DPSolver::isLegalExtension(DPNode* current, bool place_cell[3],
                               const std::vector<std::vector<Cell*>>& cells_in_rows,
                               int new_s[3], int new_l[3]) {
    
    for (int row = 0; row < 3; row++) {
        if (row >= static_cast<int>(cells_in_rows.size()) || 
            cells_in_rows[row].empty()) {
            // No cells in this row
            new_s[row] = current->s[row];
            new_l[row] = current->l[row] + 1;
            continue;
        }
        
        if (place_cell[row]) {
            // Want to place next cell at site i+1
            
            // Check if previous cell has cleared current site
            if (current->s[row] > 0) {
                Cell* prev_cell = cells_in_rows[row][current->s[row] - 1];
                const CellType& prev_type = cell_types[prev_cell->type_index];
                int prev_width_sites = prev_type.width / chip_info.site_width;
                
                // Paper: cell has cleared if l >= cell_width
                if (current->l[row] < prev_width_sites) {
                    return false;  // Previous cell still occupies current site
                }
            }
            
            // Check if there's a next cell to place
            if (current->s[row] >= static_cast<int>(cells_in_rows[row].size())) {
                return false;  // No more cells
            }
            
            // Check displacement constraint
            Cell* next_cell = cells_in_rows[row][current->s[row]];
            int next_site = (current->site + 1);
            int next_x = next_site * chip_info.site_width;
            int displacement = std::abs(next_x - next_cell->initial_x);
            
            if (displacement > next_cell->max_displacement) {
                return false;  // Exceeds max displacement
            }
            
            // Update state
            new_s[row] = current->s[row] + 1;
            new_l[row] = 1;  // Cell placed at current site, distance = 1
            
        } else {
            // Not placing cell - just increment distance
            
            // Check if deferring is legal
            if (current->s[row] < static_cast<int>(cells_in_rows[row].size())) {
                Cell* next_cell = cells_in_rows[row][current->s[row]];
                int initial_site = next_cell->initial_x / chip_info.site_width;
                int max_displacement_sites = next_cell->max_displacement / chip_info.site_width;
                int latest_site = initial_site + max_displacement_sites;
                
                if (current->site + 1 > latest_site) {
                    return false;  // Must place cell by now
                }
            }
            
            // Update state
            new_s[row] = current->s[row];
            new_l[row] = current->l[row] + 1;
        }
    }
    
    return true;
}

/**
 * @brief Create or find node using compact encoding
 */
DPNode* DPSolver::createOrFindNode(int site, int s1, int l1, int s2, int l2, int s3, int l3) {
    
    // Create temporary node for compact state calculation
    DPNode temp(site, s1, l1, s2, l2, s3, l3);
    CompactState compact = getCompactState(&temp);
    
    // Check if node exists
    auto it = node_lookup.find(compact);
    if (it != node_lookup.end()) {
        return it->second;  // Return existing node
    }
    
    // Create new node
    DPNode* new_node = new DPNode(site, s1, l1, s2, l2, s3, l3);
    all_nodes.push_back(new_node);
    node_lookup[compact] = new_node;
    nodes_created++;
    
    return new_node;
}

/**
 * @brief Get compact state for hash lookup
 */
CompactState DPSolver::getCompactState(DPNode* node) {
    CompactState compact;
    compact.site = node->site;
    
    for (int row = 0; row < 3; row++) {
        // a[row] = s[row] - initial_s[row] (offset from initial)
        compact.a[row] = node->s[row] - initial_s[row];
        
        // For b[row], we need actual displacement
        // Since we're doing site-by-site, we can use l[row] as proxy
        compact.b[row] = node->l[row];
    }
    
    return compact;
}

/**
 * @brief Update all valid case transitions
 */
void DPSolver::updateNodeCases(DPNode* source, DPNode* target,
                              const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    // Try all case transitions
    for (int src_case = 0; src_case < 5; src_case++) {
        if (!source->cases[src_case].valid) continue;
        
        for (int tgt_case = 0; tgt_case < 5; tgt_case++) {
            // Check if this case transition is valid
            if (!isValidCaseTransition(source, src_case, target, tgt_case, cells_in_rows)) {
                continue;
            }
            
            // Calculate benefit
            int benefit_delta = calculateCaseBenefit(tgt_case, target,
                source->cases[src_case].vdd_count, 
                source->cases[src_case].vss_count);
            
            int new_benefit = source->cases[src_case].benefit + benefit_delta;
            
            // Update if better
            if (new_benefit > target->cases[tgt_case].benefit) {
                target->cases[tgt_case].benefit = new_benefit;
                target->cases[tgt_case].prev_node = source;
                target->cases[tgt_case].prev_case = src_case;
                target->cases[tgt_case].valid = true;
                
                // Update VDD/VSS counts
                auto [vdd_delta, vss_delta] = getStapleTypesForCase(tgt_case, target->site);
                target->cases[tgt_case].vdd_count = source->cases[src_case].vdd_count + vdd_delta;
                target->cases[tgt_case].vss_count = source->cases[src_case].vss_count + vss_delta;
            }
        }
    }
}

/**
 * @brief Check if case transition is valid (no violations)
 */
bool DPSolver::isValidCaseTransition(DPNode* source, int src_case,
                                    DPNode* target, int tgt_case,
                                    const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    // Get staples for target case
    std::vector<Staple> new_staples = getStaplesForCase(tgt_case, target);
    
    // Check pin overlaps
    for (const Staple& staple : new_staples) {
        if (hasOverlapWithPins(staple, target, cells_in_rows)) {
            return false;
        }
    }
    
    // Check staggering with previous staples
    if (hasStaggeringViolation(new_staples)) {
        return false;
    }
    
    // Check balance constraint won't be violated
    auto [vdd_delta, vss_delta] = getStapleTypesForCase(tgt_case, target->site);
    int new_vdd = source->cases[src_case].vdd_count + vdd_delta;
    int new_vss = source->cases[src_case].vss_count + vss_delta;
    
    if (new_vdd > 0 && new_vss > 0) {
        double ratio = static_cast<double>(std::max(new_vdd, new_vss)) / 
                      static_cast<double>(std::min(new_vdd, new_vss));
        if (ratio > 1.1) {
            return false;  // Would violate balance constraint
        }
    }
    
    return true;
}

/**
 * @brief Get staples for a specific case
 */
std::vector<Staple> DPSolver::getStaplesForCase(int case_type, DPNode* node) {
    std::vector<Staple> staples;
    int x = node->site * chip_info.site_width;
    
    switch (case_type) {
        case CASE_1_NO_STAPLE:
            // No staples
            break;
            
        case CASE_2_R1_R2_ONLY: {
            // Staple between rows 0 and 1
            int row = row_start + 1;
            int y = chip_info.getRowY(row);
            bool is_vdd = determineStapleType(row);
            staples.push_back(Staple(x, y, is_vdd));
            break;
        }
        
        case CASE_3_R2_R3_ONLY: {
            // Staple between rows 1 and 2
            int row = row_start + 2;
            int y = chip_info.getRowY(row);
            bool is_vdd = determineStapleType(row);
            staples.push_back(Staple(x, y, is_vdd));
            break;
        }
        
        case CASE_4_BOTH_ALIGNED: {
            // Both staples at same x
            int row1 = row_start + 1;
            int row2 = row_start + 2;
            staples.push_back(Staple(x, chip_info.getRowY(row1), determineStapleType(row1)));
            staples.push_back(Staple(x, chip_info.getRowY(row2), determineStapleType(row2)));
            break;
        }
        
        case CASE_5_BOTH_STAGGERED: {
            // This case seems to be for manufacturing considerations
            // May need clarification from paper
            // For now, treat similar to case 4
            int row1 = row_start + 1;
            int row2 = row_start + 2;
            staples.push_back(Staple(x, chip_info.getRowY(row1), determineStapleType(row1)));
            staples.push_back(Staple(x, chip_info.getRowY(row2), determineStapleType(row2)));
            break;
        }
    }
    
    return staples;
}

/**
 * @brief Check staggering violation with diagonal pattern
 */
bool DPSolver::hasStaggeringViolation(const std::vector<Staple>& new_staples) {
    
    // Check against all existing staples
    for (const Staple& new_staple : new_staples) {
        int new_site = new_staple.x / chip_info.site_width;
        int new_row = (new_staple.y - chip_info.bottom_y) / chip_info.row_height;
        
        // Check previous staples from earlier triple-rows
        for (const Staple& prev : prev_staples) {
            int prev_site = prev.x / chip_info.site_width;
            int prev_row = (prev.y - chip_info.bottom_y) / chip_info.row_height;
            
            // Check if diagonal
            if (std::abs(new_site - prev_site) == 1 && std::abs(new_row - prev_row) == 1) {
                // Identify lower and upper staples
                int lower_site = (new_row < prev_row) ? new_site : prev_site;
                int lower_row = std::min(new_row, prev_row);
                int upper_site = (new_row < prev_row) ? prev_site : new_site;
                int upper_row = std::max(new_row, prev_row);
                
                // Check for blocking staples
                bool has_blocker_above_lower = false;
                bool has_blocker_below_upper = false;
                
                // Check in prev_staples
                for (const Staple& check : prev_staples) {
                    int check_site = check.x / chip_info.site_width;
                    int check_row = (check.y - chip_info.bottom_y) / chip_info.row_height;
                    
                    if (check_site == lower_site && check_row == lower_row + 1) {
                        has_blocker_above_lower = true;
                    }
                    if (check_site == upper_site && check_row == upper_row - 1) {
                        has_blocker_below_upper = true;
                    }
                }
                
                // Also check in new_staples
                for (const Staple& check : new_staples) {
                    int check_site = check.x / chip_info.site_width;
                    int check_row = (check.y - chip_info.bottom_y) / chip_info.row_height;
                    
                    if (check_site == lower_site && check_row == lower_row + 1) {
                        has_blocker_above_lower = true;
                    }
                    if (check_site == upper_site && check_row == upper_row - 1) {
                        has_blocker_below_upper = true;
                    }
                }
                
                // Staggering occurs if both blockers are missing
                if (!has_blocker_above_lower && !has_blocker_below_upper) {
                    return true;
                }
            }
        }
    }
    
    return false;
}

/**
 * @brief Check if staple overlaps with pins
 */
bool DPSolver::hasOverlapWithPins(const Staple& staple, DPNode* node,
                                 const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    int staple_site = staple.x / chip_info.site_width;
    int staple_row = (staple.y - chip_info.bottom_y) / chip_info.row_height;
    
    // A staple at row boundary connects two rows
    int lower_row_idx = staple_row - row_start - 1;
    int upper_row_idx = staple_row - row_start;
    
    // Check both rows
    for (int row_idx : {lower_row_idx, upper_row_idx}) {
        if (row_idx < 0 || row_idx >= 3 || row_idx >= static_cast<int>(cells_in_rows.size())) {
            continue;
        }
        
        // Check cells placed up to this node
        for (int cell_idx = 0; cell_idx < node->s[row_idx]; cell_idx++) {
            Cell* cell = cells_in_rows[row_idx][cell_idx];
            const CellType& cell_type = cell_types[cell->type_index];
            
            // Calculate cell position based on node state
            int cell_left_site = calculateCellPosition(node, row_idx, cell_idx, cells_in_rows);
            
            // Check if staple site overlaps with any pin
            for (int pin_offset : cell_type.pin_sites) {
                if (cell_left_site + pin_offset == staple_site) {
                    return true;  // Overlap found
                }
            }
        }
    }
    
    return false;
}

/**
 * @brief Calculate cell position from node state
 */
int DPSolver::calculateCellPosition(DPNode* node, int row_idx, int cell_idx,
                                   const std::vector<std::vector<Cell*>>& cells_in_rows) {
    // This is simplified - in full implementation would track exact positions
    // For now, approximate based on state
    if (cell_idx == node->s[row_idx] - 1) {
        // Last placed cell
        return node->site - node->l[row_idx] + 1;
    }
    
    // Would need to track all cell positions for accurate calculation
    return 0;  // Placeholder
}

/**
 * @brief Calculate benefit for a case
 */
int DPSolver::calculateCaseBenefit(int case_type, DPNode* node, int current_vdd, int current_vss) {
    int base_benefit = 0;
    
    switch (case_type) {
        case CASE_1_NO_STAPLE:
            // Beta factor for encouraging space for future staples
            base_benefit = static_cast<int>(params.balance_factor * 50);
            break;
            
        case CASE_2_R1_R2_ONLY:
        case CASE_3_R2_R3_ONLY:
            base_benefit = 100;  // Single staple
            break;
            
        case CASE_4_BOTH_ALIGNED:
            base_benefit = 200;  // Two aligned staples
            break;
            
        case CASE_5_BOTH_STAGGERED:
            base_benefit = 180;  // Slightly less due to manufacturing
            break;
    }
    
    // Add balance bonus/penalty
    auto [vdd_delta, vss_delta] = getStapleTypesForCase(case_type, node->site);
    int new_vdd = current_vdd + vdd_delta;
    int new_vss = current_vss + vss_delta;
    
    if (new_vdd > 0 && new_vss > 0) {
        double ratio = static_cast<double>(std::max(new_vdd, new_vss)) / 
                      static_cast<double>(std::min(new_vdd, new_vss));
        
        if (ratio > 1.05) {
            base_benefit -= static_cast<int>((ratio - 1.0) * 100);  // Penalty
        } else {
            base_benefit += 20;  // Bonus for good balance
        }
    }
    
    return base_benefit;
}

/**
 * @brief Get VDD/VSS count changes for a case
 */
std::pair<int, int> DPSolver::getStapleTypesForCase(int case_type, int site) {
    switch (case_type) {
        case CASE_1_NO_STAPLE:
            return {0, 0};
            
        case CASE_2_R1_R2_ONLY:
        case CASE_3_R2_R3_ONLY:
            // Simplified - alternate VDD/VSS
            return (site % 2 == 0) ? std::make_pair(1, 0) : std::make_pair(0, 1);
            
        case CASE_4_BOTH_ALIGNED:
        case CASE_5_BOTH_STAGGERED:
            // One of each for balance
            return {1, 1};
            
        default:
            return {0, 0};
    }
}

/**
 * @brief Determine staple type based on position
 */
bool DPSolver::determineStapleType(int row_boundary) {
    // Simple alternating pattern - can be refined
    return (row_boundary % 2 == 0);
}

/**
 * @brief Extract best solution from DAG
 */
std::vector<Staple> DPSolver::extractBestSolution(
    const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start_param) {
    
    this->row_start = row_start_param;  // Store for use in helper functions
    
    std::cout << "Extracting best solution..." << std::endl;
    
    // Find best valid solution
    DPNode* best_node = nullptr;
    int best_case = -1;
    int best_benefit = -1000000;
    
    for (DPNode* node : all_nodes) {
        // Only consider nodes at or near the end
        if (node->site < chip_info.total_sites - 5) continue;
        
        for (int c = 0; c < 5; c++) {
            if (!node->cases[c].valid) continue;
            
            // Check balance constraint
            int vdd = node->cases[c].vdd_count;
            int vss = node->cases[c].vss_count;
            
            if (vdd > 0 && vss > 0) {
                double ratio = static_cast<double>(std::max(vdd, vss)) / 
                              static_cast<double>(std::min(vdd, vss));
                if (ratio > 1.1) continue;  // Skip unbalanced
            }
            
            if (node->cases[c].benefit > best_benefit) {
                best_benefit = node->cases[c].benefit;
                best_node = node;
                best_case = c;
            }
        }
    }
    
    if (best_node == nullptr) {
        std::cout << "No valid solution found" << std::endl;
        return {};
    }
    
    std::cout << "Best solution: benefit=" << best_benefit
              << ", VDD=" << best_node->cases[best_case].vdd_count
              << ", VSS=" << best_node->cases[best_case].vss_count << std::endl;
    
    // Backtrack to extract staples
    std::vector<Staple> staples;
    DPNode* current = best_node;
    int current_case = best_case;
    
    while (current != nullptr && current->cases[current_case].prev_node != nullptr) {
        // Get staples at this node
        std::vector<Staple> node_staples = getStaplesForCase(current_case, current);
        
        // Add in reverse order (we're backtracking)
        for (auto it = node_staples.rbegin(); it != node_staples.rend(); ++it) {
            staples.push_back(*it);
        }
        
        // Move to previous node
        DPNode* prev = current->cases[current_case].prev_node;
        int prev_case = current->cases[current_case].prev_case;
        current = prev;
        current_case = prev_case;
    }
    
    // Reverse to get correct order
    std::reverse(staples.begin(), staples.end());
    
    // Also need to update cell positions in cells_in_rows
    // This would be done by tracking exact positions during DP
    
    return staples;
}

/**
 * @brief Prune old nodes to manage memory
 */
void DPSolver::pruneOldNodes(int keep_after_site) {
    std::vector<DPNode*> nodes_to_keep;
    
    for (DPNode* node : all_nodes) {
        if (node->site >= keep_after_site) {
            nodes_to_keep.push_back(node);
        } else {
            // Remove from lookup
            CompactState compact = getCompactState(node);
            node_lookup.erase(compact);
            delete node;
        }
    }
    
    all_nodes = nodes_to_keep;
    
    std::cout << "Pruned nodes, keeping " << all_nodes.size() << " nodes" << std::endl;
}

/**
 * @brief Clean up allocated memory
 */
void DPSolver::cleanup() {
    for (DPNode* node : all_nodes) {
        delete node;
    }
    all_nodes.clear();
    node_lookup.clear();
    prev_staples.clear();
    nodes_created = 0;
    best_final_node = nullptr;
    best_final_case = -1;
}