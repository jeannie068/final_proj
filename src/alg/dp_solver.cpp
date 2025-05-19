/**
 * @file dp_solver.cpp
 * @brief Implementation of dynamic programming solver for power staple insertion
 */

#include "dp_solver.hpp"
#include "../Logger.hpp"

/**
 * @brief Constructor
 */
DPSolver::DPSolver(const ChipInfo& chip_info, 
                 const std::vector<CellType>& cell_types,
                 const AlgorithmParams& params)
    : chip_info(chip_info), cell_types(cell_types), params(params) {
    
    Logger::log("DPSolver initialized with parameters:");
    Logger::log("  Balance factor: " + std::to_string(params.balance_factor));
    Logger::log("  Total sites: " + std::to_string(chip_info.total_sites));
    Logger::log("  Site width: " + std::to_string(chip_info.site_width));
    Logger::log("  Row height: " + std::to_string(chip_info.row_height));
    Logger::log("  Cell types: " + std::to_string(cell_types.size()));
    
    // Initialize statistics
    max_nodes_count = 0;
    nodes_created = 0;
    nodes_processed = 0;
    best_node = nullptr;
    max_benefit = -1;
}

/**
 * @brief Destructor
 */
DPSolver::~DPSolver() {
    cleanup();
}

/**
 * @brief Solve a triple-row optimization problem
 * 
 * This is the main entry point for the dynamic programming solver.
 * It implements the DAG-based dynamic programming approach for triple-row
 * optimization as described in the paper.
 */
std::vector<Staple> DPSolver::solveTripleRow(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start,
    const std::vector<Staple>& prev_staples) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    const int time_limit_seconds = 120; // Safety limit
    
    Logger::log("====================================================");
    Logger::log(Logger::DEBUG, "Cells in row " + std::to_string(row_start) + ": " + 
                std::to_string(cells_in_rows[0].size()));
    Logger::log(Logger::DEBUG, "Cells in row " + std::to_string(row_start + 1) + ": " + 
                std::to_string(cells_in_rows[1].size()));
    Logger::log(Logger::DEBUG, "Cells in row " + std::to_string(row_start + 2) + ": " + 
                std::to_string(cells_in_rows[2].size()));
    Logger::log(Logger::DEBUG, "Previous staples: " + std::to_string(prev_staples.size()));
    
    // Initialize the DAG
    DPNode* source_node = initializeDAG();
    Logger::log(Logger::DEBUG, "DAG initialized with source node");
    
    // Clear data structures for a fresh run
    best_node = nullptr;
    max_benefit = -1;
    
    // Clear lookup table and queue
    node_lookup_table.clear();
    while (!node_queue.empty()) node_queue.pop();
    all_nodes.clear();
    
    // Add source node
    node_queue.push(source_node);
    node_lookup_table[source_node->state] = source_node;
    all_nodes.push_back(source_node);
    nodes_created = 1;
    
    Logger::log("Beginning site-by-site processing, total sites: " + 
                std::to_string(chip_info.total_sites));
    
    // Process each site from left to right
    for (int site = 0; site <= chip_info.total_sites; site++) {
        // Check time limit
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - start_time).count();
            
        if (elapsed_seconds > time_limit_seconds) {
            Logger::log(Logger::WARNING, "WARNING: Time limit (" + std::to_string(time_limit_seconds) + 
                      "s) reached, stopping processing early");
            std::cout << "WARNING: Time limit reached at site " << site << "/" 
                    << chip_info.total_sites << std::endl;
            break;
        }
        
        processSite(site, cells_in_rows, prev_staples, row_start);
        
        // Log progress periodically
        if (site % 1000 == 0 || site == chip_info.total_sites) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                current_time - start_time);
            
            Logger::log(Logger::TRACE, "Site " + std::to_string(site) + "/" + 
                        std::to_string(chip_info.total_sites) + 
                        " processed in " + std::to_string(duration.count()) + " ms, " +
                        "queue size: " + std::to_string(node_queue.size()) + 
                        ", nodes: " + std::to_string(node_lookup_table.size()) +
                        ", total created: " + std::to_string(nodes_created));
            
            // Monitor memory usage
            monitorMemoryUsage();
        }
        
        // Clear lookup table after processing each site to save memory
        // But keep the nodes at current site
        if (site < chip_info.total_sites) {
            size_t before_size = node_lookup_table.size();
            node_lookup_table.clear();
            Logger::log(Logger::DEBUG, "Lookup table cleared: " + std::to_string(before_size) + 
                        " entries removed");
        }
    }
    
    // Backtrack to get optimal solution if a best node was found
    std::vector<Staple> inserted_staples;
    if (best_node != nullptr) {
        Logger::log("Best solution found with benefit: " + std::to_string(max_benefit));
        Logger::log("VDD staples: " + std::to_string(best_node->vdd_staples[0]) +
                    ", VSS staples: " + std::to_string(best_node->vss_staples[0]));
        Logger::log("Starting backtracking to recover solution");
        inserted_staples = backtrack(cells_in_rows, row_start);
    } else {
        Logger::log("WARNING: No valid solution found!");
    }
    
    // Apply cell placement updates
    applyPlacementUpdates(cells_in_rows);
    Logger::log("Cell placement updates applied");
    
    // Final stats
    auto end_time = std::chrono::high_resolution_clock::now();
    auto total_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                end_time - start_time);
    
    double vdd_vss_ratio = 0.0;
    int vdd_count = 0, vss_count = 0;
    for (const auto& staple : inserted_staples) {
        if (staple.is_vdd) vdd_count++;
        else vss_count++;
    }
    
    if (vdd_count > 0 && vss_count > 0) {
        vdd_vss_ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                        std::min(vdd_count, vss_count);
    }
    
    Logger::log("====================================================");
    Logger::log("Triple-row optimization completed in " + 
                std::to_string(total_duration.count()) + " ms");
    Logger::log("Staples inserted: " + std::to_string(inserted_staples.size()));
    Logger::log("VDD staples: " + std::to_string(vdd_count) + 
                ", VSS staples: " + std::to_string(vss_count) + 
                ", Ratio: " + std::to_string(vdd_vss_ratio));
    Logger::log("Nodes created: " + std::to_string(nodes_created) +
                ", Nodes processed: " + std::to_string(nodes_processed));
    Logger::log("Max nodes in memory: " + std::to_string(max_nodes_count));
    Logger::log("====================================================");
    
    // Clean up memory
    cleanup();
    
    return inserted_staples;
}

/**
 * @brief Initialize the DAG with a source node
 */
DPNode* DPSolver::initializeDAG() {
    // Create initial state
    CompactState initial_state;
    initial_state.site = 0;
    
    for (int i = 0; i < 3; i++) {
        initial_state.cell_offset[i] = 0;
        initial_state.displacement[i] = 0;
        initial_state.is_flipped[i] = false;
    }
    
    // Create source node
    DPNode* source_node = new DPNode(initial_state);
    Logger::log("  Created source node at site 0");
    
    // Initialize benefit values (only NO_STAPLE case is valid initially)
    for (int i = 0; i < 4; i++) {
        source_node->benefit[i] = (i == NO_STAPLE) ? 0 : -1000000;
        source_node->vdd_staples[i] = 0;
        source_node->vss_staples[i] = 0;
    }
    
    return source_node;
}

/**
 * @brief Process all nodes at a given site
 */
void DPSolver::processSite(int site,
                        const std::vector<std::vector<Cell*>>& cells_in_rows,
                        const std::vector<Staple>& prev_staples,
                        int row_start) {
    
    auto process_start_time = std::chrono::high_resolution_clock::now();
    const int max_processing_ms = 30000; // 30 seconds max per site
    
    size_t nodes_to_process = node_queue.size();
    Logger::log(Logger::INFO, "Processing site " + std::to_string(site) + 
                ", nodes to process: " + std::to_string(nodes_to_process));
    Logger::increaseIndent();
    
    size_t extensions_count = 0;
    size_t valid_extensions = 0;
    
    for (size_t i = 0; i < nodes_to_process; i++) {
        // Check processing time limit periodically
        if (i % 1000 == 0 && i > 0) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                current_time - process_start_time).count();
                
            if (elapsed_ms > max_processing_ms) {
                Logger::log(Logger::WARNING, "Site processing taking too long (" + 
                          std::to_string(elapsed_ms/1000) + " seconds), limiting node processing");
                std::cout << "WARNING: Site " << site << " processing timeout, processed " 
                        << i << "/" << nodes_to_process << " nodes" << std::endl;
                break;
            }
        }
        
        DPNode* node = node_queue.front();
        node_queue.pop();
        nodes_processed++;
        
        // Skip if node is not at the current site
        if (node->state.site != site) {
            Logger::log(Logger::TRACE, "Node skipped, not at current site. Node site: " + 
                      std::to_string(node->state.site) + ", Current site: " + 
                      std::to_string(site));
            continue;
        }
        
        // Generate all valid extensions
        std::vector<Extension> extensions = generateExtensions(node, site, cells_in_rows, row_start);
        extensions_count += extensions.size();
        
        // Log less frequently for better performance
        if (nodes_processed % 10000 == 0 || nodes_to_process < 10) {
            Logger::log(Logger::DEBUG, "Node " + std::to_string(nodes_processed) + 
                        ": Generated " + std::to_string(extensions.size()) + 
                        " extensions");
        }
        
        // Process each extension
        for (const Extension& ext : extensions) {
            // Create the target state
            CompactState target_state = node->state;
            target_state.site = site + 1;
            
            // Update state based on extension
            for (int row = 0; row < 3; row++) {
                if (ext.new_cell_idx[row] >= 0 && ext.new_cell_idx[row] < static_cast<int>(cells_in_rows[row].size())) {
                    // Update with new cell
                    Cell* cell = cells_in_rows[row][ext.new_cell_idx[row]];
                    int initial_site = cell->getInitialSite(chip_info.site_width);
                    int current_site = site + 1;
                    
                    target_state.cell_offset[row] = ext.new_cell_idx[row];
                    target_state.displacement[row] = ext.new_displacement[row];
                    target_state.is_flipped[row] = ext.new_is_flipped[row];
                } else {
                    // No new cell, update distance
                    target_state.displacement[row] = ext.new_displacement[row];
                }
            }
            
            // Get or create target node
            DPNode* target_node = getOrCreateNode(target_state);
            
            // Update benefit
            updateBenefit(node, target_node, site, ext, prev_staples, row_start, cells_in_rows);
            valid_extensions++;
            
            // Update best node if this is the last site
            if (site + 1 == chip_info.total_sites) {
                for (int case_idx = 0; case_idx < 4; case_idx++) {
                    if (target_node->benefit[case_idx] > max_benefit) {
                        max_benefit = target_node->benefit[case_idx];
                        best_node = target_node;
                        Logger::log(Logger::INFO, "New best node found at final site: benefit = " + 
                                  std::to_string(max_benefit) + ", case = " + 
                                  std::to_string(case_idx));
                    }
                }
            }
        }
        
        // Check processing time again every 5000 extensions to catch extremely 
        // computation-intensive node processing
        if (extensions_count % 5000 == 0 && extensions_count > 0) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                current_time - process_start_time).count();
                
            if (elapsed_ms > max_processing_ms) {
                Logger::log(Logger::WARNING, "Extensions processing taking too long (" + 
                          std::to_string(elapsed_ms/1000) + " seconds), stopping early");
                break;
            }
        }
    }
    
    Logger::decreaseIndent();
    
    // Log completion summary
    if (extensions_count > 0) {
        Logger::log(Logger::INFO, "Site " + std::to_string(site) + " complete: " + 
                    std::to_string(extensions_count) + " total extensions, " + 
                    std::to_string(valid_extensions) + " valid extensions");
    }
    
    // Log memory usage only occasionally
    if (site % 500 == 0 || site == chip_info.total_sites) {
        Logger::log(Logger::INFO, "Current nodes in memory: " + std::to_string(all_nodes.size()));
    }
    
    // Update max nodes count stat
    max_nodes_count = std::max(max_nodes_count, all_nodes.size());
    
    // Report total processing time for performance monitoring
    auto current_time = std::chrono::high_resolution_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
        current_time - process_start_time).count();
        
    if (elapsed_ms > 5000) {  // Only log if took more than 5 seconds
        Logger::log(Logger::INFO, "Site " + std::to_string(site) + 
                  " processing took " + std::to_string(elapsed_ms/1000.0) + 
                  " seconds, queue size now: " + std::to_string(node_queue.size()));
    }
}

/**
 * @brief Generate all valid extensions for a node
 */
std::vector<Extension> DPSolver::generateExtensions(
    DPNode* node, 
    int site,
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start) {
    
    std::vector<Extension> extensions;
    Logger::increaseIndent();
    Logger::log(Logger::TRACE, "Generating extensions at site " + std::to_string(site));
    
    // Extract current state information
    int s1 = node->state.cell_offset[0];
    int l1 = node->state.displacement[0];
    bool f1 = node->state.is_flipped[0];
    
    int s2 = node->state.cell_offset[1];
    int l2 = node->state.displacement[1];
    bool f2 = node->state.is_flipped[1];
    
    int s3 = node->state.cell_offset[2];
    int l3 = node->state.displacement[2];
    bool f3 = node->state.is_flipped[2];
    
    Logger::log(Logger::TRACE, "Current state: s1=" + std::to_string(s1) + 
                ", l1=" + std::to_string(l1) + 
                ", f1=" + std::to_string(f1) + 
                ", s2=" + std::to_string(s2) + 
                ", l2=" + std::to_string(l2) + 
                ", f2=" + std::to_string(f2) + 
                ", s3=" + std::to_string(s3) + 
                ", l3=" + std::to_string(l3) + 
                ", f3=" + std::to_string(f3));
    
    // Helper function to create a base extension with current state
    auto createBaseExtension = [&]() {
        Extension ext(site + 1, -1);
        ext.new_cell_idx[0] = s1;
        ext.new_displacement[0] = l1 + 1;
        ext.new_is_flipped[0] = f1;
        
        ext.new_cell_idx[1] = s2;
        ext.new_displacement[1] = l2 + 1;
        ext.new_is_flipped[1] = f2;
        
        ext.new_cell_idx[2] = s3;
        ext.new_displacement[2] = l3 + 1;
        ext.new_is_flipped[2] = f3;
        
        return ext;
    };
    
    // Process row 1
    std::vector<Extension> row1_exts;
    
    // Option 1: No new cell
    if (cells_in_rows[0].empty() || 
        s1 >= static_cast<int>(cells_in_rows[0].size()) || 
        l1 < cell_types[cells_in_rows[0][s1]->type_index].cellSiteWidth) {
        Extension ext = createBaseExtension();
        row1_exts.push_back(ext);
        Logger::log("Row 1: Adding extension with no new cell");
    }
    
    // Option 2 & 3: Place next cell (with and without flipping)
    if (!cells_in_rows[0].empty() && s1 + 1 < static_cast<int>(cells_in_rows[0].size())) {
        Cell* next_cell = cells_in_rows[0][s1 + 1];
        int initial_site = next_cell->getInitialSite(chip_info.site_width);
        int current_site = site + 1;
        
        Logger::log("Row 1: Considering cell " + std::to_string(next_cell->cell_index) + 
                    " at site " + std::to_string(current_site) + 
                    " (initial site: " + std::to_string(initial_site) + 
                    ", max displacement: " + std::to_string(next_cell->max_displacement) + ")");
        
        // Check if placement is within displacement limit
        if (std::abs(current_site - initial_site) <= next_cell->max_displacement) {
            // Try without flipping
            if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                               next_cell->initial_y, false)) {
                Extension ext = createBaseExtension();
                ext.new_cell_idx[0] = s1 + 1;
                ext.new_displacement[0] = 0;
                ext.new_is_flipped[0] = false;
                row1_exts.push_back(ext);
                Logger::log("Row 1: Adding extension with new cell (not flipped)");
            } else {
                Logger::log("Row 1: Invalid placement for new cell (not flipped)");
            }
            
            // Try with flipping
            if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                               next_cell->initial_y, true)) {
                Extension ext = createBaseExtension();
                ext.new_cell_idx[0] = s1 + 1;
                ext.new_displacement[0] = 0;
                ext.new_is_flipped[0] = true;
                row1_exts.push_back(ext);
                Logger::log("Row 1: Adding extension with new cell (flipped)");
            } else {
                Logger::log("Row 1: Invalid placement for new cell (flipped)");
            }
        } else {
            Logger::log("Row 1: Cell exceeds displacement limit");
        }
    }
    
    Logger::log("Row 1 extensions: " + std::to_string(row1_exts.size()));
    
    // Process row 2
    std::vector<Extension> row2_exts;
    
    for (const Extension& row1_ext : row1_exts) {
        // Option 1: No new cell
        if (cells_in_rows[1].empty() || 
            s2 >= static_cast<int>(cells_in_rows[1].size()) || 
            l2 < cell_types[cells_in_rows[1][s2]->type_index].cellSiteWidth) {
            Extension ext = row1_ext;
            row2_exts.push_back(ext);
            // Don't log every extension to avoid excessive logs
        }
        
        // Option 2 & 3: Place next cell (with and without flipping)
        if (!cells_in_rows[1].empty() && s2 + 1 < static_cast<int>(cells_in_rows[1].size())) {
            Cell* next_cell = cells_in_rows[1][s2 + 1];
            int initial_site = next_cell->getInitialSite(chip_info.site_width);
            int current_site = site + 1;
            
            // Only log if this is the first extension from row 1
            if (row2_exts.empty()) {
                Logger::log("Row 2: Considering cell " + std::to_string(next_cell->cell_index) + 
                            " at site " + std::to_string(current_site));
            }
            
            // Check if placement is within displacement limit
            if (std::abs(current_site - initial_site) <= next_cell->max_displacement) {
                // Try without flipping
                if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                                   next_cell->initial_y, false)) {
                    Extension ext = row1_ext;
                    ext.new_cell_idx[1] = s2 + 1;
                    ext.new_displacement[1] = 0;
                    ext.new_is_flipped[1] = false;
                    row2_exts.push_back(ext);
                }
                
                // Try with flipping
                if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                                   next_cell->initial_y, true)) {
                    Extension ext = row1_ext;
                    ext.new_cell_idx[1] = s2 + 1;
                    ext.new_displacement[1] = 0;
                    ext.new_is_flipped[1] = true;
                    row2_exts.push_back(ext);
                }
            }
        }
    }
    
    Logger::log("Row 2 extensions: " + std::to_string(row2_exts.size()));
    
    // Process row 3
    for (const Extension& row2_ext : row2_exts) {
        // Option 1: No new cell
        if (cells_in_rows[2].empty() || 
            s3 >= static_cast<int>(cells_in_rows[2].size()) || 
            l3 < cell_types[cells_in_rows[2][s3]->type_index].cellSiteWidth) {
            Extension ext = row2_ext;
            extensions.push_back(ext);
        }
        
        // Option 2 & 3: Place next cell (with and without flipping)
        if (!cells_in_rows[2].empty() && s3 + 1 < static_cast<int>(cells_in_rows[2].size())) {
            Cell* next_cell = cells_in_rows[2][s3 + 1];
            int initial_site = next_cell->getInitialSite(chip_info.site_width);
            int current_site = site + 1;
            
            // Only log if this is the first extension
            if (extensions.empty()) {
                Logger::log("Row 3: Considering cell " + std::to_string(next_cell->cell_index) + 
                            " at site " + std::to_string(current_site));
            }
            
            // Check if placement is within displacement limit
            if (std::abs(current_site - initial_site) <= next_cell->max_displacement) {
                // Try without flipping
                if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                                   next_cell->initial_y, false)) {
                    Extension ext = row2_ext;
                    ext.new_cell_idx[2] = s3 + 1;
                    ext.new_displacement[2] = 0;
                    ext.new_is_flipped[2] = false;
                    extensions.push_back(ext);
                }
                
                // Try with flipping
                if (isValidPlacement(next_cell, current_site * chip_info.site_width, 
                                   next_cell->initial_y, true)) {
                    Extension ext = row2_ext;
                    ext.new_cell_idx[2] = s3 + 1;
                    ext.new_displacement[2] = 0;
                    ext.new_is_flipped[2] = true;
                    extensions.push_back(ext);
                }
            }
        }
    }
    
    Logger::log("Row 3 extensions: " + std::to_string(extensions.size()));
    
    // Calculate staple benefits for each extension
    int staple_possibilities = 0;
    for (Extension& ext : extensions) {
        // Calculate potential staple positions
        bool can_insert_r1_r2 = canInsertStaple(site + 1, 0, cells_in_rows, node);
        bool can_insert_r2_r3 = canInsertStaple(site + 1, 1, cells_in_rows, node);
        
        if (can_insert_r1_r2) staple_possibilities++;
        if (can_insert_r2_r3) staple_possibilities++;
        
        // Set staple benefits
        if (can_insert_r1_r2 && can_insert_r2_r3) {
            ext.staple_benefit = 2;
            ext.vdd_count = isVDDRow(row_start) ? 1 : 0;
            ext.vdd_count += isVDDRow(row_start + 1) ? 1 : 0;
            ext.vss_count = 2 - ext.vdd_count;
        } else if (can_insert_r1_r2) {
            ext.staple_benefit = 1;
            ext.vdd_count = isVDDRow(row_start) ? 1 : 0;
            ext.vss_count = 1 - ext.vdd_count;
        } else if (can_insert_r2_r3) {
            ext.staple_benefit = 1;
            ext.vdd_count = isVDDRow(row_start + 1) ? 1 : 0;
            ext.vss_count = 1 - ext.vdd_count;
        }
    }
    
    Logger::log("Staple possibilities: " + std::to_string(staple_possibilities) + 
                " (R1-R2: " + (canInsertStaple(site + 1, 0, cells_in_rows, node) ? "Yes" : "No") + 
                ", R2-R3: " + (canInsertStaple(site + 1, 1, cells_in_rows, node) ? "Yes" : "No") + ")");
    Logger::decreaseIndent();
    
    return extensions;
}

/**
 * @brief Check if a staple can be inserted at a given site
 */
bool DPSolver::canInsertStaple(int site, 
                             int row, 
                             const std::vector<std::vector<Cell*>>& cells_in_rows,
                             DPNode* node) {
    // Check if the two adjacent rows have empty space at the given site
    Logger::increaseIndent();
    
    // First row to check
    int row1 = row;
    if (row1 >= static_cast<int>(cells_in_rows.size())) {
        // Logger::log("Cannot insert staple: Invalid row1 index " + std::to_string(row1));
        // Logger::decreaseIndent();
        return false;
    }
    
    int s1 = node->state.cell_offset[row1];
    int l1 = node->state.displacement[row1];
    bool f1 = node->state.is_flipped[row1];
    
    // Check if there's a cell crossing the site in row1
    if (!cells_in_rows[row1].empty() && 
        s1 < static_cast<int>(cells_in_rows[row1].size())) {
        Cell* cell1 = cells_in_rows[row1][s1];
        const CellType& type1 = cell_types[cell1->type_index];
        int site_width = type1.cellSiteWidth;
        
        if (l1 < site_width) {
            // Cell crosses the site, check if there's a pin
            int relative_site = site - (cell1->initial_x / chip_info.site_width + l1);
            if (type1.hasPinAt(relative_site, f1)) {
                // Logger::log("Cannot insert staple at site " + std::to_string(site) + 
                //           ", row " + std::to_string(row) + ": pin collision in row1");
                // Logger::decreaseIndent();
                return false;
            }
        }
    }
    
    // Second row to check
    int row2 = row + 1;
    if (row2 >= static_cast<int>(cells_in_rows.size())) {
        // Logger::log("Cannot insert staple: Invalid row2 index " + std::to_string(row2));
        // Logger::decreaseIndent();
        return false;
    }
    
    int s2 = node->state.cell_offset[row2];
    int l2 = node->state.displacement[row2];
    bool f2 = node->state.is_flipped[row2];
    
    // Check if there's a cell crossing the site in row2
    if (!cells_in_rows[row2].empty() && 
        s2 < static_cast<int>(cells_in_rows[row2].size())) {
        Cell* cell2 = cells_in_rows[row2][s2];
        const CellType& type2 = cell_types[cell2->type_index];
        int site_width = type2.cellSiteWidth;
        
        if (l2 < site_width) {
            // Cell crosses the site, check if there's a pin
            int relative_site = site - (cell2->initial_x / chip_info.site_width + l2);
            if (type2.hasPinAt(relative_site, f2)) {
                // Logger::log("Cannot insert staple at site " + std::to_string(site) + 
                //           ", row " + std::to_string(row) + ": pin collision in row2");
                // Logger::decreaseIndent();
                return false;
            }
        }
    }
    
    // If we reach here, there's no obstruction
    // Logger::decreaseIndent();
    return true;
}

/**
 * @brief Check for staggering violations
 */
bool DPSolver::hasStaggeringViolation(int site,
                                    int staple_case,
                                    int prev_case,
                                    const std::vector<Staple>& prev_staples,
                                    int row_start) {
    // Check for staggering between current staples
    if ((staple_case == R1_R2_STAPLE && prev_case == R2_R3_STAPLE) ||
        (staple_case == R2_R3_STAPLE && prev_case == R1_R2_STAPLE)) {
        Logger::log("Staggering violation: Between R1-R2 and R2-R3 staples");
        return true;  // Staggering between R1-R2 and R2-R3
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
                Logger::log("Staggering violation: With previous staple at site " + 
                          std::to_string(site) + ", row " + std::to_string(prev_row));
                return true;
            }
        }
        else if (prev_staple.x == next_site_x) {
            // Previous staple at next site
            int prev_row = prev_staple.y / chip_info.row_height;
            
            if ((staple_case == R1_R2_STAPLE && prev_row == row_start + 3) ||
                (staple_case == R2_R3_STAPLE && prev_row == row_start - 1)) {
                Logger::log("Staggering violation: With previous staple at site " + 
                          std::to_string(site + 1) + ", row " + std::to_string(prev_row));
                return true;
            }
        }
    }
    
    return false;
}

/**
 * @brief Get or create a node with given state
 */
DPNode* DPSolver::getOrCreateNode(const CompactState& state) {
    auto it = node_lookup_table.find(state);
    if (it != node_lookup_table.end()) {
        // Log only occasionally to avoid excessive output
        if (nodes_processed % 1000 == 0) {
            Logger::log("Node found in lookup table: site " + std::to_string(state.site));
        }
        return it->second;
    }
    
    // Create new node
    DPNode* new_node = new DPNode(state);
    
    // Log only occasionally to avoid excessive output
    if (nodes_created % 1000 == 0 || nodes_created < 10) {
        Logger::log("Creating new node #" + std::to_string(nodes_created) + 
                  " at site " + std::to_string(state.site) + 
                  " with offsets (" + 
                  std::to_string(state.cell_offset[0]) + "," + 
                  std::to_string(state.cell_offset[1]) + "," + 
                  std::to_string(state.cell_offset[2]) + ")");
    }
    
    // Add to lookup table and queue
    node_lookup_table[state] = new_node;
    node_queue.push(new_node);
    all_nodes.push_back(new_node);
    nodes_created++;
    
    return new_node;
}

/**
 * @brief Update node benefit based on extension
 */
void DPSolver::updateBenefit(DPNode* from_node,
                           DPNode* to_node,
                           int site,
                           const Extension& extension,
                           const std::vector<Staple>& prev_staples,
                           int row_start,
                           const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    // Try all combinations of staple insertion cases
    for (int from_case = 0; from_case < 4; from_case++) {
        // Skip invalid source cases
        if (from_node->benefit[from_case] < -900000) {
            continue;
        }
        
        // Determine valid target cases based on potential staple insertions
        std::vector<int> valid_target_cases;
        valid_target_cases.push_back(NO_STAPLE); // Always consider no staple case
        
        bool can_insert_r1_r2 = canInsertStaple(site + 1, 0, cells_in_rows, to_node);
        bool can_insert_r2_r3 = canInsertStaple(site + 1, 1, cells_in_rows, to_node);
        
        if (can_insert_r1_r2) {
            valid_target_cases.push_back(R1_R2_STAPLE);
        }
        
        if (can_insert_r2_r3) {
            valid_target_cases.push_back(R2_R3_STAPLE);
        }
        
        if (can_insert_r1_r2 && can_insert_r2_r3) {
            valid_target_cases.push_back(BOTH_STAPLES);
        }
        
        // Log the valid cases for the first few iterations
        if (nodes_processed < 5 || nodes_processed % 1000 == 0) {
            std::string valid_cases_str = "Valid target cases:";
            for (int case_idx : valid_target_cases) {
                valid_cases_str += " " + std::to_string(case_idx);
            }
            Logger::log(Logger::TRACE, valid_cases_str);
        }
        
        // Process each valid target case
        for (int to_case : valid_target_cases) {
            // Check for staggering violations
            if (to_case != NO_STAPLE && 
                hasStaggeringViolation(site + 1, to_case, from_case, prev_staples, row_start)) {
                continue;  // Staggering violation, skip this case
            }
            
            // Calculate staple benefits based on the case
            int staple_benefit = 0;
            int vdd_staples = 0;
            int vss_staples = 0;
            
            switch (to_case) {
                case NO_STAPLE:
                    // No staples inserted
                    // Check if we should apply balance factor for reserving space
                    if (site + 1 < chip_info.total_sites) {
                        // Calculate a heuristic benefit for leaving space
                        // This would encourage balanced staple insertion
                        double current_vdd = from_node->vdd_staples[from_case];
                        double current_vss = from_node->vss_staples[from_case];
                        
                        // Only apply if there's significant imbalance
                        if (current_vdd > 0 && current_vss > 0) {
                            double ratio = std::max(current_vdd, current_vss) / 
                                         std::min(current_vdd, current_vss);
                            
                            if (ratio > 1.05) {
                                // Apply balance factor
                                staple_benefit = 0;
                                
                                // Add a tiny bonus to encourage reservation for the less frequent staple type
                                if ((current_vdd > current_vss && can_insert_r1_r2 && !isVDDRow(row_start)) ||
                                    (current_vdd > current_vss && can_insert_r2_r3 && !isVDDRow(row_start + 1)) ||
                                    (current_vss > current_vdd && can_insert_r1_r2 && isVDDRow(row_start)) ||
                                    (current_vss > current_vdd && can_insert_r2_r3 && isVDDRow(row_start + 1))) {
                                    staple_benefit = params.balance_factor;
                                    Logger::log("Applied balance factor for NO_STAPLE case");
                                }
                            }
                        }
                    }
                    break;
                
                case R1_R2_STAPLE:
                    // Insert staple between R1 and R2
                    staple_benefit = 1;
                    if (isVDDRow(row_start)) {
                        vdd_staples = 1;
                    } else {
                        vss_staples = 1;
                    }
                    break;
                
                case R2_R3_STAPLE:
                    // Insert staple between R2 and R3
                    staple_benefit = 1;
                    if (isVDDRow(row_start + 1)) {
                        vdd_staples = 1;
                    } else {
                        vss_staples = 1;
                    }
                    break;
                
                case BOTH_STAPLES:
                    // Insert staples between R1-R2 and R2-R3
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
                    break;
            }
            
            // Calculate new benefit
            int new_benefit = from_node->benefit[from_case] + staple_benefit;
            int new_vdd = from_node->vdd_staples[from_case] + vdd_staples;
            int new_vss = from_node->vss_staples[from_case] + vss_staples;
            
            // Apply dynamic staple balance adjustment
            if (new_vdd > 0 && new_vss > 0) {
                double ratio = static_cast<double>(std::max(new_vdd, new_vss)) / 
                              std::min(new_vdd, new_vss);
                
                // Apply stronger balance encouragement as ratio approaches the limit
                if (ratio > 1.05) {
                    // Apply balance encouragement with dynamic weight
                    double balance_weight = params.balance_factor * (1.0 + (ratio - 1.05) * 2.0);
                    
                    if (new_vdd > new_vss) {
                        // Encourage VSS staples
                        if (to_case == R1_R2_STAPLE && !isVDDRow(row_start)) {
                            new_benefit += balance_weight;
                            Logger::log("Encouraged VSS staple: +" + std::to_string(balance_weight));
                        }
                        else if (to_case == R2_R3_STAPLE && !isVDDRow(row_start + 1)) {
                            new_benefit += balance_weight;
                            Logger::log("Encouraged VSS staple: +" + std::to_string(balance_weight));
                        }
                    }
                    else if (new_vss > new_vdd) {
                        // Encourage VDD staples
                        if (to_case == R1_R2_STAPLE && isVDDRow(row_start)) {
                            new_benefit += balance_weight;
                            Logger::log("Encouraged VDD staple: +" + std::to_string(balance_weight));
                        }
                        else if (to_case == R2_R3_STAPLE && isVDDRow(row_start + 1)) {
                            new_benefit += balance_weight;
                            Logger::log("Encouraged VDD staple: +" + std::to_string(balance_weight));
                        }
                    }
                }
            }
            
            // Update benefit if better
            if (new_benefit > to_node->benefit[to_case]) {
                // Log significant benefit improvements
                if ((nodes_processed < 10 && staple_benefit > 0) || 
                    (new_benefit > to_node->benefit[to_case] + 5)) {
                    Logger::log(Logger::DEBUG, "Significant benefit improvement: case " + std::to_string(to_case) + 
                                " from " + std::to_string(to_node->benefit[to_case]) + 
                                " to " + std::to_string(new_benefit) + 
                                " (VDD: " + std::to_string(new_vdd) + 
                                ", VSS: " + std::to_string(new_vss) + ")");
                }
                
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
 * @brief Backtrack to get the optimal solution
 */
std::vector<Staple> DPSolver::backtrack(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start) {
    
    std::vector<Staple> staples;
    if (!best_node) {
        Logger::log(Logger::WARNING, "No solution found, nothing to backtrack");
        return staples;
    }
    
    // Find the case with maximum benefit
    int best_case = 0;
    int max_benefit = best_node->benefit[0];
    
    for (int i = 1; i < 4; i++) {
        if (best_node->benefit[i] > max_benefit) {
            max_benefit = best_node->benefit[i];
            best_case = i;
        }
    }
    
    Logger::log(Logger::INFO, "Backtracking from site " + std::to_string(best_node->state.site) + 
                " with best case " + std::to_string(best_case) + 
                " (benefit: " + std::to_string(max_benefit) + 
                ", VDD: " + std::to_string(best_node->vdd_staples[best_case]) + 
                ", VSS: " + std::to_string(best_node->vss_staples[best_case]) + ")");
    
    // Prepare to backtrack to update cell positions
    std::vector<std::pair<Cell*, std::pair<int, bool>>> cell_updates; // <cell, <new_x, flipped>>
    
    // Backtrack from best node
    DPNode* current = best_node;
    int current_case = best_case;
    int site = chip_info.total_sites;
    
    Logger::increaseIndent();
    Logger::log("Starting backtracking from site " + std::to_string(site) + 
                " with best case " + std::to_string(best_case) + 
                " (benefit: " + std::to_string(max_benefit) + ")");
    
    int step = 0;
    while (current != nullptr && current->prev_node[current_case] != nullptr) {
        step++;
        Logger::log("Backtrack step " + std::to_string(step) + ", Site: " + std::to_string(site) + 
                    ", Case: " + std::to_string(current_case));
        
        // Get previous node and case
        DPNode* prev_node = current->prev_node[current_case];
        int prev_case = current->case_from_prev[current_case];
        
        // Insert staples based on the case
        switch (current_case) {
            case R1_R2_STAPLE: {
                // Insert staple between R1 and R2
                int x = site * chip_info.site_width;
                int y = (row_start + 1) * chip_info.row_height;
                bool is_vdd = isVDDRow(row_start);
                staples.push_back(Staple(x, y, is_vdd));
                Logger::log("Inserting R1-R2 staple at site " + std::to_string(site) + 
                            ", position (" + std::to_string(x) + "," + std::to_string(y) + "), " +
                            (is_vdd ? "VDD" : "VSS") + " type");
                break;
            }
            
            case R2_R3_STAPLE: {
                // Insert staple between R2 and R3
                int x = site * chip_info.site_width;
                int y = (row_start + 2) * chip_info.row_height;
                bool is_vdd = isVDDRow(row_start + 1);
                staples.push_back(Staple(x, y, is_vdd));
                Logger::log("Inserting R2-R3 staple at site " + std::to_string(site) + 
                            ", position (" + std::to_string(x) + "," + std::to_string(y) + "), " +
                            (is_vdd ? "VDD" : "VSS") + " type");
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
                
                Logger::log("Inserting both staples at site " + std::to_string(site) + ": " +
                            "R1-R2 (" + (is_vdd1 ? "VDD" : "VSS") + "), " +
                            "R2-R3 (" + (is_vdd2 ? "VDD" : "VSS") + ")");
                break;
            }
            
            case NO_STAPLE:
                // No staple inserted
                Logger::log("No staple at site " + std::to_string(site));
                break;
        }
        
        // Update cell positions based on state change
        for (int row = 0; row < 3; row++) {
            // Check if cell assignment changed
            int prev_offset = prev_node->state.cell_offset[row];
            int curr_offset = current->state.cell_offset[row];
            bool prev_flipped = prev_node->state.is_flipped[row];
            bool curr_flipped = current->state.is_flipped[row];
            
            if (prev_offset != curr_offset || prev_flipped != curr_flipped) {
                // Cell assignment changed, need to update
                if (!cells_in_rows[row].empty() && 
                    curr_offset >= 0 && curr_offset < static_cast<int>(cells_in_rows[row].size())) {
                    Cell* cell = cells_in_rows[row][curr_offset];
                    int new_x = site * chip_info.site_width;
                    
                    // Record cell update
                    cell_updates.push_back(std::make_pair(cell, std::make_pair(new_x, curr_flipped)));
                    
                    Logger::log("Updating cell " + std::to_string(cell->cell_index) + 
                                " in row " + std::to_string(row_start + row) + 
                                " to position " + std::to_string(new_x) + 
                                (curr_flipped ? " (flipped)" : ""));
                }
            }
        }
        
        // Move to previous node
        site--;
        current = prev_node;
        current_case = prev_case;
    }
    
    // Apply cell position updates
    int total_updates = 0;
    for (const auto& update : cell_updates) {
        Cell* cell = update.first;
        int new_x = update.second.first;
        bool flipped = update.second.second;
        
        // Update cell position only if it's different
        if (cell->current_x != new_x || cell->is_flipped != flipped) {
            Logger::log("Cell " + std::to_string(cell->cell_index) + ": " +
                        "Position " + std::to_string(cell->current_x) + "->" + std::to_string(new_x) + 
                        ", Flipped " + std::to_string(cell->is_flipped) + "->" + std::to_string(flipped));
            cell->current_x = new_x;
            cell->is_flipped = flipped;
            total_updates++;
        }
    }
    
    Logger::log("Cell updates complete: " + std::to_string(total_updates) + " cells modified");
    
    // Reverse staples vector since we added them in reverse order
    std::reverse(staples.begin(), staples.end());
    
    // Log staple balance
    int vdd_count = 0, vss_count = 0;
    for (const Staple& staple : staples) {
        if (staple.is_vdd) vdd_count++;
        else vss_count++;
    }
    
    double ratio = 0.0;
    if (vdd_count > 0 && vss_count > 0) {
        ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                std::min(vdd_count, vss_count);
    }
    
    Logger::log("Backtracking complete: " + std::to_string(staples.size()) + " staples inserted");
    Logger::log("Final staple balance: VDD=" + std::to_string(vdd_count) + 
                ", VSS=" + std::to_string(vss_count) + 
                ", Ratio=" + std::to_string(ratio));
    Logger::decreaseIndent();
    
    return staples;
}

/**
 * @brief Apply final cell placement updates
 */
void DPSolver::applyPlacementUpdates(
    const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    Logger::log("Applying final cell placement updates");
    Logger::increaseIndent();
    
    int total_cells = 0;
    int cells_moved = 0;
    int cells_flipped = 0;
    
    // Count total cells and check for any final updates needed
    for (size_t row = 0; row < cells_in_rows.size(); row++) {
        for (size_t i = 0; i < cells_in_rows[row].size(); i++) {
            Cell* cell = cells_in_rows[row][i];
            total_cells++;
            
            // Check if cell has been moved
            if (cell->current_x != cell->initial_x) {
                cells_moved++;
                
                // Calculate displacement in sites
                int displacement = std::abs(cell->current_x - cell->initial_x) / chip_info.site_width;
                if (displacement > 0) {
                    Logger::log("Cell " + std::to_string(cell->cell_index) + 
                              " moved by " + std::to_string(displacement) + 
                              " sites (max allowed: " + std::to_string(cell->max_displacement) + ")");
                }
            }
            
            // Check if cell has been flipped
            if (cell->is_flipped) {
                cells_flipped++;
            }
        }
    }
    
    Logger::log("Final placement stats:");
    Logger::log("  Total cells: " + std::to_string(total_cells));
    Logger::log("  Cells moved: " + std::to_string(cells_moved) + 
                " (" + std::to_string(100.0 * cells_moved / total_cells) + "%)");
    Logger::log("  Cells flipped: " + std::to_string(cells_flipped) + 
                " (" + std::to_string(100.0 * cells_flipped / total_cells) + "%)");
    
    Logger::decreaseIndent();
}

/**
 * @brief Check if cell placement is valid
 */
bool DPSolver::isValidPlacement(const Cell* cell, int x, int y, bool is_flipped) {
    // Check if the placement is within displacement limit
    int displacement = std::abs(x - cell->initial_x);
    int displacement_in_sites = displacement / chip_info.site_width;
    
    // Log only occasionally to avoid excessive output
    bool should_log = (nodes_processed < 10) || (nodes_processed % 1000 == 0 && nodes_processed < 10000);
    
    if (displacement_in_sites > cell->max_displacement) {
        if (should_log) {
            Logger::log("Invalid placement: Cell " + std::to_string(cell->cell_index) + 
                      " displacement " + std::to_string(displacement_in_sites) + 
                      " exceeds max " + std::to_string(cell->max_displacement));
        }
        return false;
    }
    
    // Check if the placement is in the same row
    if (y != cell->initial_y) {
        if (should_log) {
            Logger::log("Invalid placement: Cell " + std::to_string(cell->cell_index) + 
                      " row changed from " + std::to_string(cell->initial_y) + 
                      " to " + std::to_string(y));
        }
        return false;
    }
    
    // Additional checks can be added here if needed
    
    return true;
}

/**
 * @brief Get compact state encoding
 */
CompactState DPSolver::getCompactState(
    int site,
    int s1, int l1, bool f1,
    int s2, int l2, bool f2,
    int s3, int l3, bool f3,
    const std::vector<std::vector<Cell*>>& cells_in_rows) {
    
    CompactState state;
    state.site = site;
    
    // Calculate compact encoding for each row
    for (int row = 0; row < 3; row++) {
        int s = (row == 0) ? s1 : (row == 1) ? s2 : s3;
        int l = (row == 0) ? l1 : (row == 1) ? l2 : l3;
        bool f = (row == 0) ? f1 : (row == 1) ? f2 : f3;
        
        // Calculate cell offset relative to initial placement
        if (cells_in_rows[row].empty()) {
            state.cell_offset[row] = 0;
            state.displacement[row] = 0;
            state.is_flipped[row] = false;
        } else {
            // Find reference cell in initial placement
            int initial_s = 0;
            for (size_t i = 0; i < cells_in_rows[row].size(); i++) {
                Cell* cell = cells_in_rows[row][i];
                int cell_start_site = cell->initial_x / chip_info.site_width;
                int cell_end_site = cell_start_site + cell_types[cell->type_index].cellSiteWidth;
                
                if (cell_start_site <= site && site < cell_end_site) {
                    initial_s = i;
                    break;
                }
                
                if (cell_start_site > site) {
                    if (i > 0) initial_s = i - 1;
                    break;
                }
            }
            
            // Calculate relative offset
            state.cell_offset[row] = s - initial_s;
            state.displacement[row] = l;
            state.is_flipped[row] = f;
        }
    }
    
    return state;
}

/**
 * @brief Check if row is VDD or VSS
 */
bool DPSolver::isVDDRow(int row_idx) const {
    // Even rows (0, 2, 4, ...) are VDD, odd rows are VSS
    bool is_vdd = (row_idx % 2 == 0);
    
    // Only log during the first few iterations to avoid excessive output
    if (nodes_processed < 5) {
        Logger::log("Row " + std::to_string(row_idx) + " is " + 
                  (is_vdd ? "VDD" : "VSS") + " row");
    }
    
    return is_vdd;
}

/**
 * @brief Calculate staple benefit with balance constraints
 */
double DPSolver::calculateStapleBenefit(
    int staple_case,
    int current_vdd,
    int current_vss,
    int row_start) {
    
    double base_benefit = 1.0;
    bool should_log = (nodes_processed < 20) || (nodes_processed % 10000 == 0);
    
    // If no staples yet, use base benefit
    if (current_vdd == 0 || current_vss == 0) {
        if (should_log) {
            Logger::log("Base benefit " + std::to_string(base_benefit) + 
                        " (no balance adjustment: " + 
                        std::to_string(current_vdd) + " VDD, " + 
                        std::to_string(current_vss) + " VSS)");
        }
        return base_benefit;
    }
    
    // Calculate current ratio
    double ratio = static_cast<double>(std::max(current_vdd, current_vss)) / 
                  std::min(current_vdd, current_vss);
    
    // If ratio is approaching the limit, adjust benefit
    if (ratio > 1.05) {
        double adjustment = (ratio - 1.05) * 2.0;
        
        // Check if this staple helps balance
        bool helps_balance = false;
        
        if (current_vdd > current_vss) {
            // Need more VSS staples
            switch (staple_case) {
                case R1_R2_STAPLE:
                    helps_balance = !isVDDRow(row_start);
                    break;
                case R2_R3_STAPLE:
                    helps_balance = !isVDDRow(row_start + 1);
                    break;
                case BOTH_STAPLES:
                    // Check how many VSS staples would be added
                    int new_vss = 0;
                    if (!isVDDRow(row_start)) new_vss++;
                    if (!isVDDRow(row_start + 1)) new_vss++;
                    helps_balance = (new_vss > 0);
                    break;
            }
        } else {
            // Need more VDD staples
            switch (staple_case) {
                case R1_R2_STAPLE:
                    helps_balance = isVDDRow(row_start);
                    break;
                case R2_R3_STAPLE:
                    helps_balance = isVDDRow(row_start + 1);
                    break;
                case BOTH_STAPLES:
                    // Check how many VDD staples would be added
                    int new_vdd = 0;
                    if (isVDDRow(row_start)) new_vdd++;
                    if (isVDDRow(row_start + 1)) new_vdd++;
                    helps_balance = (new_vdd > 0);
                    break;
            }
        }
        
        // Adjust benefit
        double adjusted_benefit;
        if (helps_balance) {
            adjusted_benefit = base_benefit * (1.0 + adjustment);
            if (should_log) {
                Logger::log("Increased benefit: " + std::to_string(base_benefit) + " -> " + 
                           std::to_string(adjusted_benefit) + " (helps balance, ratio: " + 
                           std::to_string(ratio) + ")");
            }
        } else {
            adjusted_benefit = base_benefit * (1.0 - adjustment * 0.5);
            if (should_log) {
                Logger::log("Decreased benefit: " + std::to_string(base_benefit) + " -> " + 
                           std::to_string(adjusted_benefit) + " (hurts balance, ratio: " + 
                           std::to_string(ratio) + ")");
            }
        }
        return adjusted_benefit;
    }
    
    if (should_log) {
        Logger::log("Using base benefit " + std::to_string(base_benefit) + 
                   " (balanced: " + std::to_string(current_vdd) + " VDD, " + 
                   std::to_string(current_vss) + " VSS, ratio: " + 
                   std::to_string(ratio) + ")");
    }
    return base_benefit;
}

/**
 * @brief Clean up allocated memory
 */
void DPSolver::cleanup() {
    size_t nodes_to_delete = all_nodes.size();
    Logger::log("Cleaning up " + std::to_string(nodes_to_delete) + " nodes");
    
    // Delete all allocated nodes
    for (DPNode* node : all_nodes) {
        delete node;
    }
    
    // Clear containers
    all_nodes.clear();
    node_lookup_table.clear();
    while (!node_queue.empty()) {
        node_queue.pop();
    }
    
    // Reset node pointers
    best_node = nullptr;
    
    Logger::log("Memory cleanup complete");
}

/**
 * @brief Monitor memory usage and report statistics
 */
void DPSolver::monitorMemoryUsage() {
    size_t node_count = all_nodes.size();
    size_t lookup_table_size = node_lookup_table.size();
    size_t queue_size = node_queue.size();
    size_t estimated_memory_bytes = node_count * sizeof(DPNode);
    double estimated_memory_mb = static_cast<double>(estimated_memory_bytes) / (1024 * 1024);
    
    Logger::log("Memory stats:");
    Logger::log("  Nodes in memory: " + std::to_string(node_count));
    Logger::log("  Lookup table size: " + std::to_string(lookup_table_size));
    Logger::log("  Queue size: " + std::to_string(queue_size));
    Logger::log("  Estimated memory usage: " + std::to_string(estimated_memory_mb) + " MB");
    
    // Warning if memory usage is high
    if (estimated_memory_bytes > 12ULL * 1024 * 1024 * 1024) {
        Logger::log("WARNING: Memory usage approaching limit of 16GB!");
        std::cout << "WARNING: Memory usage approaching limit of 16GB!" << std::endl;
    }
}