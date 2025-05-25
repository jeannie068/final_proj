/**
 * @file dp_solver.cpp
 * @brief Implementation of DAG-based dynamic programming solver for MATRO
 */

#include "dp_solver.hpp"
#include "../Logger.hpp"
#include "../memory_usage.hpp"
#include <algorithm>
#include <limits>
#include <thread>

/**
 * @brief Constructor
 */
DPSolver::DPSolver(const ChipInfo& chip_info, 
                   const std::vector<CellType>& cell_types,
                   const AlgorithmParams& params)
    : chip_info(chip_info), cell_types(cell_types), params(params),
      best_final_node(nullptr), best_final_case(-1), 
      nodes_created(0), max_nodes_in_memory(0) {
    
    Logger::log(Logger::INFO, "DPSolver initialized");
    Logger::log(Logger::INFO, "  Total sites: " + std::to_string(chip_info.total_sites));
    Logger::log(Logger::INFO, "  Balance factor: " + std::to_string(params.balance_factor));
}

/**
 * @brief Destructor
 */
DPSolver::~DPSolver() {
    cleanup();
}

// /**
//  * @brief Robust solver for triple-row optimization with timeout protection and fallback
//  */
// std::vector<Staple> DPSolver::solveTripleRow(
//     const std::vector<std::vector<Cell*>>& cells_in_rows,
//     int row_start,
//     const std::vector<Staple>& prev_staples) {
    
//     auto start_time = std::chrono::high_resolution_clock::now();
//     size_t start_memory = getCurrentMemoryUsage();
    
//     // === INITIALIZATION PHASE ===
//     std::cout << "=== STARTING TRIPLE-ROW " << row_start << "-" << (row_start+2) 
//               << " (Memory: " << start_memory << "MB) ===" << std::endl;
    
//     Logger::log(Logger::INFO, "Starting triple-row optimization for rows " + 
//                 std::to_string(row_start) + "-" + std::to_string(row_start + 2));
//     Logger::log(Logger::INFO, "Memory at start: " + std::to_string(start_memory) + " MB");
//     Logger::log(Logger::INFO, "Previous staples: " + std::to_string(prev_staples.size()));
    
//     // 強制完全清理 - 確保clean state
//     forceCompleteCleanup();
    
//     // 設定timeout - 每個subproblem最多2分鐘
//     const int MAX_SUBPROBLEM_TIME_SECONDS = 120;
//     const int PROGRESS_REPORT_INTERVAL = 20; // 每20個sites報告一次
    
//     // Initialize tracking variables
//     bool timeout_occurred = false;
//     bool algorithm_stuck = false;
//     int last_reported_site = -1;
//     size_t max_memory_this_run = start_memory;
    
//     try {
//         // === SETUP PHASE ===
//         computeInitialS(cells_in_rows);
        
//         // Log cell information
//         for (int row = 0; row < 3; row++) {
//             std::cout << "Row " << (row_start + row) << ": " << cells_in_rows[row].size() 
//                       << " cells" << std::endl;
//             Logger::log(Logger::INFO, "Row " + std::to_string(row_start + row) + 
//                         ": " + std::to_string(cells_in_rows[row].size()) + " cells");
//         }
        
//         // Create source node
//         DPNode* source = createSourceNode();
//         node_queue.push(source);
//         all_nodes.push_back(source);
//         nodes_created = 1;
        
//         best_final_node = nullptr;
//         best_final_case = -1;
        
//         std::cout << "Processing sites 0-" << chip_info.total_sites << "..." << std::endl;
        
//         // === MAIN PROCESSING LOOP ===
//         for (int site = 0; site <= chip_info.total_sites; site++) {
//             // === TIMEOUT CHECK ===
//             auto current_time = std::chrono::high_resolution_clock::now();
//             auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
//                 current_time - start_time).count();
            
//             if (elapsed_seconds > MAX_SUBPROBLEM_TIME_SECONDS) {
//                 timeout_occurred = true;
//                 std::cout << "TIMEOUT: Subproblem exceeded " << MAX_SUBPROBLEM_TIME_SECONDS 
//                           << "s at site " << site << std::endl;
//                 Logger::log(Logger::ERROR, "TIMEOUT at site " + std::to_string(site) + 
//                            " after " + std::to_string(elapsed_seconds) + " seconds");
//                 break;
//             }
            
//             // === PROGRESS REPORTING ===
//             if (site % PROGRESS_REPORT_INTERVAL == 0 || site == chip_info.total_sites) {
//                 size_t current_memory = getCurrentMemoryUsage();
//                 max_memory_this_run = std::max(max_memory_this_run, current_memory);
                
//                 size_t queue_size = node_queue.size();
//                 std::cout << "Site " << site << "/" << chip_info.total_sites 
//                           << " (" << queue_size << " nodes, " << current_memory 
//                           << "MB, " << elapsed_seconds << "s)" << std::endl;
                
//                 Logger::log(Logger::INFO, "Progress: site " + std::to_string(site) + 
//                            ", nodes: " + std::to_string(queue_size) + 
//                            ", memory: " + std::to_string(current_memory) + "MB");
                
//                 last_reported_site = site;
//             }
            
//             // === PROCESS CURRENT SITE ===
//             size_t nodes_at_site = node_queue.size();
            
//             if (nodes_at_site == 0) {
//                 // No nodes to process - check if we're stuck
//                 if (site > 50) { // Only worry if we're well into the process
//                     algorithm_stuck = true;
//                     std::cout << "STUCK: No nodes at site " << site << std::endl;
//                     Logger::log(Logger::ERROR, "Algorithm stuck at site " + std::to_string(site));
//                     break;
//                 }
//                 continue; // Early sites might naturally have no nodes
//             }
            
//             // === MEMORY PROTECTION ===
//             size_t current_memory = getCurrentMemoryUsage();
//             if (current_memory > start_memory + 1000) { // More than 1GB increase
//                 std::cout << "MEMORY WARNING: Usage increased by " 
//                           << (current_memory - start_memory) << "MB" << std::endl;
//                 Logger::log(Logger::WARNING, "High memory usage: " + std::to_string(current_memory) + "MB");
                
//                 // Force aggressive pruning
//                 if (all_nodes.size() > 1000) {
//                     aggressivePruning(site);
//                 }
//             }
            
//             // Process nodes at this site with timeout per site
//             auto site_start = std::chrono::high_resolution_clock::now();
//             processNodesAtSiteWithTimeout(site, cells_in_rows, prev_staples, row_start, 10); // 10s per site max
//             auto site_end = std::chrono::high_resolution_clock::now();
            
//             auto site_duration = std::chrono::duration_cast<std::chrono::seconds>(site_end - site_start);
//             if (site_duration.count() > 5) {
//                 Logger::log(Logger::WARNING, "Site " + std::to_string(site) + 
//                            " took " + std::to_string(site_duration.count()) + " seconds");
//             }
            
//             // Update memory tracking
//             max_nodes_in_memory = std::max(max_nodes_in_memory, all_nodes.size());
//         }
        
//         // === SOLUTION EXTRACTION ===
//         std::vector<Staple> solution;
        
//         if (timeout_occurred || algorithm_stuck) {
//             std::cout << "Using fallback solution generation..." << std::endl;
//             solution = generateFallbackSolution(cells_in_rows, row_start, prev_staples);
//         } else {
//             // Normal solution extraction
//             solution = extractOptimalSolution(cells_in_rows, row_start);
//         }
        
//         // === FINAL REPORTING ===
//         auto end_time = std::chrono::high_resolution_clock::now();
//         auto total_duration = std::chrono::duration_cast<std::chrono::seconds>(
//             end_time - start_time);
        
//         size_t end_memory = getCurrentMemoryUsage();
        
//         std::cout << "=== TRIPLE-ROW " << row_start << "-" << (row_start+2) 
//                   << " COMPLETED ===" << std::endl;
//         std::cout << "Time: " << total_duration.count() << "s, "
//                   << "Memory: " << start_memory << "->" << end_memory << "MB, "
//                   << "Staples: " << solution.size() << std::endl;
        
//         Logger::log(Logger::INFO, "Triple-row " + std::to_string(row_start) + "-" + 
//                    std::to_string(row_start+2) + " completed:");
//         Logger::log(Logger::INFO, "  Duration: " + std::to_string(total_duration.count()) + "s");
//         Logger::log(Logger::INFO, "  Memory: " + std::to_string(start_memory) + " -> " + 
//                    std::to_string(end_memory) + " MB");
//         Logger::log(Logger::INFO, "  Peak memory: " + std::to_string(max_memory_this_run) + "MB");
//         Logger::log(Logger::INFO, "  Staples: " + std::to_string(solution.size()));
//         Logger::log(Logger::INFO, "  Nodes created: " + std::to_string(nodes_created));
//         Logger::log(Logger::INFO, "  Max nodes in memory: " + std::to_string(max_nodes_in_memory));
        
//         if (timeout_occurred) {
//             Logger::log(Logger::WARNING, "  TIMEOUT occurred");
//         }
//         if (algorithm_stuck) {
//             Logger::log(Logger::WARNING, "  ALGORITHM STUCK");
//         }
        
//         // Apply cell placements (if any)
//         applyCellPlacements(cells_in_rows);
        
//         return solution;
        
//     } catch (const std::exception& e) {
//         std::cout << "EXCEPTION in triple-row solver: " << e.what() << std::endl;
//         Logger::log(Logger::ERROR, "Exception in solveTripleRow: " + std::string(e.what()));
        
//         // Force cleanup and return fallback solution
//         forceCompleteCleanup();
//         return generateFallbackSolution(cells_in_rows, row_start, prev_staples);
        
//     } catch (...) {
//         std::cout << "UNKNOWN EXCEPTION in triple-row solver" << std::endl;
//         Logger::log(Logger::ERROR, "Unknown exception in solveTripleRow");
        
//         // Force cleanup and return fallback solution
//         forceCompleteCleanup();
//         return generateFallbackSolution(cells_in_rows, row_start, prev_staples);
//     }
// }


/**
 * @brief Create initial source node
 */
DPNode* DPSolver::createSourceNode() {
    Logger::log(Logger::INFO, "Creating source node at site 0");
    
    // Paper's source node: (0, 0, 0, 0, 0) for double-row
    // For triple-row: (0, 0, 0, 0, 0, 0, 0)
    DPNode* source = new DPNode(0, 0, 0, 0, 0, 0, 0);
    source->benefit[CASE_1_NO_STAPLE] = 0;
    
    return source;
}

Cell* DPSolver::getNextCellToPlace(int row, int s_j, const std::vector<Cell*>& row_cells) {
    // s_j is the number of cells already placed in this row
    // So the next cell to consider is at index s_j
    if (s_j < static_cast<int>(row_cells.size())) {
        return row_cells[s_j];
    }
    return nullptr;
}

Cell* DPSolver::getLastPlacedCell(int row, int s_j, const std::vector<Cell*>& row_cells) {
    // The last placed cell is at index s_j - 1
    if (s_j > 0 && s_j <= static_cast<int>(row_cells.size())) {
        return row_cells[s_j - 1];
    }
    return nullptr;
}

/**
 * @brief Process all nodes at current site
 */
void DPSolver::processNodesAtSite(int site,
                                 const std::vector<std::vector<Cell*>>& cells_in_rows,
                                 const std::vector<Staple>& prev_staples,
                                 int row_start) {
    auto site_start_time = std::chrono::high_resolution_clock::now();
    
    // 原有的node processing邏輯...
    std::vector<DPNode*> nodes_at_site;
    size_t queue_size = node_queue.size();
    
    for (size_t i = 0; i < queue_size; i++) {
        DPNode* node = node_queue.front();
        node_queue.pop();
        
        if (node->site == site) {
            nodes_at_site.push_back(node);
        } else if (node->site > site) {
            node_queue.push(node);
        }
    }
    
    // Generate extensions for each node
    for (DPNode* node : nodes_at_site) {
        generateExtensions(node, cells_in_rows, prev_staples, row_start);
    }
    
    // 新增：定期pruning
    if (site % PRUNING_FREQUENCY == 0 && site > 0) {
        pruneNodes(site);
    }
    
    // 新增：監控單個site的處理時間
    auto site_end_time = std::chrono::high_resolution_clock::now();
    auto site_duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        site_end_time - site_start_time);
    
    if (site_duration.count() > 5000) {  // 超過5秒的site
        Logger::log(Logger::WARNING, "Site " + std::to_string(site) + 
                    " took " + std::to_string(site_duration.count()) + 
                    "ms with " + std::to_string(all_nodes.size()) + " nodes");
        
        // 強制aggressive pruning
        if (all_nodes.size() > 500) {
            pruneNodes(site);
        }
    }
    
    // 更新memory tracking
    max_nodes_in_memory = std::max(max_nodes_in_memory, all_nodes.size());
}


void DPSolver::pruneNodes(int current_site) {
    // 更保守的pruning參數
    const int MAX_NODES = 5000;  // 增加最大節點數
    const int BENEFIT_THRESHOLD = 8;  // 增加threshold，更保守的pruning
    
    if (all_nodes.size() <= MAX_NODES) {
        return;
    }
    
    Logger::log(Logger::WARNING, "Pruning nodes at site " + std::to_string(current_site) + 
               ", current count: " + std::to_string(all_nodes.size()));
    
    // 找到當前最佳收益
    int best_benefit = -1000000;
    for (DPNode* node : all_nodes) {
        for (int i = 0; i < DPNode::NUM_CASES; i++) {
            best_benefit = std::max(best_benefit, node->benefit[i]);
        }
    }
    
    // 更保守的pruning：只移除收益顯著差的nodes
    std::vector<DPNode*> nodes_to_keep;
    int pruned_count = 0;
    
    for (DPNode* node : all_nodes) {
        int node_best_benefit = -1000000;
        for (int i = 0; i < DPNode::NUM_CASES; i++) {
            node_best_benefit = std::max(node_best_benefit, node->benefit[i]);
        }
        
        // 只有當收益顯著低於最佳時才移除
        if (node_best_benefit < best_benefit - BENEFIT_THRESHOLD && nodes_to_keep.size() > 100) {
            // 從lookup table移除
            CompactState compact;
            compact.site = node->site;
            for (int i = 0; i < 3; i++) {
                compact.a[i] = node->s[i] - initial_s[i];
                compact.b[i] = node->l[i];
            }
            node_lookup.erase(compact);
            delete node;
            pruned_count++;
        } else {
            nodes_to_keep.push_back(node);
        }
    }
    
    // 確保至少保留一些nodes
    if (nodes_to_keep.size() < 50) {
        Logger::log(Logger::ERROR, "Pruning would remove too many nodes, keeping all");
        return;
    }
    
    all_nodes = std::move(nodes_to_keep);
    
    Logger::log(Logger::WARNING, "Pruned " + std::to_string(pruned_count) + 
               " nodes, remaining: " + std::to_string(all_nodes.size()));
}

/**
 * @brief Generate extensions using binary decision tree
 */
void DPSolver::generateExtensions(DPNode* node,
                                const std::vector<std::vector<Cell*>>& cells_in_rows,
                                const std::vector<Staple>& prev_staples,
                                int row_start) {
    
    Logger::log(Logger::INFO, "=== GENERATING EXTENSIONS ===");
    Logger::log(Logger::INFO, "Site: " + std::to_string(node->site));
    Logger::log(Logger::INFO, "State: (" + 
               std::to_string(node->s[0]) + "," + std::to_string(node->l[0]) + "," +
               std::to_string(node->s[1]) + "," + std::to_string(node->l[1]) + "," +
               std::to_string(node->s[2]) + "," + std::to_string(node->l[2]) + ")");
    
    // 新增：顯示每個row的cell width information
    for (int row = 0; row < 3; row++) {
        Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + " detailed info:");
        if (node->s[row] > 0 && node->s[row] <= static_cast<int>(cells_in_rows[row].size())) {
            Cell* last_placed = cells_in_rows[row][node->s[row] - 1];
            const CellType& last_cell_type = cell_types[last_placed->type_index];
            int last_width_sites = last_cell_type.width / chip_info.site_width;
            
            Logger::log(Logger::DEBUG, "  Last placed cell " + std::to_string(last_placed->cell_index) + 
                       " width: " + std::to_string(last_cell_type.width) + 
                       " (" + std::to_string(last_width_sites) + " sites)");
            Logger::log(Logger::DEBUG, "  Current l=" + std::to_string(node->l[row]) + 
                       ", cell_width=" + std::to_string(last_width_sites) + 
                       ", cleared=" + (last_width_sites <= node->l[row] ? "YES" : "NO"));
        }
    }
    
    int valid_extensions = 0;
    
    // Try 8 extension combinations
    for (int mask = 0; mask < 8; mask++) {
        bool place_r1 = (mask & 0x4) != 0;
        bool place_r2 = (mask & 0x2) != 0;
        bool place_r3 = (mask & 0x1) != 0;
        
        Logger::log(Logger::DEBUG, "Testing extension " + std::to_string(mask) + ": " +
                   (place_r1 ? "P" : "D") + (place_r2 ? "P" : "D") + (place_r3 ? "P" : "D"));
        
        // Check extension valid
        bool valid = true;
        for (int row = 0; row < 3; row++) {
            bool place_row = (row == 0) ? place_r1 : (row == 1) ? place_r2 : place_r3;
            bool can_do = place_row ? 
                canPlaceNextCell(row, node->s[row], node->l[row], node->site, cells_in_rows[row]) :
                canDeferNextCell(row, node->s[row], node->l[row], node->site, cells_in_rows[row]);
            
            if (!can_do) {
                valid = false;
                Logger::log(Logger::DEBUG, "  Row " + std::to_string(row) + ": FAIL");
                break;
            } else {
                Logger::log(Logger::DEBUG, "  Row " + std::to_string(row) + ": OK");
            }
        }
        
        if (valid) {
            valid_extensions++;
            Logger::log(Logger::DEBUG, "Extension " + std::to_string(mask) + ": VALID");
            
            // 修正：正確計算新狀態的l值
            int new_s[3], new_l[3];
            for (int row = 0; row < 3; row++) {
                bool place_row = (row == 0) ? place_r1 : (row == 1) ? place_r2 : place_r3;
                
                if (place_row) {
                    new_s[row] = node->s[row] + 1;
                    
                    // 關鍵修正：計算實際的l值
                    if (new_s[row] <= static_cast<int>(cells_in_rows[row].size())) {
                        Cell* placed_cell = cells_in_rows[row][new_s[row] - 1];
                        const CellType& cell_type = cell_types[placed_cell->type_index];
                        int cell_width_sites = cell_type.width / chip_info.site_width;
                        
                        // l是從下一個site到cell左邊界的距離
                        // 如果cell在site i+1被placed，那麼在site i+1，距離是cell_width_sites
                        new_l[row] = cell_width_sites;
                        
                        Logger::log(Logger::DEBUG, "  Placed cell " + std::to_string(placed_cell->cell_index) + 
                                   " width=" + std::to_string(cell_width_sites) + " sites, new_l=" + std::to_string(new_l[row]));
                    } else {
                        new_l[row] = 1;
                    }
                } else {
                    new_s[row] = node->s[row];
                    new_l[row] = node->l[row] + 1;
                }
            }
            
            // 創建新節點
            DPNode* target = createOrFindNode(node, cells_in_rows, node->site + 1, 
                                             new_s[0], new_l[0],
                                             new_s[1], new_l[1], 
                                             new_s[2], new_l[2]);
            
            if (target != nullptr) {
                updateBenefit(node, target, cells_in_rows, prev_staples, row_start);
            }
        } else {
            Logger::log(Logger::DEBUG, "Extension " + std::to_string(mask) + ": INVALID");
        }
    }
    
    Logger::log(Logger::INFO, "Total valid extensions: " + std::to_string(valid_extensions));
    
    // 如果沒有valid extensions，顯示詳細的blocking原因
    if (valid_extensions == 0) {
        Logger::log(Logger::WARNING, "ALGORITHM STUCK - No valid extensions at site " + std::to_string(node->site));
        Logger::log(Logger::WARNING, "Debugging information:");
        
        for (int row = 0; row < 3; row++) {
            if (node->s[row] > 0) {
                Cell* last_placed = cells_in_rows[row][node->s[row] - 1];
                const CellType& cell_type = cell_types[last_placed->type_index];
                int cell_width_sites = cell_type.width / chip_info.site_width;
                
                Logger::log(Logger::WARNING, "Row " + std::to_string(row) + 
                           ": Last placed cell " + std::to_string(last_placed->cell_index) + 
                           " spans " + std::to_string(cell_width_sites) + " sites" +
                           ", current l=" + std::to_string(node->l[row]) + 
                           ", blocking=" + (cell_width_sites > node->l[row] ? "YES" : "NO"));
            }
        }
    }
    
    Logger::log(Logger::DEBUG, "=== END EXTENSION GENERATION ===");
}


/**
 * @brief Check if we can defer placing the next cell
 */
bool DPSolver::canDeferNextCell(int row, int s_j, int l_j, int site,
                               const std::vector<Cell*>& row_cells) {
    
    // 檢查最後放置的cell是否已經清空current site
    if (s_j > 0 && s_j <= static_cast<int>(row_cells.size())) {
        Cell* last_placed = row_cells[s_j - 1];
        const CellType& cell_type = cell_types[last_placed->type_index];
        int cell_width_sites = cell_type.width / chip_info.site_width;
        
        Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + " defer check:");
        Logger::log(Logger::DEBUG, "  Last placed cell " + std::to_string(last_placed->cell_index) + 
                   " width: " + std::to_string(cell_width_sites) + " sites");
        Logger::log(Logger::DEBUG, "  Current l_j: " + std::to_string(l_j));
        Logger::log(Logger::DEBUG, "  Clearance check: " + std::to_string(cell_width_sites) + 
                   " <= " + std::to_string(l_j) + " = " + (cell_width_sites <= l_j ? "OK" : "FAIL"));
        
        // Paper條件：w(c_j_sj) <= l_j (cell已經清空site)
        if (cell_width_sites > l_j) {
            Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + 
                       ": Cannot defer - last placed cell still occupies site");
            return false;
        }
    }
    
    // 如果沒有更多cells要處理，可以defer
    Cell* next_cell = getNextCellToPlace(row, s_j, row_cells);
    if (next_cell == nullptr) {
        Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + 
                   ": Can defer - no more cells");
        return true;
    }
    
    // 檢查下一個cell是否可以在其displacement限制內被放置
    int initial_site = next_cell->initial_x / chip_info.site_width;
    int max_displacement_sites = next_cell->max_displacement / chip_info.site_width;
    int latest_site = initial_site + max_displacement_sites;
    
    bool can_defer = (site < latest_site);
    
    Logger::log(Logger::DEBUG, "  Next cell " + std::to_string(next_cell->cell_index) + 
               " latest placement site: " + std::to_string(latest_site) + 
               ", current site: " + std::to_string(site) + 
               ", can defer: " + (can_defer ? "YES" : "NO"));
    
    return can_defer;
}


/**
 * @brief Check if we can place the next cell at site + 1
 */
bool DPSolver::canPlaceNextCell(int row, int s_j, int l_j, int site,
                               const std::vector<Cell*>& row_cells) {
    // Check if the last placed cell has cleared the current site
    Cell* last_placed = getLastPlacedCell(row, s_j, row_cells);
    if (last_placed != nullptr) {
        const CellType& cell_type = cell_types[last_placed->type_index];
        int cell_width_sites = cell_type.width / chip_info.site_width;
        
        // Paper condition: w(c_j_sj) <= l_j (cell has cleared)
        if (cell_width_sites > l_j) {
            Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + 
                       ": Cannot place - last cell still occupies site");
            return false;
        }
    }
    
    // Get the next cell to place
    Cell* next_cell = getNextCellToPlace(row, s_j, row_cells);
    if (next_cell == nullptr) {
        Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + 
                   ": Cannot place - no more cells");
        return false;
    }
    
    // Check if we can place it at site + 1
    int target_site = site + 1;
    int target_x = target_site * chip_info.site_width;
    int displacement = std::abs(target_x - next_cell->initial_x);
    
    bool can_place = displacement <= next_cell->max_displacement;
    
    Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + " place check:");
    Logger::log(Logger::DEBUG, "  Next cell " + std::to_string(next_cell->cell_index) + 
               " initial_x=" + std::to_string(next_cell->initial_x));
    Logger::log(Logger::DEBUG, "  Target site=" + std::to_string(target_site) + 
               " (x=" + std::to_string(target_x) + ")");
    Logger::log(Logger::DEBUG, "  Displacement=" + std::to_string(displacement) + 
               ", Max=" + std::to_string(next_cell->max_displacement));
    std::string result_str = can_place ? "OK" : "FAIL";
    Logger::log(Logger::DEBUG, "  Can place: " + result_str);
    
    return can_place;
}


/**
 * @brief Create or find node with given state
 */
DPNode* DPSolver::createOrFindNode(DPNode* node, const std::vector<std::vector<Cell*>>& cells_in_rows, 
                                int site, int s1, int l1, int s2, int l2, int s3, int l3) {
    // Create compact state for lookup
    CompactState compact;
    compact.site = site;
    compact.a[0] = s1 - initial_s[0];  // Cell offset
    compact.a[1] = s2 - initial_s[1];
    compact.a[2] = s3 - initial_s[2];

    // Calculate actual cell displacements (b values)
    for (int row = 0; row < 3; row++) {
        if (node->s[row] > 0 && node->s[row] <= static_cast<int>(cells_in_rows[row].size())) {
            Cell* cell = cells_in_rows[row][node->s[row] - 1];
            int current_pos = site * chip_info.site_width - node->l[row] * chip_info.site_width;
            compact.b[row] = current_pos - cell->initial_x;
        } else {
            compact.b[row] = 0;
        }
    }
    
    // Check if node already exists
    auto it = node_lookup.find(compact);
    if (it != node_lookup.end()) {
        return it->second;
    }
    
    // Create new node
    DPNode* new_node = new DPNode(site, s1, l1, s2, l2, s3, l3);
    node_lookup[compact] = new_node;
    node_queue.push(new_node);
    all_nodes.push_back(new_node);
    nodes_created++;
    
    return new_node;
}

/**
 * @brief Update benefit for target node based on source node
 */
void DPSolver::updateBenefit(DPNode* source, DPNode* target,
                           const std::vector<std::vector<Cell*>>& cells_in_rows,
                           const std::vector<Staple>& prev_staples,
                           int row_start) {
    
    // Get valid staple cases for target site
    std::vector<int> valid_cases = getValidStapleCases(
        target->site, target->s, target->l, cells_in_rows);
    
    // Try all combinations of source and target cases
    for (int src_case = 0; src_case < DPNode::NUM_CASES; src_case++) {
        // Skip if source case has invalid benefit
        if (source->benefit[src_case] < -900000) continue;
        
        for (int tgt_case : valid_cases) {
            // Check anti-parallel line-ends constraint
            if (hasAntiParallelViolation(target->site, tgt_case, src_case, 
                                        prev_staples, row_start)) {
                continue;
            }
            
            // Calculate staple benefit for this case
            int staple_benefit = 0;
            int vdd_added = 0, vss_added = 0;
            
            switch (tgt_case) {
                case CASE_1_NO_STAPLE:
                    // No staple, might apply balance factor
                    if (params.balance_factor > 0 && 
                        source->vdd_staples[src_case] != source->vss_staples[src_case]) {
                        staple_benefit = static_cast<int>(params.balance_factor * 100);
                    }
                    break;
                    
                case CASE_2_R1_R2:
                    staple_benefit = 100;  // Scale for integer math
                    if (isVDDRow(row_start)) vdd_added = 1;
                    else vss_added = 1;
                    break;
                    
                case CASE_3_R2_R3:
                    staple_benefit = 100;
                    if (isVDDRow(row_start + 1)) vdd_added = 1;
                    else vss_added = 1;
                    break;
                    
                case CASE_4_BOTH:
                    staple_benefit = 200;
                    if (isVDDRow(row_start)) vdd_added++;
                    else vss_added++;
                    if (isVDDRow(row_start + 1)) vdd_added++;
                    else vss_added++;
                    break;
                    
                case CASE_5_SPECIAL:
                    // Special case handling based on configuration
                    staple_benefit = 150;  // Example value
                    break;
            }
            
            // Apply balance adjustment
            int total_vdd = source->vdd_staples[src_case] + vdd_added;
            int total_vss = source->vss_staples[src_case] + vss_added;
            staple_benefit = calculateStapleBenefit(tgt_case, total_vdd, total_vss, row_start);
            
            // Update if better
            int new_benefit = source->benefit[src_case] + staple_benefit;
            if (new_benefit > target->benefit[tgt_case]) {
                target->benefit[tgt_case] = new_benefit;
                target->prev_node[tgt_case] = source;
                target->prev_case[tgt_case] = src_case;
                target->vdd_staples[tgt_case] = total_vdd;
                target->vss_staples[tgt_case] = total_vss;
            }
        }
    }
}

/**
 * @brief Check if a staple can be inserted (no pin obstruction)
 */
bool DPSolver::canInsertStaple(int site, int between_rows, 
                              int s[3], int l[3],
                              const std::vector<std::vector<Cell*>>& cells_in_rows) {
    int row1 = between_rows;
    int row2 = between_rows + 1;
    
    if (row1 < 0 || row2 >= 3) return false;
    
    // Check row1
    if (s[row1] > 0 && s[row1] <= static_cast<int>(cells_in_rows[row1].size())) {
        Cell* cell = cells_in_rows[row1][s[row1] - 1];
        const CellType& cell_type = cell_types[cell->type_index];
        int cell_width_in_sites = cell_type.width / chip_info.site_width;
        
        if (l[row1] < cell_width_in_sites) {
            // Cell overlaps site, check for pin
            int site_in_cell = cell_width_in_sites - l[row1];
            if (site_in_cell >= 0 && site_in_cell < cell_width_in_sites) {
                if (cell_type.hasPinAt(site_in_cell, cell->is_flipped)) {
                    return false;
                }
            }
        }
    }
    
    // Check row2 (similar logic)
    if (s[row2] > 0 && s[row2] <= static_cast<int>(cells_in_rows[row2].size())) {
        Cell* cell = cells_in_rows[row2][s[row2] - 1];
        const CellType& cell_type = cell_types[cell->type_index];
        int cell_width_in_sites = cell_type.width / chip_info.site_width;
        
        if (l[row2] < cell_width_in_sites) {
            int site_in_cell = cell_width_in_sites - l[row2];
            if (site_in_cell >= 0 && site_in_cell < cell_width_in_sites) {
                if (cell_type.hasPinAt(site_in_cell, cell->is_flipped)) {
                    return false;
                }
            }
        }
    }
    
    return true;
}

/**
 * @brief Get valid staple cases at current site
 */
std::vector<int> DPSolver::getValidStapleCases(int site, int s[3], int l[3],
                                              const std::vector<std::vector<Cell*>>& cells_in_rows) {
    std::vector<int> valid_cases;
    
    // Case 1: No staple (always valid)
    valid_cases.push_back(CASE_1_NO_STAPLE);
    
    // Case 2: R1-R2 staple
    if (canInsertStaple(site, 0, s, l, cells_in_rows)) {
        valid_cases.push_back(CASE_2_R1_R2);
    }
    
    // Case 3: R2-R3 staple
    if (canInsertStaple(site, 1, s, l, cells_in_rows)) {
        valid_cases.push_back(CASE_3_R2_R3);
    }
    
    // Case 4: Both staples
    if (canInsertStaple(site, 0, s, l, cells_in_rows) &&
        canInsertStaple(site, 1, s, l, cells_in_rows)) {
        valid_cases.push_back(CASE_4_BOTH);
    }
    
    // Case 5: Special configuration (requires specific conditions)
    // This would depend on the specific implementation requirements
    
    return valid_cases;
}

/**
 * @brief Check for anti-parallel line-ends violation
 */
bool DPSolver::hasAntiParallelViolation(int site, int staple_case, int prev_case,
                                       const std::vector<Staple>& prev_staples,
                                       int row_start) {
    // Check violations within current three rows (as shown in Figure 8 of paper)
    
    // Case-specific checks based on the paper
    if (staple_case == CASE_2_R1_R2 && prev_case == CASE_3_R2_R3) {
        // Potential staggering pattern
        return true;
    }
    
    if (staple_case == CASE_3_R2_R3 && prev_case == CASE_2_R1_R2) {
        // Potential staggering pattern
        return true;
    }
    
    // Check violations with previous round staples (Figure 9 in paper)
    int site_x = site * chip_info.site_width;
    
    for (const Staple& prev_staple : prev_staples) {
        // Check if previous staple creates anti-parallel line-ends
        if (std::abs(prev_staple.x - site_x) == chip_info.site_width) {
            // Adjacent sites, check for violation
            int prev_row = prev_staple.y / chip_info.row_height;
            
            if (staple_case == CASE_2_R1_R2 && prev_row == row_start - 1) {
                return true;
            }
            
            if (staple_case == CASE_3_R2_R3 && prev_row == row_start + 3) {
                return true;
            }
        }
    }
    
    return false;
}

/**
 * @brief Calculate staple benefit with balance factor
 */
int DPSolver::calculateStapleBenefit(int staple_case, int vdd_count, int vss_count,
                                   int row_start) {
    int base_benefit = 0;
    
    switch (staple_case) {
        case CASE_1_NO_STAPLE:
            base_benefit = 0;
            break;
        case CASE_2_R1_R2:
        case CASE_3_R2_R3:
            base_benefit = 100;
            break;
        case CASE_4_BOTH:
            base_benefit = 200;
            break;
        case CASE_5_SPECIAL:
            base_benefit = 150;
            break;
    }
    
    // Apply balance factor as described in Section 2 of the paper
    if (vdd_count > 0 && vss_count > 0) {
        double ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                      std::min(vdd_count, vss_count);
        
        if (ratio > 1.1) {
            // Apply balance penalty/bonus
            if ((vdd_count > vss_count && staple_case == CASE_3_R2_R3 && !isVDDRow(row_start + 1)) ||
                (vss_count > vdd_count && staple_case == CASE_2_R1_R2 && isVDDRow(row_start))) {
                // This staple helps balance
                base_benefit += static_cast<int>(params.balance_factor * 100);
            }
        }
    }
    
    return base_benefit;
}

/**
 * @brief Check if row is VDD or VSS
 */
bool DPSolver::isVDDRow(int row_idx) const {
    // Even rows are VDD, odd rows are VSS
    return (row_idx % 2 == 0);
}

/**
 * @brief Backtrack from best node to recover solution
 */
std::vector<Staple> DPSolver::backtrackSolution(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start) {
    
    std::vector<Staple> solution;
    
    if (!best_final_node || best_final_case < 0) {
        Logger::log(Logger::WARNING, "No valid solution to backtrack");
        return solution;
    }
    
    Logger::log(Logger::INFO, "Backtracking from final node with case " + 
                std::to_string(best_final_case));
    
    // Trace back through the DAG
    DPNode* current = best_final_node;
    int current_case = best_final_case;
    
    while (current != nullptr && current->prev_node[current_case] != nullptr) {
        // Insert staples based on the case
        int site_x = current->site * chip_info.site_width;
        
        switch (current_case) {
            case CASE_2_R1_R2: {
                int y = chip_info.getRowY(row_start + 1);
                solution.push_back(Staple(site_x, y, isVDDRow(row_start)));
                break;
            }
            
            case CASE_3_R2_R3: {
                int y = chip_info.getRowY(row_start + 2);
                solution.push_back(Staple(site_x, y, isVDDRow(row_start + 1)));
                break;
            }
            
            case CASE_4_BOTH: {
                int y1 = chip_info.getRowY(row_start + 1);
                int y2 = chip_info.getRowY(row_start + 2);
                solution.push_back(Staple(site_x, y1, isVDDRow(row_start)));
                solution.push_back(Staple(site_x, y2, isVDDRow(row_start + 1)));
                break;
            }
            
            case CASE_5_SPECIAL:
                // Handle special configuration
                break;
        }
        
        // Move to previous node
        DPNode* prev = current->prev_node[current_case];
        int prev_case = current->prev_case[current_case];
        current = prev;
        current_case = prev_case;
    }
    
    // Reverse to get correct order
    std::reverse(solution.begin(), solution.end());
    
    Logger::log(Logger::INFO, "Backtracking complete, " + 
                std::to_string(solution.size()) + " staples inserted");
    
    return solution;
}

/**
 * @brief Apply cell placement updates from backtracking
 */
void DPSolver::applyCellPlacements(const std::vector<std::vector<Cell*>>& cells_in_rows) {
    // This would update the current_x positions of cells based on the final placement
    // For now, keeping cells at their initial positions
    // Cell flipping would be handled as a post-processing step
    
    Logger::log(Logger::INFO, "Cell placement updates applied");
}

/**
 * @brief Compute initial s values for compact encoding
 */
void DPSolver::computeInitialS(const std::vector<std::vector<Cell*>>& cells_in_rows) {
    Logger::log(Logger::DEBUG, "Computing initial s values...");
    // For each row, find cells that cross or are to the left of site 0
    for (int row = 0; row < 3; row++) {
        initial_s[row] = 0;  // Initial no cell be processed
        
        Logger::log(Logger::DEBUG, "Row " + std::to_string(row) + ":");
        Logger::log(Logger::DEBUG, "  Total cells: " + std::to_string(cells_in_rows[row].size()));
        
        if (!cells_in_rows[row].empty()) {
            Cell* first_cell = cells_in_rows[row][0];
            int first_site = first_cell->initial_x / chip_info.site_width;
            Logger::log(Logger::DEBUG, "  First cell at site: " + std::to_string(first_site));
        }
        
        Logger::log(Logger::DEBUG, "  Initial s[" + std::to_string(row) + "] = " + std::to_string(initial_s[row]));
    }
}



/**
 * @brief Clean up allocated memory
 */
void DPSolver::cleanup() {
    Logger::log(Logger::DEBUG, "Starting DPSolver cleanup");
    
    size_t initial_nodes = all_nodes.size();
    size_t initial_lookup = node_lookup.size();
    
    // 清理所有節點
    for (DPNode* node : all_nodes) {
        if (node != nullptr) {
            delete node;
        }
    }
    all_nodes.clear();
    node_lookup.clear();
    
    // 清理隊列
    while (!node_queue.empty()) {
        node_queue.pop();
    }
    
    // Reset所有狀態變數
    nodes_created = 0;
    max_nodes_in_memory = 0;
    best_final_node = nullptr;
    best_final_case = -1;
    
    // Reset initial_s values
    for (int i = 0; i < 3; i++) {
        initial_s[i] = 0;
    }
    
    Logger::log(Logger::DEBUG, "DPSolver cleanup completed: deleted " + 
               std::to_string(initial_nodes) + " nodes, " + 
               std::to_string(initial_lookup) + " lookup entries");
}

/**
 * @brief Force complete cleanup of all data structures
 */
void DPSolver::forceCompleteCleanup() {
    Logger::log(Logger::DEBUG, "Force complete cleanup starting...");
    
    // 清理所有節點
    for (DPNode* node : all_nodes) {
        if (node != nullptr) {
            delete node;
        }
    }
    all_nodes.clear();
    
    // 清理lookup table
    node_lookup.clear();
    
    // 清理隊列
    while (!node_queue.empty()) {
        node_queue.pop();
    }
    
    // Reset所有狀態
    nodes_created = 0;
    max_nodes_in_memory = 0;
    best_final_node = nullptr;
    best_final_case = -1;
    
    // Reset initial_s
    for (int i = 0; i < 3; i++) {
        initial_s[i] = 0;
    }
    
    Logger::log(Logger::DEBUG, "Force complete cleanup finished");
}

/**
 * @brief Process nodes at site with timeout protection
 */
void DPSolver::processNodesAtSiteWithTimeout(int site,
                                           const std::vector<std::vector<Cell*>>& cells_in_rows,
                                           const std::vector<Staple>& prev_staples,
                                           int row_start,
                                           int max_seconds) {
    
    auto site_start = std::chrono::high_resolution_clock::now();
    
    // 收集當前site的所有nodes
    std::vector<DPNode*> nodes_at_site;
    size_t queue_size = node_queue.size();
    
    for (size_t i = 0; i < queue_size; i++) {
        DPNode* node = node_queue.front();
        node_queue.pop();
        
        if (node->site == site) {
            nodes_at_site.push_back(node);
        } else if (node->site > site) {
            node_queue.push(node); // Put back for later
        }
        // Nodes with site < current site are discarded (should not happen)
    }
    
    Logger::log(Logger::DEBUG, "Processing " + std::to_string(nodes_at_site.size()) + 
               " nodes at site " + std::to_string(site));
    
    // Process each node with timeout check
    for (size_t i = 0; i < nodes_at_site.size(); i++) {
        // Check timeout every 100 nodes
        if (i % 100 == 0) {
            auto current = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(current - site_start);
            
            if (elapsed.count() > max_seconds) {
                Logger::log(Logger::WARNING, "Site " + std::to_string(site) + 
                           " timeout after processing " + std::to_string(i) + " nodes");
                break;
            }
        }
        
        generateExtensions(nodes_at_site[i], cells_in_rows, prev_staples, row_start);
    }
    
    // Pruning after processing
    if (site % PRUNING_FREQUENCY == 0 && site > 0) {
        pruneNodes(site);
    }
}

/**
 * @brief Aggressive pruning to reduce memory usage
 */
void DPSolver::aggressivePruning(int current_site) {
    Logger::log(Logger::WARNING, "Applying aggressive pruning at site " + std::to_string(current_site));
    
    const int AGGRESSIVE_MAX_NODES = 500;
    const int AGGRESSIVE_THRESHOLD = 2;
    
    if (all_nodes.size() <= AGGRESSIVE_MAX_NODES) {
        return;
    }
    
    // 找最佳收益
    int best_benefit = -1000000;
    for (DPNode* node : all_nodes) {
        for (int i = 0; i < DPNode::NUM_CASES; i++) {
            best_benefit = std::max(best_benefit, node->benefit[i]);
        }
    }
    
    // 只保留最好的nodes
    std::vector<DPNode*> nodes_to_keep;
    int pruned_count = 0;
    
    for (DPNode* node : all_nodes) {
        int node_best = -1000000;
        for (int i = 0; i < DPNode::NUM_CASES; i++) {
            node_best = std::max(node_best, node->benefit[i]);
        }
        
        if (node_best >= best_benefit - AGGRESSIVE_THRESHOLD && 
            nodes_to_keep.size() < AGGRESSIVE_MAX_NODES) {
            nodes_to_keep.push_back(node);
        } else {
            // Remove from lookup table
            CompactState compact;
            compact.site = node->site;
            for (int i = 0; i < 3; i++) {
                compact.a[i] = node->s[i] - initial_s[i];
                compact.b[i] = node->l[i];
            }
            node_lookup.erase(compact);
            delete node;
            pruned_count++;
        }
    }
    
    all_nodes = std::move(nodes_to_keep);
    
    Logger::log(Logger::WARNING, "Aggressive pruning: removed " + std::to_string(pruned_count) + 
               " nodes, kept " + std::to_string(all_nodes.size()));
}

/**
 * @brief Extract optimal solution from DAG
 */
std::vector<Staple> DPSolver::extractOptimalSolution(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start) {
    
    // Find best final node
    int max_benefit = -1000000;
    for (DPNode* node : all_nodes) {
        if (node->site == chip_info.total_sites) {
            for (int c = 0; c < DPNode::NUM_CASES; c++) {
                if (node->benefit[c] > max_benefit) {
                    max_benefit = node->benefit[c];
                    best_final_node = node;
                    best_final_case = c;
                }
            }
        }
    }
    
    if (best_final_node != nullptr) {
        Logger::log(Logger::INFO, "Optimal solution found with benefit: " + std::to_string(max_benefit));
        return backtrackSolution(cells_in_rows, row_start);
    } else {
        Logger::log(Logger::WARNING, "No valid final solution found");
        return generateFallbackSolution(cells_in_rows, row_start, {});
    }
}

/**
 * @brief Generate fallback solution when main algorithm fails
 */
std::vector<Staple> DPSolver::generateFallbackSolution(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start,
    const std::vector<Staple>& prev_staples) {
    
    Logger::log(Logger::WARNING, "Generating fallback solution for rows " + 
               std::to_string(row_start) + "-" + std::to_string(row_start+2));
    
    std::vector<Staple> fallback_solution;
    
    // Simple greedy approach: insert staples where possible
    for (int site = 5; site < chip_info.total_sites - 5; site += 8) { // Every 8 sites
        int x = site * chip_info.site_width;
        
        // Try to insert staples between adjacent rows
        // R1-R2 staple
        if (row_start + 1 < chip_info.num_rows) {
            int y = chip_info.getRowY(row_start + 1);
            bool is_vdd = isVDDRow(row_start);
            fallback_solution.push_back(Staple(x, y, is_vdd));
        }
        
        // R2-R3 staple  
        if (row_start + 2 < chip_info.num_rows) {
            int y = chip_info.getRowY(row_start + 2);
            bool is_vdd = isVDDRow(row_start + 1);
            fallback_solution.push_back(Staple(x, y, is_vdd));
        }
        
        // Stop if we have enough staples
        if (fallback_solution.size() >= 50) break;
    }
    
    Logger::log(Logger::INFO, "Fallback solution generated: " + 
               std::to_string(fallback_solution.size()) + " staples");
    
    return fallback_solution;
}