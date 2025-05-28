/**
 * @file dp_solver.cpp
 * @brief Enhanced implementation of DAG-based dynamic programming solver for MATRO
 */

#include "dp_solver.hpp"
#include "../Logger.hpp"
#include "../memory_usage.hpp"
#include <algorithm>
#include <limits>
#include <cmath>

/**
 * @brief Constructor
 */
DPSolver::DPSolver(const ChipInfo& chip_info, 
                   const std::vector<CellType>& cell_types,
                   const AlgorithmParams& params)
    : chip_info(chip_info), cell_types(cell_types), params(params),
      best_final_node(nullptr), best_final_case(-1), nodes_created(0) {
    
    Logger::log(Logger::INFO, "Enhanced DPSolver initialized");
    Logger::log(Logger::INFO, "  Enhanced compact state with displacement tracking");
    Logger::log(Logger::INFO, "  Intelligent pruning: BEAM_WIDTH=" + std::to_string(BEAM_WIDTH));
}

/**
 * @brief Destructor
 */
DPSolver::~DPSolver() {
    cleanup();
}

/**
 * @brief Enhanced MATRO implementation with proper constraint checking
 */
std::vector<Staple> DPSolver::solveTripleRow(
    const std::vector<std::vector<Cell*>>& cells_in_rows,
    int row_start,
    const std::vector<Staple>& prev_staples) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    size_t start_memory = getCurrentMemoryUsage();
    
    std::cout << "=== ENHANCED TRIPLE-ROW " << row_start << "-" << (row_start+2) 
              << " (Memory: " << start_memory << "MB) ===" << std::endl;
    
    // Clean up previous state
    cleanup();
    
    // Initialize cell positions for displacement tracking
    for (int row = 0; row < 3; row++) {
        initial_s[row] = 0;
        if (row < static_cast<int>(cells_in_rows.size())) {
            for (size_t i = 0; i < cells_in_rows[row].size() && i < 1000; i++) {
                initial_cell_x[row][i] = cells_in_rows[row][i]->initial_x;
            }
        }
    }
    
    try {
        // Create enhanced source node
        EnhancedDPNode* source = new EnhancedDPNode(0, 0, 0, 0, 0, 0, 0);
        source->benefit[CASE_1_NO_STAPLE] = 0;
        
        std::queue<EnhancedDPNode*> Q;
        Q.push(source);
        all_nodes.push_back(source);
        
        // Create enhanced compact state
        EnhancedCompactState compact;
        compact.site = 0;
        for (int i = 0; i < 3; i++) {
            compact.a[i] = 0;
            compact.b[i] = 0;
            compact.flipped[i] = false;
        }
        node_lookup[compact] = source;
        
        std::cout << "Enhanced source created with displacement tracking" << std::endl;
        
        // Main DAG construction loop
        int processed_nodes = 0;
        int last_pruning_site = -1;
        
        while (!Q.empty()) {
            auto current_time = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
                current_time - start_time);
            
            if (elapsed.count() > 45) {
                std::cout << "Enhanced solver timeout at " << elapsed.count() << "s" << std::endl;
                break;
            }
            
            EnhancedDPNode* u = Q.front();
            Q.pop();
            processed_nodes++;
            
            int i = u->site;
            
            if (i >= chip_info.total_sites) {
                continue;
            }
            
            // Intelligent pruning every PRUNING_FREQUENCY sites
            if (i > last_pruning_site + PRUNING_FREQUENCY) {
                pruneNodesIntelligently(i);
                last_pruning_site = i;
            }
            
            // Global node limit check
            if (all_nodes.size() > MAX_TOTAL_NODES) {
                std::cout << "Enhanced global node limit reached (" << all_nodes.size() 
                          << "), applying aggressive pruning" << std::endl;
                pruneNodesIntelligently(i);
            }
            
            // Generate extensions (including flipping possibilities)
            int valid_extensions = 0;
            
            for (int mask = 0; mask < 8; mask++) {
                bool place_r1 = (mask & 4) != 0;
                bool place_r2 = (mask & 2) != 0; 
                bool place_r3 = (mask & 1) != 0;
                
                // For each placement mask, try different flipping combinations
                int flip_combinations = 1;
                for (int row = 0; row < 3; row++) {
                    if ((row == 0 && place_r1) || (row == 1 && place_r2) || (row == 2 && place_r3)) {
                        flip_combinations = 2; // Try both original and flipped
                        break;
                    }
                }
                
                for (int flip_mask = 0; flip_mask < flip_combinations; flip_mask++) {
                    bool flip_r1 = place_r1 && (flip_mask & 1);
                    bool flip_r2 = place_r2 && (flip_mask & 1);
                    bool flip_r3 = place_r3 && (flip_mask & 1);
                    
                    int new_s[3], new_l[3], new_disp[3];
                    bool new_flipped[3] = {flip_r1, flip_r2, flip_r3};
                    bool placement_mask[3];
                    
                    // Check extension legality with enhanced validation
                    if (isLegalExtension(u, i + 1, new_s, new_l, new_disp, new_flipped,
                                      placement_mask, cells_in_rows, prev_staples, row_start)) {
                        
                        valid_extensions++;
                        
                        // Create or find target node
                        EnhancedDPNode* v = createOrFindEnhancedNode(
                            i + 1, new_s[0], new_l[0], new_s[1], new_l[1], new_s[2], new_l[2],
                            new_disp[0], new_disp[1], new_disp[2],
                            new_flipped[0], new_flipped[1], new_flipped[2]);
                        
                        if (v != nullptr) {
                            Q.push(v);
                            // Enhanced benefit update with early constraint checking
                            updateBenefitEnhanced(u, v, cells_in_rows, prev_staples, row_start);
                        }
                    }
                }
            }
            
            if (processed_nodes % 1000 == 0) {
                std::cout << "Enhanced DP: processed " << processed_nodes 
                          << " nodes, site " << i << ", valid ext: " << valid_extensions << std::endl;
            }
        }
        
        std::cout << "Enhanced DAG completed. Nodes: " << all_nodes.size() 
                  << ", Processed: " << processed_nodes << std::endl;
        
        return extractBestSolutionEnhanced(cells_in_rows, row_start);
        
    } catch (const std::exception& e) {
        std::cout << "Enhanced solver exception: " << e.what() << std::endl;
        cleanup();
        return generateEnhancedFallback(cells_in_rows, row_start);
    }
}

/**
 * @brief Enhanced legality checking with integrated staggering validation
 */
bool DPSolver::isLegalExtension(EnhancedDPNode* source, int target_site,
                                   int new_s[3], int new_l[3], int new_disp[3], bool new_flipped[3],
                                   bool placement_mask[3], // 新增參數：明確指示每行是否放置單元
                                   const std::vector<std::vector<Cell*>>& cells_in_rows,
                                   const std::vector<Staple>& prev_staples, int row_start) {
    
    // Calculate new state based on placement decisions
    for (int row = 0; row < 3; row++) {
        if (row >= static_cast<int>(cells_in_rows.size())) {
            // No cells in this row
            new_s[row] = source->s[row];
            new_l[row] = source->l[row] + 1;
            new_disp[row] = source->displacement[row];
            new_flipped[row] = source->is_flipped[row];
            continue;
        }
        
        if (placement_mask[row]) {
            // === 放置單元的情況 ===
            
            // 1. 檢查是否有更多單元可以放置
            if (source->s[row] >= static_cast<int>(cells_in_rows[row].size())) {
                return false; // 沒有更多單元
            }
            
            // 2. 檢查最後放置的單元是否已清空當前 site
            if (source->s[row] > 0) {
                Cell* last_cell = cells_in_rows[row][source->s[row] - 1];
                const CellType& cell_type = cell_types[last_cell->type_index];
                int cell_width_sites = cell_type.width / chip_info.site_width;
                
                if (cell_width_sites > source->l[row]) {
                    return false; // 最後單元仍占據 site
                }
            }
            
            // 3. 檢查要放置的單元的位移約束
            Cell* next_cell = cells_in_rows[row][source->s[row]];
            int target_x = target_site * chip_info.site_width;
            int displacement = std::abs(target_x - next_cell->initial_x);
            
            if (displacement > next_cell->max_displacement) {
                return false; // 位移過大
            }
            
            // 4. 更新狀態：放置單元
            new_s[row] = source->s[row] + 1;
            const CellType& cell_type = cell_types[next_cell->type_index];
            int cell_width_sites = cell_type.width / chip_info.site_width;
            new_l[row] = cell_width_sites;
            new_disp[row] = displacement; // 正確的位移追蹤
            // new_flipped[row] 由調用者設定
            
        } else {
            // === 延遲單元的情況 ===
            
            // 1. 檢查最後放置的單元是否已清空當前 site
            if (source->s[row] > 0) {
                Cell* last_cell = cells_in_rows[row][source->s[row] - 1];
                const CellType& cell_type = cell_types[last_cell->type_index];
                int cell_width_sites = cell_type.width / chip_info.site_width;
                
                if (cell_width_sites > source->l[row]) {
                    return false; // 最後單元仍占據 site
                }
            }
            
            // 2. 檢查下一個單元是否仍可在約束內被放置
            if (source->s[row] < static_cast<int>(cells_in_rows[row].size())) {
                Cell* next_cell = cells_in_rows[row][source->s[row]];
                int initial_site = next_cell->initial_x / chip_info.site_width;
                int max_displacement_sites = next_cell->max_displacement / chip_info.site_width;
                int latest_site = initial_site + max_displacement_sites;
                
                if (target_site > latest_site) {
                    return false; // 會違反下一個單元的約束
                }
            }
            
            // 3. 更新狀態：延遲單元
            new_s[row] = source->s[row];
            new_l[row] = source->l[row] + 1;
            new_disp[row] = source->displacement[row];
            new_flipped[row] = source->is_flipped[row];
        }
    }
    
    return true; // 擴展合法
}

/**
 * @brief Enhanced benefit update with early staggering and overlap checking
 */
void DPSolver::updateBenefitEnhanced(EnhancedDPNode* source, EnhancedDPNode* target,
                                    const std::vector<std::vector<Cell*>>& cells_in_rows,
                                    const std::vector<Staple>& prev_staples, int row_start) {
    
    for (int src_case = 0; src_case < EnhancedDPNode::NUM_CASES; src_case++) {
        if (source->benefit[src_case] < -900000) continue;
        
        for (int tgt_case = 0; tgt_case < 5; tgt_case++) {
            EnhancedStapleCase case_type = static_cast<EnhancedStapleCase>(tgt_case);
            
            // Calculate staples for this case
            std::vector<Staple> new_staples = calculateStaplesForCase(case_type, target, row_start);
            
            // Early constraint checking
            bool case_valid = true;
            
            // Check for staggering violations
            if (hasStaggeringViolationEnhanced(new_staples, prev_staples, target->site, row_start)) {
                case_valid = false;
            }
            
            // Check for pin and staple overlaps
            for (const Staple& staple : new_staples) {
                if (hasOverlapWithPinsOrStaples(staple, cells_in_rows, prev_staples, row_start)) {
                    case_valid = false;
                    break;
                }
            }
            
            if (!case_valid) {
                continue; // Skip this case due to constraint violations
            }
            
            // Calculate benefits
            int base_benefit = calculateEnhancedCaseBenefit(case_type, 
                source->vdd_staples[src_case], source->vss_staples[src_case]);
            
            // Count VDD/VSS staples for this case
            int vdd_added = 0, vss_added = 0;
            for (const Staple& staple : new_staples) {
                if (staple.is_vdd) vdd_added++;
                else vss_added++;
            }
            
            // Balance adjustment
            int new_vdd = source->vdd_staples[src_case] + vdd_added;
            int new_vss = source->vss_staples[src_case] + vss_added;
            
            int balance_adjustment = 0;
            if (new_vdd > 0 && new_vss > 0) {
                double ratio = static_cast<double>(std::max(new_vdd, new_vss)) / 
                              static_cast<double>(std::min(new_vdd, new_vss));
                
                if (ratio > 1.1) {
                    // Penalty for violating balance constraint
                    balance_adjustment -= static_cast<int>((ratio - 1.1) * 300);
                    if (ratio > 1.5) {
                        continue; // Reject severely unbalanced solutions
                    }
                } else if (ratio <= 1.05) {
                    balance_adjustment += 150; // Bonus for excellent balance
                } else if (ratio <= 1.1) {
                    balance_adjustment += 75;  // Bonus for good balance
                }
            }
            
            // First-time bonus for creating missing type
            if (source->vdd_staples[src_case] == 0 && vdd_added > 0) {
                balance_adjustment += 200;
            }
            if (source->vss_staples[src_case] == 0 && vss_added > 0) {
                balance_adjustment += 200;
            }
            
            // Final benefit calculation
            int final_benefit = base_benefit + balance_adjustment;
            int new_total_benefit = source->benefit[src_case] + final_benefit;
            
            // Update if better
            if (new_total_benefit > target->benefit[tgt_case]) {
                target->benefit[tgt_case] = new_total_benefit;
                target->prev_node[tgt_case] = source;
                target->prev_case[tgt_case] = src_case;
                target->vdd_staples[tgt_case] = new_vdd;
                target->vss_staples[tgt_case] = new_vss;
                target->staples_this_step[tgt_case] = new_staples; // Store for consistency
            }
        }
    }
}

/**
 * @brief Calculate staples for a specific enhanced case
 */
std::vector<Staple> DPSolver::calculateStaplesForCase(EnhancedStapleCase case_type, 
                                                     EnhancedDPNode* node, int row_start) {
    std::vector<Staple> staples;
    int x = node->site * chip_info.site_width;
    
    switch (case_type) {
        case CASE_1_NO_STAPLE:
            // No staples
            break;
            
        case CASE_2_R1_R2_ONLY: {
            int y = chip_info.getRowY(row_start + 1);
            bool is_vdd = determineStapleTypeAdvanced(0, 0, 0, node->site); // Will be refined
            staples.push_back(Staple(x, y, is_vdd));
            break;
        }
        
        case CASE_3_R2_R3_ONLY: {
            int y = chip_info.getRowY(row_start + 2);
            bool is_vdd = determineStapleTypeAdvanced(1, 0, 0, node->site); // Will be refined
            staples.push_back(Staple(x, y, is_vdd));
            break;
        }
        
        case CASE_4_BOTH_ALIGNED: {
            // Both staples at same x position
            int y1 = chip_info.getRowY(row_start + 1);
            int y2 = chip_info.getRowY(row_start + 2);
            
            // Use complementary types for balance
            bool first_is_vdd = determineStapleTypeAdvanced(0, 0, 0, node->site);
            staples.push_back(Staple(x, y1, first_is_vdd));
            staples.push_back(Staple(x, y2, !first_is_vdd));
            break;
        }
        
        case CASE_5_BOTH_STAGGERED: {
            // Staggered staples (slightly offset to avoid manufacturing issues)
            int y1 = chip_info.getRowY(row_start + 1);
            int y2 = chip_info.getRowY(row_start + 2);
            int x2 = x + chip_info.site_width / 2; // Half-site offset
            
            bool first_is_vdd = determineStapleTypeAdvanced(0, 0, 0, node->site);
            staples.push_back(Staple(x, y1, first_is_vdd));
            staples.push_back(Staple(x2, y2, !first_is_vdd));
            break;
        }
    }
    
    return staples;
}

/**
 * @brief Enhanced staggering violation check
 */
bool DPSolver::hasStaggeringViolationEnhanced(const std::vector<Staple>& new_staples,
                                             const std::vector<Staple>& prev_staples,
                                             int current_site, int row_start) {
    
    // Check each new staple against existing staples
    for (const Staple& new_staple : new_staples) {
        int new_site = new_staple.x / chip_info.site_width;
        int new_row = (new_staple.y - chip_info.bottom_y) / chip_info.row_height;
        
        // Check against previous staples in a window
        for (const Staple& prev : prev_staples) {
            int prev_site = prev.x / chip_info.site_width;
            int prev_row = (prev.y - chip_info.bottom_y) / chip_info.row_height;
            
            // Check for staggering pattern: adjacent sites with alternating rows
            if (std::abs(prev_site - new_site) <= 2 && prev_row != new_row) {
                // Look for the forbidden staggering pattern from paper Figure 1(c)
                if (std::abs(prev_site - new_site) == 1) {
                    // Adjacent sites with different rows = staggering
                    return true;
                }
            }
        }
        
        // Check against other new staples
        for (const Staple& other : new_staples) {
            if (&other == &new_staple) continue;
            
            int other_site = other.x / chip_info.site_width;
            int other_row = (other.y - chip_info.bottom_y) / chip_info.row_height;
            
            if (std::abs(other_site - new_site) == 1 && other_row != new_row) {
                return true; // Staggering within new staples
            }
        }
    }
    
    return false;
}

/**
 * @brief Detailed overlap checking with pins and other staples
 */
bool DPSolver::hasOverlapWithPinsOrStaples(const Staple& new_staple,
                                          const std::vector<std::vector<Cell*>>& cells_in_rows,
                                          const std::vector<Staple>& existing_staples,
                                          int row_start) {
    
    // Check overlap with existing staples
    if (overlapWithStaples(new_staple, existing_staples)) {
        return true;
    }
    
    // Check overlap with cell pins
    if (overlapWithPins(new_staple, cells_in_rows, row_start)) {
        return true;
    }
    
    return false;
}

/**
 * @brief Check if staple overlaps with cell pins
 */
bool DPSolver::overlapWithPins(const Staple& staple, 
                              const std::vector<std::vector<Cell*>>& cells_in_rows, 
                              int row_start) {
    
    int staple_site = staple.x / chip_info.site_width;
    int staple_row = (staple.y - chip_info.bottom_y) / chip_info.row_height;
    
    // Check the two rows that this staple connects
    int lower_row_index = staple_row - row_start - 1;
    int upper_row_index = staple_row - row_start;
    
    for (int row_idx : {lower_row_index, upper_row_index}) {
        if (row_idx < 0 || row_idx >= static_cast<int>(cells_in_rows.size())) continue;
        
        // Find cells that might have pins at the staple location
        for (Cell* cell : cells_in_rows[row_idx]) {
            int cell_left_site = cell->current_x / chip_info.site_width;
            const CellType& cell_type = cell_types[cell->type_index];
            int cell_width_sites = cell_type.width / chip_info.site_width;
            int cell_right_site = cell_left_site + cell_width_sites - 1;
            
            // Check if staple site is within cell boundaries
            if (staple_site >= cell_left_site && staple_site <= cell_right_site) {
                int site_offset_in_cell = staple_site - cell_left_site;
                
                // Check if this site has a pin (considering flipping)
                if (cell_type.hasPinAt(site_offset_in_cell, cell->is_flipped)) {
                    std::cout << "Staple inserted at (" << staple.x << "," << staple.y 
                              << ") overlaps with pin(s) or other staple(s)!" << std::endl;
                    return true;
                }
            }
        }
    }
    
    return false;
}

/**
 * @brief Check if staple overlaps with other staples
 */
bool DPSolver::overlapWithStaples(const Staple& staple, const std::vector<Staple>& existing_staples) {
    
    for (const Staple& existing : existing_staples) {
        // Check for exact position overlap
        if (staple.x == existing.x && staple.y == existing.y) {
            std::cout << "Staple inserted at (" << staple.x << "," << staple.y 
                      << ") overlaps with pin(s) or other staple(s)!" << std::endl;
            return true;
        }
        
        // Check for too close proximity (within half a site)
        if (std::abs(staple.x - existing.x) < chip_info.site_width / 2 && 
            std::abs(staple.y - existing.y) < chip_info.row_height / 2) {
            return true;
        }
    }
    
    return false;
}

/**
 * @brief Intelligent node pruning strategy
 */
void DPSolver::pruneNodesIntelligently(int current_site) {
    
    // Collect nodes at recent sites for pruning
    std::vector<std::pair<EnhancedDPNode*, int>> candidates; // (node, best_case_benefit)
    
    for (EnhancedDPNode* node : all_nodes) {
        if (node->site >= current_site - 3) { // Keep recent nodes
            int best_benefit = -1000000;
            int best_case = -1;
            
            for (int c = 0; c < EnhancedDPNode::NUM_CASES; c++) {
                if (node->benefit[c] > best_benefit) {
                    best_benefit = node->benefit[c];
                    best_case = c;
                }
            }
            
            if (best_case >= 0) {
                candidates.push_back({node, best_benefit});
            }
        }
    }
    
    // Sort by benefit (descending)
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) {
                  return a.second > b.second;
              });
    
    // Keep only top BEAM_WIDTH nodes
    if (candidates.size() > BEAM_WIDTH) {
        for (size_t i = BEAM_WIDTH; i < candidates.size(); i++) {
            EnhancedDPNode* node_to_remove = candidates[i].first;
            
            // Remove from lookup table
            EnhancedCompactState compact;
            compact.site = node_to_remove->site;
            for (int j = 0; j < 3; j++) {
                compact.a[j] = node_to_remove->s[j] - initial_s[j];
                compact.b[j] = node_to_remove->displacement[j];
                compact.flipped[j] = node_to_remove->is_flipped[j];
            }
            node_lookup.erase(compact);
            
            // Remove from all_nodes
            auto it = std::find(all_nodes.begin(), all_nodes.end(), node_to_remove);
            if (it != all_nodes.end()) {
                all_nodes.erase(it);
            }
            
            delete node_to_remove;
        }
        
        std::cout << "Pruned " << (candidates.size() - BEAM_WIDTH) 
                  << " nodes, keeping " << BEAM_WIDTH << " best" << std::endl;
    }
}

/**
 * @brief Create or find enhanced node with proper compact encoding
 */
EnhancedDPNode* DPSolver::createOrFindEnhancedNode(int site, int s1, int l1, int s2, int l2, int s3, int l3,
                                                  int disp1, int disp2, int disp3,
                                                  bool flip1, bool flip2, bool flip3) {
    
    // Create enhanced compact state with proper displacement
    EnhancedCompactState compact;
    compact.site = site;
    compact.a[0] = s1 - initial_s[0];
    compact.a[1] = s2 - initial_s[1];
    compact.a[2] = s3 - initial_s[2];
    compact.b[0] = disp1;  // Actual displacement, not l value
    compact.b[1] = disp2;
    compact.b[2] = disp3;
    compact.flipped[0] = flip1;
    compact.flipped[1] = flip2;
    compact.flipped[2] = flip3;
    
    // Check if node already exists
    auto it = node_lookup.find(compact);
    if (it != node_lookup.end()) {
        return it->second;
    }
    
    // Create new enhanced node
    EnhancedDPNode* new_node = new EnhancedDPNode(site, s1, l1, s2, l2, s3, l3);
    new_node->displacement[0] = disp1;
    new_node->displacement[1] = disp2;
    new_node->displacement[2] = disp3;
    new_node->is_flipped[0] = flip1;
    new_node->is_flipped[1] = flip2;
    new_node->is_flipped[2] = flip3;
    
    node_lookup[compact] = new_node;
    all_nodes.push_back(new_node);
    
    return new_node;
}

/**
 * @brief Calculate enhanced case benefit based on paper
 */
int DPSolver::calculateEnhancedCaseBenefit(EnhancedStapleCase case_type, int current_vdd, int current_vss) {
    switch (case_type) {
        case CASE_1_NO_STAPLE:
            // β factor for balance encouragement
            if (current_vdd != current_vss && params.balance_factor > 0) {
                return static_cast<int>(params.balance_factor * 50);
            }
            return 0;
            
        case CASE_2_R1_R2_ONLY:
        case CASE_3_R2_R3_ONLY:
            return 100; // Single staple benefit
            
        case CASE_4_BOTH_ALIGNED:
            return 190; // Slightly less than 2x single due to potential constraints
            
        case CASE_5_BOTH_STAGGERED:
            return 180; // Less than aligned due to manufacturing complexity
            
        default:
            return 0;
    }
}

/**
 * @brief Determine staple type with advanced balancing
 */
bool DPSolver::determineStapleTypeAdvanced(int between_rows, int current_vdd, int current_vss, int site) {
    // Priority 1: Avoid having zero of either type
    if (current_vdd == 0) return true;   // Generate VDD
    if (current_vss == 0) return false;  // Generate VSS
    
    // Priority 2: Balance existing ratio
    if (current_vdd > 0 && current_vss > 0) {
        double ratio = static_cast<double>(std::max(current_vdd, current_vss)) / 
                      static_cast<double>(std::min(current_vdd, current_vss));
        
        if (ratio > 1.05) {
            return (current_vdd < current_vss); // Generate minority type
        }
    }
    
    // Priority 3: Spatial distribution
    bool spatial_vdd = ((site + between_rows) % 2 == 0);
    return spatial_vdd;
}

/**
 * @brief Enhanced solution extraction with proper consistency
 */
std::vector<Staple> DPSolver::extractBestSolutionEnhanced(
    const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start) {
    
    std::cout << "\n=== ENHANCED Solution Extraction ===" << std::endl;
    
    // Find best solution considering all constraints
    std::vector<std::tuple<EnhancedDPNode*, int, int, int, int, double>> candidates;
    
    for (EnhancedDPNode* node : all_nodes) {
        for (int c = 0; c < EnhancedDPNode::NUM_CASES; c++) {
            if (node->benefit[c] > -900000) {
                int vdd = node->vdd_staples[c];
                int vss = node->vss_staples[c];
                int total = vdd + vss;
                
                if (total > 0) {
                    double ratio = (vdd > 0 && vss > 0) ? 
                                  static_cast<double>(std::max(vdd, vss)) / 
                                  static_cast<double>(std::min(vdd, vss)) : 999.0;
                    
                    candidates.push_back({node, c, node->benefit[c], vdd, vss, ratio});
                }
            }
        }
    }
    
    if (candidates.empty()) {
        return generateEnhancedFallback(cells_in_rows, row_start);
    }
    
    // Sort by: 1) Balance constraint satisfaction, 2) Total staples, 3) Benefit
    std::sort(candidates.begin(), candidates.end(),
              [](const auto& a, const auto& b) {
                  double ratio_a = std::get<5>(a);
                  double ratio_b = std::get<5>(b);
                  bool valid_a = (ratio_a <= 1.1);
                  bool valid_b = (ratio_b <= 1.1);
                  
                  if (valid_a != valid_b) return valid_a > valid_b;
                  
                  int total_a = std::get<3>(a) + std::get<4>(a);
                  int total_b = std::get<3>(b) + std::get<4>(b);
                  if (total_a != total_b) return total_a > total_b;
                  
                  return std::get<2>(a) > std::get<2>(b);
              });
    
    auto [best_node, best_case, benefit, vdd_count, vss_count, ratio] = candidates[0];
    
    std::cout << "Enhanced solution: site=" << best_node->site 
              << ", case=" << best_case << ", benefit=" << benefit
              << ", total=" << (vdd_count + vss_count)
              << " (VDD:" << vdd_count << ", VSS:" << vss_count 
              << "), ratio=" << std::fixed << std::setprecision(3) << ratio << std::endl;
    
    // Extract staples using stored information for consistency
    std::vector<Staple> staples;
    EnhancedDPNode* current = best_node;
    int current_case = best_case;
    
    while (current != nullptr && current->prev_node[current_case] != nullptr) {
        // Use stored staple information for consistency
        if (!current->staples_this_step[current_case].empty()) {
            for (const Staple& staple : current->staples_this_step[current_case]) {
                staples.push_back(staple);
            }
        }
        
        EnhancedDPNode* prev = current->prev_node[current_case];
        int prev_case = current->prev_case[current_case];
        current = prev;
        current_case = prev_case;
    }
    
    std::reverse(staples.begin(), staples.end());
    
    std::cout << "Enhanced extraction: " << staples.size() << " staples" << std::endl;
    
    return staples;
}

/**
 * @brief Generate enhanced fallback solution
 */
std::vector<Staple> DPSolver::generateEnhancedFallback(
    const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start) {
    
    std::cout << "Generating enhanced fallback solution..." << std::endl;
    
    std::vector<Staple> staples;
    int vdd_count = 0, vss_count = 0;
    
    // Use wider spacing to avoid overlaps
    int spacing = 12;
    bool next_should_be_vdd = true;
    
    for (int site = spacing; site < chip_info.total_sites - spacing; site += spacing) {
        int x = site * chip_info.site_width;
        
        // R1-R2 staple
        if (row_start + 1 < chip_info.num_rows) {
            int y = chip_info.getRowY(row_start + 1);
            
            Staple candidate(x, y, next_should_be_vdd);
            
            // Check for overlaps before adding
            if (!hasOverlapWithPinsOrStaples(candidate, cells_in_rows, staples, row_start)) {
                staples.push_back(candidate);
                
                if (next_should_be_vdd) vdd_count++;
                else vss_count++;
            }
        }
        
        // Alternate type for next staple
        next_should_be_vdd = !next_should_be_vdd;
        
        // Stop when we have enough balanced staples
        if (staples.size() >= 100) {
            double current_ratio = (vdd_count > 0 && vss_count > 0) ? 
                                  static_cast<double>(std::max(vdd_count, vss_count)) / 
                                  static_cast<double>(std::min(vdd_count, vss_count)) : 1.0;
            
            if (current_ratio <= 1.1) {
                break;
            }
        }
    }
    
    std::cout << "Enhanced fallback: " << staples.size() << " staples (VDD:" << vdd_count 
              << ", VSS:" << vss_count << ")" << std::endl;
    
    return staples;
}

/**
 * @brief Check if row is VDD or VSS
 */
bool DPSolver::isVDDRow(int row_idx) const {
    return (row_idx % 2 == 0);
}

/**
 * @brief Clean up allocated memory
 */
void DPSolver::cleanup() {
    Logger::log(Logger::DEBUG, "Enhanced DPSolver cleanup");
    
    for (EnhancedDPNode* node : all_nodes) {
        if (node != nullptr) {
            delete node;
        }
    }
    all_nodes.clear();
    node_lookup.clear();
    
    while (!node_queue.empty()) {
        node_queue.pop();
    }
    
    nodes_created = 0;
    best_final_node = nullptr;
    best_final_case = -1;
    
    for (int i = 0; i < 3; i++) {
        initial_s[i] = 0;
        for (int j = 0; j < 1000; j++) {
            initial_cell_x[i][j] = 0;
        }
    }
    
    Logger::log(Logger::DEBUG, "Enhanced DPSolver cleanup completed");
}