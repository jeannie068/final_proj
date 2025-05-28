#include "optimizer.hpp"
#include "../Logger.hpp"
#include "initial_staple_generator.hpp"

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
    auto global_start_time = std::chrono::high_resolution_clock::now();
    Logger::log(Logger::INFO, "Starting optimization process");
    
    // Option 1: Use simple initial solution (no cell refinement)
    bool USE_SIMPLE_INITIAL = true;  // Toggle this for testing
    
    if (USE_SIMPLE_INITIAL) {
        std::cout << "Using simple initial staple generation (no cell refinement)" << std::endl;
        
        // Generate initial staples without moving cells
        inserted_staples = generateSimpleInitialSolution(chip_info, cell_types, cells);
        
        // Create solution with original cell positions
        Solution solution;
        solution.refined_cells = cells;  // Keep original positions
        solution.inserted_staples = inserted_staples;
        // solution.inserted_staples = generateSimpleInitialStaplesFixed(chip_info, cell_types, cells);
        solution.updateStats();
        
        auto end_time = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - global_start_time);
        
        std::cout << "Initial solution generated in " << duration.count() << " ms" << std::endl;
        std::cout << "Total staples: " << solution.total_staples << std::endl;
        std::cout << "Balance constraint satisfied: " 
                  << (solution.isBalanceConstraintSatisfied() ? "Yes" : "No") << std::endl;
        
        return solution;
    }
    
    // Option 2: Original complex DP algorithm
    // 清空staples列表 - 確保乾淨開始
    inserted_staples.clear();
    // Solve the problem by processing triple-row subproblems
    std::vector<Staple> prev_staples;
    Logger::log("Processing " + std::to_string((chip_info.num_rows + 1) / 2) + " triple-row subproblems");
    
    int total_new_staples = 0;  // 新增：追蹤總的新增staples

    // Use increaseIndent() for nested logs
    Logger::increaseIndent();
    // for (int row = 0; row < 12; row += 2) {
    for (int row = 0; row < chip_info.num_rows - 2; row += 2) {
        // Check total execution time
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - global_start_time).count();
            
        if (elapsed_seconds > 500) {
            Logger::log("WARNING: Approaching time limit (500s), stopping optimization early");
            std::cout << "WARNING: Approaching time limit, stopping optimization early" << std::endl;
            break;
        }
        
        Logger::log("Processing triple-row subproblem for rows " + std::to_string(row) + 
                    " to " + std::to_string(row+2));
        int row_end = std::min(row + 3, chip_info.num_rows);
        
        // 解決subproblem並獲得這一輪的新staples
        std::vector<Staple> new_staples_this_round = solveTripleRow(row, row_end, prev_staples);
        
        // 修正：正確累積staples
        inserted_staples.insert(inserted_staples.end(), 
                               new_staples_this_round.begin(), 
                               new_staples_this_round.end());
        
        total_new_staples += new_staples_this_round.size();
        
        Logger::log("Cumulative staples so far: " + std::to_string(inserted_staples.size()));
        std::cout << "Cumulative total staples: " << inserted_staples.size() << std::endl;
        
        // 更新prev_staples for next iteration - 這裡可能需要包含所有之前的staples
        prev_staples.insert(prev_staples.end(), 
                           new_staples_this_round.begin(), 
                           new_staples_this_round.end());
    }
    Logger::decreaseIndent();
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - global_start_time);
    
    std::cout << "Optimization completed in " << duration.count() << " ms" << std::endl;
    Logger::log("Optimization completed with " + std::to_string(inserted_staples.size()) + " staples");
    
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
    Logger::log("Starting triple-row optimization for rows " + std::to_string(row_start) + 
                " to " + std::to_string(row_end-1));
    
    // 提取cells for each row in the triple-row problem
    std::vector<std::vector<Cell*>> cells_in_rows;
    for (int r = row_start; r < row_end; r++) {
        if (r < static_cast<int>(cells_by_row.size())) {
            cells_in_rows.push_back(cells_by_row[r]);
        } else {
            cells_in_rows.push_back(std::vector<Cell*>());
        }
    }
    
    // Ensure we have exactly 3 rows
    while (cells_in_rows.size() < 3) {
        cells_in_rows.push_back(std::vector<Cell*>());
    }
    
    // Create DPSolver instance - 關鍵：每次創建新的instance
    DPSolver dp_solver(chip_info, cell_types, params);
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // 解決這個triple-row subproblem
    std::vector<Staple> new_staples_this_round = dp_solver.solveTripleRow(
        cells_in_rows, row_start, prev_staples);

    // std::vector<Staple> new_staples_this_round = dp_solver.solveTripleRowMinimal(
    //    cells_in_rows, row_start, prev_staples);
    
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);
    
    // 修正：顯示THIS round的實際新增staples
    std::cout << "Triple-row optimization completed in " << duration.count() 
              << " ms. Inserted " << new_staples_this_round.size() 
              << " NEW staples (this round only)." << std::endl;
    
    Logger::log("Triple-row subproblem completed: " + 
                std::to_string(new_staples_this_round.size()) + " NEW staples inserted");
    
    // 重要：只返回這一輪新增的staples
    return new_staples_this_round;
}

/**
 * @brief Determine if a position is a VDD or VSS row
 */
bool Optimizer::isVDDRow(int row_idx) const {
    // Alternate rows for power and ground
    // Even rows (0, 2, 4, ...) are VDD, odd rows are VSS
    return (row_idx % 2 == 0);
}

/**
 * @brief Simple standalone function with enhanced overlap checking
 * Direct replacement for generateSimpleInitialStaples
 */
/**
 * @brief Rewritten function for aggressive staple generation
 */
std::vector<Staple> generateSimpleInitialStaplesFixed(const ChipInfo& chip_info,
                                                     const std::vector<CellType>& cell_types,
                                                     const std::vector<Cell>& cells) {
    std::vector<Staple> staples;
    
    std::cout << "\n=== Aggressive Initial Staple Generation ===" << std::endl;
    std::cout << "Chip: " << chip_info.num_rows << " rows, " 
              << chip_info.total_sites << " sites/row" << std::endl;
    
    // Step 1: Build pin occupation map
    std::vector<std::vector<bool>> has_pin(chip_info.num_rows, 
                                          std::vector<bool>(chip_info.total_sites, false));
    
    int total_pins_marked = 0;
    for (const Cell& cell : cells) {
        int row = cell.initial_y / chip_info.row_height;
        if (row < 0 || row >= chip_info.num_rows) continue;
        
        const CellType& ctype = cell_types[cell.type_index];
        int start_site = cell.initial_x / chip_info.site_width;
        
        // Only mark actual pin sites, not entire cell
        for (int pin_offset : ctype.pin_sites) {
            int pin_site = start_site + pin_offset;
            if (pin_site >= 0 && pin_site < chip_info.total_sites) {
                has_pin[row][pin_site] = true;
                total_pins_marked++;
            }
        }
    }
    std::cout << "Total pins marked: " << total_pins_marked << std::endl;
    
    // Calculate available sites per row
    std::vector<int> available_sites_per_row(chip_info.num_rows - 1, 0);
    for (int row_boundary = 1; row_boundary < chip_info.num_rows; row_boundary++) {
        for (int site = 0; site < chip_info.total_sites; site++) {
            if (!has_pin[row_boundary-1][site] && !has_pin[row_boundary][site]) {
                available_sites_per_row[row_boundary-1]++;
            }
        }
    }
    
    int total_available = 0;
    for (int avail : available_sites_per_row) {
        total_available += avail;
    }
    std::cout << "Total available positions for staples: " << total_available << std::endl;
    
    // Step 2: Collect ALL possible staple positions
    struct StapleCandidate {
        int site;
        int row_boundary;
        int x;
        int y;
    };
    std::vector<StapleCandidate> all_candidates;
    
    for (int row_boundary = 1; row_boundary < chip_info.num_rows; row_boundary++) {
        for (int site = 0; site < chip_info.total_sites; site++) {
            if (!has_pin[row_boundary-1][site] && !has_pin[row_boundary][site]) {
                StapleCandidate cand;
                cand.site = site;
                cand.row_boundary = row_boundary;
                cand.x = site * chip_info.site_width;
                cand.y = chip_info.getRowY(row_boundary);
                all_candidates.push_back(cand);
            }
        }
    }
    
    std::cout << "Found " << all_candidates.size() << " candidate positions" << std::endl;
    
    // Step 3: Greedy placement with staggering check
    std::set<std::pair<int, int>> staple_locs; // (site, row_boundary)
    int vdd_count = 0, vss_count = 0;
    int rejected_stagger = 0;
    int rejected_balance = 0;
    
    // Sort candidates to spread them out evenly
    std::sort(all_candidates.begin(), all_candidates.end(), 
              [](const StapleCandidate& a, const StapleCandidate& b) {
                  // Sort by (site + row) to get diagonal distribution
                  return (a.site + a.row_boundary) < (b.site + b.row_boundary);
              });
    
    // Try to place as many staples as possible
    for (size_t i = 0; i < all_candidates.size(); i++) {
        const StapleCandidate& cand = all_candidates[i];
        
        // Determine type for balance
        bool is_vdd;
        if (vdd_count == 0) {
            is_vdd = true;
        } else if (vss_count == 0) {
            is_vdd = false;
        } else {
            double current_ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                                  static_cast<double>(std::min(vdd_count, vss_count));
            if (current_ratio > 1.05) {
                // Force minority type
                is_vdd = (vdd_count < vss_count);
            } else {
                // Alternate based on position for spatial distribution
                is_vdd = ((cand.site + cand.row_boundary) % 2 == 0);
            }
        }
        
        // Check future balance
        int future_vdd = vdd_count + (is_vdd ? 1 : 0);
        int future_vss = vss_count + (is_vdd ? 0 : 1);
        if (future_vdd > 0 && future_vss > 0) {
            double future_ratio = static_cast<double>(std::max(future_vdd, future_vss)) / 
                                 static_cast<double>(std::min(future_vdd, future_vss));
            if (future_ratio > 1.095) { // Very tight margin
                rejected_balance++;
                // Try opposite type
                is_vdd = !is_vdd;
                future_vdd = vdd_count + (is_vdd ? 1 : 0);
                future_vss = vss_count + (is_vdd ? 0 : 1);
                future_ratio = static_cast<double>(std::max(future_vdd, future_vss)) / 
                              static_cast<double>(std::min(future_vdd, future_vss));
                if (future_ratio > 1.095) {
                    continue; // Skip this position
                }
            }
        }
        
        // Check for staggering violations
        bool has_stagger = false;
        std::vector<std::pair<int, int>> diagonal_neighbors;
        
        // Find all diagonal neighbors
        for (const auto& existing : staple_locs) {
            if (std::abs(cand.site - existing.first) == 1 && 
                std::abs(cand.row_boundary - existing.second) == 1) {
                diagonal_neighbors.push_back(existing);
            }
        }
        
        // Check each diagonal neighbor for staggering
        for (const auto& neighbor : diagonal_neighbors) {
            int lower_site = (cand.row_boundary < neighbor.second) ? cand.site : neighbor.first;
            int lower_row = std::min(cand.row_boundary, neighbor.second);
            int upper_site = (cand.row_boundary < neighbor.second) ? neighbor.first : cand.site;
            int upper_row = std::max(cand.row_boundary, neighbor.second);
            
            // Check for blockers
            bool blocker_above = staple_locs.count({lower_site, lower_row + 1}) > 0;
            bool blocker_below = staple_locs.count({upper_site, upper_row - 1}) > 0;
            
            if (!blocker_above && !blocker_below) {
                has_stagger = true;
                break;
            }
        }
        
        if (has_stagger) {
            rejected_stagger++;
            // Try to find a nearby position that doesn't create staggering
            bool found_alternative = false;
            
            // Try adjacent rows first
            for (int row_offset : {-1, 1}) {
                int alt_row = cand.row_boundary + row_offset;
                if (alt_row >= 1 && alt_row < chip_info.num_rows &&
                    !has_pin[alt_row-1][cand.site] && !has_pin[alt_row][cand.site] &&
                    staple_locs.count({cand.site, alt_row}) == 0) {
                    
                    // Check if this alternative creates staggering
                    bool alt_has_stagger = false;
                    for (const auto& existing : staple_locs) {
                        if (std::abs(cand.site - existing.first) == 1 && 
                            std::abs(alt_row - existing.second) == 1) {
                            // Check blockers for alternative position
                            int l_site = (alt_row < existing.second) ? cand.site : existing.first;
                            int l_row = std::min(alt_row, existing.second);
                            int u_site = (alt_row < existing.second) ? existing.first : cand.site;
                            int u_row = std::max(alt_row, existing.second);
                            
                            bool b_above = staple_locs.count({l_site, l_row + 1}) > 0;
                            bool b_below = staple_locs.count({u_site, u_row - 1}) > 0;
                            
                            if (!b_above && !b_below) {
                                alt_has_stagger = true;
                                break;
                            }
                        }
                    }
                    
                    if (!alt_has_stagger) {
                        // Use alternative position
                        Staple alt_staple(cand.site * chip_info.site_width, 
                                        chip_info.getRowY(alt_row), is_vdd);
                        staples.push_back(alt_staple);
                        staple_locs.insert({cand.site, alt_row});
                        if (is_vdd) vdd_count++;
                        else vss_count++;
                        found_alternative = true;
                        break;
                    }
                }
            }
            
            if (!found_alternative) {
                continue; // Skip this position
            }
        } else {
            // No staggering, place the staple
            Staple new_staple(cand.x, cand.y, is_vdd);
            staples.push_back(new_staple);
            staple_locs.insert({cand.site, cand.row_boundary});
            if (is_vdd) vdd_count++;
            else vss_count++;
        }
        
        // Progress report every 1000 staples
        if (staples.size() % 1000 == 0 && staples.size() > 0) {
            std::cout << "Progress: " << staples.size() << " staples placed (VDD:" 
                      << vdd_count << ", VSS:" << vss_count << ")" << std::endl;
        }
    }
    
    // Step 4: Post-processing to fix any remaining issues
    std::cout << "\nPost-processing to maximize staple count..." << std::endl;
    
    // Try to fill gaps where possible
    int added_in_gaps = 0;
    for (int site = 0; site < chip_info.total_sites; site++) {
        for (int row_boundary = 1; row_boundary < chip_info.num_rows; row_boundary++) {
            // Skip if already has staple or has pins
            if (staple_locs.count({site, row_boundary}) > 0 ||
                has_pin[row_boundary-1][site] || has_pin[row_boundary][site]) {
                continue;
            }
            
            // Check if we can add without creating staggering
            bool can_add = true;
            for (const auto& existing : staple_locs) {
                if (std::abs(site - existing.first) == 1 && 
                    std::abs(row_boundary - existing.second) == 1) {
                    // Check for blockers
                    int l_site = (row_boundary < existing.second) ? site : existing.first;
                    int l_row = std::min(row_boundary, existing.second);
                    int u_site = (row_boundary < existing.second) ? existing.first : site;
                    int u_row = std::max(row_boundary, existing.second);
                    
                    bool b_above = staple_locs.count({l_site, l_row + 1}) > 0;
                    bool b_below = staple_locs.count({u_site, u_row - 1}) > 0;
                    
                    if (!b_above && !b_below) {
                        can_add = false;
                        break;
                    }
                }
            }
            
            if (can_add) {
                // Determine type to maintain balance
                bool is_vdd = (vdd_count <= vss_count);
                
                // Check balance constraint
                int future_vdd = vdd_count + (is_vdd ? 1 : 0);
                int future_vss = vss_count + (is_vdd ? 0 : 1);
                if (future_vdd > 0 && future_vss > 0) {
                    double ratio = static_cast<double>(std::max(future_vdd, future_vss)) / 
                                  static_cast<double>(std::min(future_vdd, future_vss));
                    if (ratio <= 1.1) {
                        Staple gap_staple(site * chip_info.site_width, 
                                        chip_info.getRowY(row_boundary), is_vdd);
                        staples.push_back(gap_staple);
                        staple_locs.insert({site, row_boundary});
                        if (is_vdd) vdd_count++;
                        else vss_count++;
                        added_in_gaps++;
                    }
                }
            }
        }
    }
    
    // Final statistics
    std::cout << "\nGeneration completed:" << std::endl;
    std::cout << "  Total candidates: " << all_candidates.size() << std::endl;
    std::cout << "  Staples placed: " << staples.size() << std::endl;
    std::cout << "  Rejected for staggering: " << rejected_stagger << std::endl;
    std::cout << "  Rejected for balance: " << rejected_balance << std::endl;
    std::cout << "  Added in gaps: " << added_in_gaps << std::endl;
    std::cout << "  VDD: " << vdd_count << ", VSS: " << vss_count << std::endl;
    if (vdd_count > 0 && vss_count > 0) {
        double ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                      static_cast<double>(std::min(vdd_count, vss_count));
        std::cout << "  Balance ratio: " << std::fixed << std::setprecision(3) << ratio << std::endl;
    }
    
    return staples;
}