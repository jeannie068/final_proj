#include "optimizer.hpp"
#include "../Logger.hpp"

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
    
    // 清空staples列表 - 確保乾淨開始
    inserted_staples.clear();
    // Solve the problem by processing triple-row subproblems
    std::vector<Staple> prev_staples;
    Logger::log("Processing " + std::to_string((chip_info.num_rows + 1) / 2) + " triple-row subproblems");
    
    int total_new_staples = 0;  // 新增：追蹤總的新增staples

    // Use increaseIndent() for nested logs
    Logger::increaseIndent();
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