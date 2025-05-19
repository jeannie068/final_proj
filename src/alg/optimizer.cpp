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
    
    Logger::init("optimizer_log.txt");
    
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
    Logger::log("Starting optimization process");
    
    // Solve the problem by processing triple-row subproblems
    std::vector<Staple> prev_staples;
    Logger::log("Processing " + std::to_string((chip_info.num_rows + 1) / 2) + " triple-row subproblems");
    
    // Use increaseIndent() for nested logs
    Logger::increaseIndent();
    for (int row = 0; row < chip_info.num_rows - 2; row += 2) {
        // Check total execution time
        auto current_time = std::chrono::high_resolution_clock::now();
        auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(
            current_time - global_start_time).count();
            
        if (elapsed_seconds > 570) { // 590 seconds safety limit
            Logger::log("WARNING: Approaching time limit (570s), stopping optimization early");
            std::cout << "WARNING: Approaching time limit, stopping optimization early" << std::endl;
            break;
        }
        
        Logger::log("Processing triple-row subproblem for rows " + std::to_string(row) + 
                    " to " + std::to_string(row+2));
        int row_end = std::min(row + 3, chip_info.num_rows);
        std::vector<Staple> new_staples = solveTripleRow(row, row_end, prev_staples);
        
        // Add new staples to the list
        inserted_staples.insert(inserted_staples.end(), new_staples.begin(), new_staples.end());
        Logger::log("Subproblem completed: " + std::to_string(new_staples.size()) + " staples inserted");
        
        // Update prev_staples for the next iteration
        prev_staples = new_staples;
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
    
    // Extract cells for each row in the triple-row problem
    std::vector<std::vector<Cell*>> cells_in_rows;
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
    
    // Create DPSolver instance
    DPSolver dp_solver(chip_info, cell_types, params);
    
    // Delegate to the solver
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<Staple> inserted_staples = dp_solver.solveTripleRow(
        cells_in_rows, row_start, prev_staples);
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        end_time - start_time);
    
    std::cout << "Triple-row optimization completed in " << duration.count() 
              << " ms. Inserted " << inserted_staples.size() 
              << " staples." << std::endl;
    
    return inserted_staples;
}

/**
 * @brief Determine if a position is a VDD or VSS row
 */
bool Optimizer::isVDDRow(int row_idx) const {
    // Alternate rows for power and ground
    // Even rows (0, 2, 4, ...) are VDD, odd rows are VSS
    return (row_idx % 2 == 0);
}