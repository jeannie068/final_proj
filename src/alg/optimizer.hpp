/**
 * @file optimizer.hpp
 * @brief Core algorithm for power staple insertion optimization
 * 
 * This file contains the implementation of the DAG-based dynamic programming
 * approach for power staple insertion optimization as described in the paper
 * "Manufacturing-Aware Power Staple Insertion Optimization by Enhanced Multi-Row
 * Detailed Placement Refinement".
 */

#ifndef OPTIMIZER_HPP
#define OPTIMIZER_HPP

#include "../data_structure/data_structure.hpp"
#include "dp_solver.hpp"
#include <unordered_map>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <memory>

std::vector<Staple> generateSimpleInitialStaplesFixed(const ChipInfo& chip_info,
                                                     const std::vector<CellType>& cell_types,
                                                     const std::vector<Cell>& cells);

/**
 * @brief Main optimizer class for power staple insertion
 */
class Optimizer {
public:
    /**
     * @brief Constructor
     * @param chip_info Chip information
     * @param cell_types Vector of cell types
     * @param cells Vector of cells
     * @param params Algorithm parameters
     */
    Optimizer(const ChipInfo& chip_info,
              const std::vector<CellType>& cell_types,
              const std::vector<Cell>& cells,
              const AlgorithmParams& params = AlgorithmParams());
    
    /**
     * @brief Run the optimization algorithm
     * @return Solution containing refined placement and inserted staples
     */
    Solution run();
    
private:
    // Input data
    ChipInfo chip_info;
    std::vector<CellType> cell_types;
    std::vector<Cell> cells;
    AlgorithmParams params;
    
    // Problem state
    std::vector< std::vector<Cell*> > cells_by_row;  // Cells organized by row
    std::vector<Staple> inserted_staples;          // All inserted staples
    
    /**
     * @brief Solve a triple-row optimization problem
     * @param row_start Start row index
     * @param row_end End row index (exclusive)
     * @param prev_staples Staples inserted in previous rounds
     * @return Vector of staples inserted in this round
     */
    std::vector<Staple> solveTripleRow(int row_start, int row_end, 
                                     const std::vector<Staple>& prev_staples);
    
    /**
     * @brief Determine if a position is a VDD or VSS row
     * @param row_idx Row index
     * @return true if it's a VDD row, false if it's a VSS row
     */
    bool isVDDRow(int row_idx) const;
};

#endif // OPTIMIZER_HPP