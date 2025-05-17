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
#include <unordered_map>
#include <queue>
#include <vector>
#include <algorithm>
#include <iostream>
#include <chrono>
#include <memory>

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
     * @brief Initialize the DAG for a triple-row problem
     * @param row_start Start row index
     * @param row_end End row index (exclusive)
     * @param cells_in_rows Cells in each of the three rows
     * @return Source node of the DAG
     */
    DPNode* initializeDAG(int row_start, int row_end, 
                        const std::vector< std::vector<Cell*> >& cells_in_rows);
    
    /**
     * @brief Process all nodes at a given site
     * @param site Current site index
     * @param queue Queue of nodes to process
     * @param node_lookup_table Lookup table for nodes
     * @param cells_in_rows Cells in each of the three rows
     * @param prev_staples Staples inserted in previous rounds
     * @param row_start Start row index
     */
    void processSite(int site, std::queue<DPNode*>& queue,
                    std::unordered_map<CompactState, DPNode*, CompactStateHasher>& node_lookup_table,
                    const std::vector< std::vector<Cell*> >& cells_in_rows,
                    const std::vector<Staple>& prev_staples,
                    int row_start);
    
    /**
     * @brief Generate all valid extensions for a node
     * @param node Current node
     * @param site Current site
     * @param cells_in_rows Cells in each of the three rows
     * @return Vector of valid extensions
     */
    std::vector<Extension> generateExtensions(DPNode* node, int site,
                                           const std::vector< std::vector<Cell*> >& cells_in_rows);
    
    /**
     * @brief Check if a staple can be inserted at a given position
     * @param site Site index
     * @param row Row index (0-based, relative to the triple-row problem)
     * @param cells_in_rows Cells in each row
     * @param node Current node
     * @return true if a staple can be inserted, false otherwise
     */
    bool canInsertStaple(int site, int row, 
                        const std::vector< std::vector<Cell*> >& cells_in_rows,
                        DPNode* node);
    
    /**
     * @brief Check if there's a staggering violation at the given site
     * @param site Site index
     * @param staple_case Staple case being considered
     * @param prev_case Previous staple case
     * @param prev_staples Staples inserted in previous rounds
     * @param row_start Start row index
     * @return true if there's a staggering violation, false otherwise
     */
    bool hasStaggeringViolation(int site, int staple_case, int prev_case,
                              const std::vector<Staple>& prev_staples,
                              int row_start);
    
    /**
     * @brief Create or get a node with the given configuration
     * @param state State of the node
     * @param queue Queue of nodes to process
     * @param node_lookup_table Lookup table for nodes
     * @return Pointer to the new or existing node
     */
    DPNode* getOrCreateNode(const CompactState& state,
                          std::queue<DPNode*>& queue,
                          std::unordered_map<CompactState, DPNode*, CompactStateHasher>& node_lookup_table);
    
    /**
     * @brief Update benefit of a node based on an extension
     * @param from_node Source node
     * @param to_node Target node
     * @param site Current site
     * @param extension Extension being considered
     * @param prev_staples Staples inserted in previous rounds
     * @param row_start Start row index
     */
    void updateBenefit(DPNode* from_node, DPNode* to_node, int site,
                    const Extension& extension,
                    const std::vector<Staple>& prev_staples,
                    int row_start,
                    const std::vector< std::vector<Cell*> >& cells_in_rows);
    
    /**
     * @brief Backtrack through the DAG to get the optimal solution
     * @param best_node Node with the best benefit
     * @param cells_in_rows Cells in each row
     * @param row_start Start row index
     * @return Vector of staples to be inserted
     */
    std::vector<Staple> backtrack(DPNode* best_node,
                                const std::vector< std::vector<Cell*> >& cells_in_rows,
                                int row_start);
    
    /**
     * @brief Check if a cell placement is valid
     * @param cell Cell to check
     * @param x X-coordinate to place the cell
     * @param y Y-coordinate to place the cell
     * @param is_flipped Whether the cell is flipped
     * @return true if the placement is valid, false otherwise
     */
    bool isValidPlacement(const Cell* cell, int x, int y, bool is_flipped);
    
    /**
     * @brief Get the compact state encoding for a configuration
     * @param site Current site
     * @param s1, l1, f1 Row 1 configuration (cell index, distance, flip flag)
     * @param s2, l2, f2 Row 2 configuration
     * @param s3, l3, f3 Row 3 configuration
     * @param cells_in_rows Cells in each row
     * @return Compact state encoding
     */
    CompactState getCompactState(int site, 
                              int s1, int l1, bool f1,
                              int s2, int l2, bool f2,
                              int s3, int l3, bool f3,
                              const std::vector< std::vector<Cell*> >& cells_in_rows);
    
    /**
     * @brief Determine if a position is a VDD or VSS row
     * @param row_idx Row index
     * @return true if it's a VDD row, false if it's a VSS row
     */
    bool isVDDRow(int row_idx) const;
};

#endif // OPTIMIZER_HPP