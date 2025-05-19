/**
 * @file dp_solver.hpp
 * @brief Dynamic programming solver for power staple insertion optimization
 * 
 * This file contains the declarations for the DAG-based dynamic programming
 * solver used for power staple insertion optimization, implementing the
 * approach described in "Manufacturing-Aware Power Staple Insertion Optimization
 * by Enhanced Multi-Row Detailed Placement Refinement".
 */

#ifndef DP_SOLVER_HPP
#define DP_SOLVER_HPP

#include "../data_structure/data_structure.hpp"
#include <unordered_map>
#include <queue>
#include <vector>
#include <algorithm>
#include <memory>
#include <chrono>
#include <iostream>

/**
 * @brief Dynamic programming solver for power staple insertion optimization
 */
class DPSolver {
public:
    /**
     * @brief Constructor
     * @param chip_info Chip information
     * @param cell_types Vector of cell types
     * @param params Algorithm parameters
     */
    DPSolver(const ChipInfo& chip_info, 
             const std::vector<CellType>& cell_types,
             const AlgorithmParams& params = AlgorithmParams());
    
    /**
     * @brief Destructor
     */
    ~DPSolver();
    
    /**
     * @brief Solve a triple-row optimization problem
     * @param cells_in_rows Cells in each of the three rows being optimized
     * @param row_start Start row index in the global layout
     * @param prev_staples Staples inserted in previous optimization rounds
     * @return Vector of newly inserted staples
     */
    std::vector<Staple> solveTripleRow(
        const std::vector<std::vector<Cell*>>& cells_in_rows,
        int row_start,
        const std::vector<Staple>& prev_staples);
    
private:
    // Input data
    ChipInfo chip_info;
    std::vector<CellType> cell_types;
    AlgorithmParams params;
    
    // Memory management
    std::unordered_map<CompactState, DPNode*, CompactStateHasher> node_lookup_table;
    std::queue<DPNode*> node_queue;
    std::vector<DPNode*> all_nodes;  // For cleanup
    
    // Best node tracking
    DPNode* best_node;
    int max_benefit;
    
    // Stats for memory monitoring
    size_t max_nodes_count;
    size_t nodes_created;
    size_t nodes_processed;
    
    /**
     * @brief Initialize the DAG with a source node
     * @return Source node for the DAG
     */
    DPNode* initializeDAG();
    
    /**
     * @brief Process all nodes at a given site
     * @param site Current site
     * @param cells_in_rows Cells in each row
     * @param prev_staples Staples from previous rounds
     * @param row_start Global start row index
     */
    void processSite(int site,
                     const std::vector<std::vector<Cell*>>& cells_in_rows,
                     const std::vector<Staple>& prev_staples,
                     int row_start);
    
    /**
     * @brief Generate all valid extensions for a node
     * @param node Current node
     * @param site Current site
     * @param cells_in_rows Cells in each row
     * @param row_start Global start row index
     * @return Vector of valid extensions
     */
    std::vector<Extension> generateExtensions(
        DPNode* node, 
        int site,
        const std::vector<std::vector<Cell*>>& cells_in_rows,
        int row_start);
    
    /**
     * @brief Check if a staple can be inserted at a given site
     * @param site Site index
     * @param row Row index (0-based within triple-row problem)
     * @param cells_in_rows Cells in each row
     * @param node Current node
     * @return true if staple can be inserted, false otherwise
     */
    bool canInsertStaple(int site, 
                         int row, 
                         const std::vector<std::vector<Cell*>>& cells_in_rows,
                         DPNode* node);
    
    /**
     * @brief Check for staggering violations
     * @param site Site index
     * @param staple_case Current staple insertion case
     * @param prev_case Previous staple insertion case
     * @param prev_staples Staples from previous rounds
     * @param row_start Global start row index
     * @return true if staggering violation exists, false otherwise
     */
    bool hasStaggeringViolation(int site,
                               int staple_case,
                               int prev_case,
                               const std::vector<Staple>& prev_staples,
                               int row_start);
    
    /**
     * @brief Get or create a node with given state
     * @param state Compact state encoding
     * @return Pointer to new or existing node
     */
    DPNode* getOrCreateNode(const CompactState& state);
    
    /**
     * @brief Update node benefit based on extension
     * @param from_node Source node
     * @param to_node Target node
     * @param site Current site
     * @param extension Current extension
     * @param prev_staples Staples from previous rounds
     * @param row_start Global start row index
     * @param cells_in_rows Cells in each row
     */
    void updateBenefit(DPNode* from_node,
                      DPNode* to_node,
                      int site,
                      const Extension& extension,
                      const std::vector<Staple>& prev_staples,
                      int row_start,
                      const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    /**
     * @brief Backtrack to get the optimal solution
     * @param cells_in_rows Cells in each row
     * @param row_start Global start row index
     * @return Vector of inserted staples
     */
    std::vector<Staple> backtrack(
        const std::vector<std::vector<Cell*>>& cells_in_rows,
        int row_start);

    void pruneNodes(int site);

    std::vector<Staple> solveTripleRowWithWindows(const std::vector<std::vector<Cell *>> &cells_in_rows, int row_start, const std::vector<Staple> &prev_staples);

    std::vector<Staple> solveTripleRowWindow(const std::vector<std::vector<Cell *>> &window_cells, int row_start, const std::vector<Staple> &window_prev_staples, int window_start, int window_end);

    /**
     * @brief Apply final cell placement updates
     * @param cells_in_rows Cells in each row
     */
    void applyPlacementUpdates(
        const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    /**
     * @brief Check if cell placement is valid
     * @param cell Cell to check
     * @param x X-coordinate
     * @param y Y-coordinate
     * @param is_flipped Whether cell is flipped
     * @return true if placement is valid, false otherwise
     */
    bool isValidPlacement(const Cell* cell, int x, int y, bool is_flipped);
    
    /**
     * @brief Get compact state encoding
     * @param site Current site
     * @param s1,l1,f1 Row 1 cell state (index, distance, flip)
     * @param s2,l2,f2 Row 2 cell state
     * @param s3,l3,f3 Row 3 cell state
     * @param cells_in_rows Cells in each row
     * @return Compact state encoding
     */
    CompactState getCompactState(
        int site,
        int s1, int l1, bool f1,
        int s2, int l2, bool f2,
        int s3, int l3, bool f3,
        const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    /**
     * @brief Check if row is VDD or VSS
     * @param row_idx Row index
     * @return true if VDD row, false if VSS row
     */
    bool isVDDRow(int row_idx) const;
    
    /**
     * @brief Calculate staple benefit with balance constraints
     * @param staple_case Staple insertion case
     * @param current_vdd Current VDD staple count
     * @param current_vss Current VSS staple count
     * @param row_start Global start row index
     * @return Adjusted benefit value
     */
    double calculateStapleBenefit(
        int staple_case,
        int current_vdd,
        int current_vss,
        int row_start);
    
    /**
     * @brief Clean up allocated memory
     */
    void cleanup();

    /**
     * @brief Monitor memory usage and report statistics
     */
    void monitorMemoryUsage();

    bool shouldApplyBalanceFactor(int current_vdd, int current_vss, int row_start, 
                             bool can_insert_r1_r2, bool can_insert_r2_r3);
    void applyStapleBalanceAdjustment(int& benefit, int vdd_count, int vss_count, 
                                 int staple_case, int row_start);

    bool isPromisingSolution(const CompactState &state);

    std::vector<std::vector<Cell *>> extractCellsForWindow(const std::vector<std::vector<Cell *>> &cells_in_rows, int window_start, int window_end);

    std::vector<Staple> extractStaplesForWindow(const std::vector<Staple> &staples, int window_start, int window_end);

    void releaseNodesBeforeSite(int site);
};

#endif // DP_SOLVER_HPP