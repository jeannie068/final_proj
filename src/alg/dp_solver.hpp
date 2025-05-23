/**
 * @file dp_solver.hpp
 * @brief Dynamic programming solver for manufacturing-aware power staple insertion
 * 
 * This file implements the DAG-based dynamic programming approach for triple-row
 * optimization as described in "Manufacturing-Aware Power Staple Insertion 
 * Optimization by Enhanced Multi-Row Detailed Placement Refinement" (ASP-DAC 2021).
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
 * @brief Dynamic programming solver for triple-row optimization
 */
class DPSolver {
public:
    DPSolver(const ChipInfo& chip_info, 
             const std::vector<CellType>& cell_types,
             const AlgorithmParams& params = AlgorithmParams());
    
    ~DPSolver();
    
    /**
     * @brief Solve triple-row optimization problem using MATRO algorithm
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
    
    // DAG structure
    std::queue<DPNode*> node_queue;
    std::unordered_map<CompactState, DPNode*, CompactStateHasher> node_lookup;
    std::vector<DPNode*> all_nodes;  // For memory cleanup
    
    // Tracking
    DPNode* best_final_node;
    int best_final_case;
    size_t nodes_created;
    size_t max_nodes_in_memory;
    
    // Initial cell positions for compact encoding
    int initial_s[3];

    // Pruning
    static const int MAX_NODES_PER_SITE = 2000;
    static const int PRUNING_FREQUENCY = 10;  // Every 10 sites pruning one time
    static const int BENEFIT_THRESHOLD = 3;   // 收益差距閾值
    
    /**
     * @brief Create initial DAG source node
     */
    DPNode *createSourceNode();

    Cell *getNextCellToPlace(int row, int s_j, const std::vector<Cell *> &row_cells);

    Cell *getLastPlacedCell(int row, int s_j, const std::vector<Cell *> &row_cells);

    /**
     * @brief Process nodes at current site and generate extensions
     */
    void processNodesAtSite(int site,
                           const std::vector<std::vector<Cell*>>& cells_in_rows,
                           const std::vector<Staple>& prev_staples,
                           int row_start);

    void pruneNodes(int current_site);

    /**
     * @brief Generate extensions from a node using binary decision tree
     * 
     * Following Section 3.1: up to 8 extensions (2^3) based on place/defer decisions
     */
    void generateExtensions(DPNode* node,
                           const std::vector<std::vector<Cell*>>& cells_in_rows,
                           const std::vector<Staple>& prev_staples,
                           int row_start);
    
    /**
     * @brief Check if we can defer placing the next cell beyond site i+1
     */
    bool canDeferNextCell(int row, int s_j, int l_j, int site,
                         const std::vector<Cell*>& row_cells);
    
    /**
     * @brief Check if we can place the next cell at site i+1
     */
    bool canPlaceNextCell(int row, int s_j, int l_j, int site,
                         const std::vector<Cell*>& row_cells);
    
    /**
     * @brief Create or find node with given state
     */
    DPNode* createOrFindNode(DPNode* node, const std::vector<std::vector<Cell*>>& cells_in_rows, int site, int s1, int l1, int s2, int l2, 
                            int s3, int l3);
    
    /**
     * @brief Update benefit for target node based on source node
     */
    void updateBenefit(DPNode* source, DPNode* target,
                      const std::vector<std::vector<Cell*>>& cells_in_rows,
                      const std::vector<Staple>& prev_staples,
                      int row_start);
    
    /**
     * @brief Check if a staple can be inserted (no pin obstruction)
     */
    bool canInsertStaple(int site, int between_rows, 
                        int s[3], int l[3],
                        const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    /**
     * @brief Check for anti-parallel line-ends violation
     */
    bool hasAntiParallelViolation(int site, int staple_case, int prev_case,
                                 const std::vector<Staple>& prev_staples,
                                 int row_start);
    
    /**
     * @brief Determine which of the 5 staple cases are valid at this site
     */
    std::vector<int> getValidStapleCases(int site, int s[3], int l[3],
                                       const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    /**
     * @brief Calculate staple benefit with balance factor
     */
    int calculateStapleBenefit(int staple_case, int vdd_count, int vss_count,
                              int row_start);
    
    /**
     * @brief Check if row is VDD or VSS
     */
    bool isVDDRow(int row_idx) const;
    
    /**
     * @brief Backtrack from best node to recover solution
     */
    std::vector<Staple> backtrackSolution(
        const std::vector<std::vector<Cell*>>& cells_in_rows,
        int row_start);
    
    /**
     * @brief Apply cell placement updates from backtracking
     */
    void applyCellPlacements(const std::vector<std::vector<Cell*>>& cells_in_rows);

    /**
     * @brief Clean up allocated memory
     */
    void cleanup();
    
    /**
     * @brief Get initial s values for compact encoding
     */
    void computeInitialS(const std::vector<std::vector<Cell*>>& cells_in_rows);
};

class StuckDetector {
private:
    std::chrono::high_resolution_clock::time_point last_progress_time;
    int last_site_processed;
    static const int STUCK_TIMEOUT_SECONDS = 30;  // 30秒沒進展視為stuck
    
public:
    StuckDetector() {
        reset();
    }
    
    void reset() {
        last_progress_time = std::chrono::high_resolution_clock::now();
        last_site_processed = -1;
    }
    
    void updateProgress(int current_site) {
        if (current_site > last_site_processed) {
            last_progress_time = std::chrono::high_resolution_clock::now();
            last_site_processed = current_site;
        }
    }
    
    bool isStuck() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now - last_progress_time);
        return duration.count() > STUCK_TIMEOUT_SECONDS;
    }
    
    int getStuckDuration() {
        auto now = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(now - last_progress_time);
        return duration.count();
    }
};

#endif // DP_SOLVER_HPP