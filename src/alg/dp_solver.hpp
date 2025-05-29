/**
 * @file dp_solver.hpp
 * @brief Corrected dynamic programming solver for manufacturing-aware power staple insertion
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
#include <iomanip>

/**
 * @brief Five staple cases from paper Figure 7
 */
enum StapleCase {
    CASE_1_NO_STAPLE = 0,
    CASE_2_R1_R2_ONLY = 1,
    CASE_3_R2_R3_ONLY = 2,
    CASE_4_BOTH_ALIGNED = 3,
    CASE_5_BOTH_STAGGERED = 4
};

/**
 * @brief DAG node for triple-row DP with 5 cases
 */
struct DPNode {
    // State: (i, s1, l1, s2, l2, s3, l3)
    int site;      // Current site i
    int s[3];      // Last placed cell index in each row
    int l[3];      // Distance from site to cell's left boundary
    
    // For each of 5 cases
    struct CaseInfo {
        int benefit;           // Maximum accumulated benefit
        DPNode* prev_node;     // Previous node in optimal path
        int prev_case;         // Which case in previous node
        int vdd_count;         // VDD staples so far
        int vss_count;         // VSS staples so far
        bool valid;            // Is this case valid?
        
        CaseInfo() : benefit(-1000000), prev_node(nullptr), prev_case(-1), 
                     vdd_count(0), vss_count(0), valid(false) {}
    } cases[5];
    
    // Constructor
    DPNode(int i, int s1, int l1, int s2, int l2, int s3, int l3) : site(i) {
        s[0] = s1; s[1] = s2; s[2] = s3;
        l[0] = l1; l[1] = l2; l[2] = l3;
    }
};

/**
 * @brief Compact state for hash lookup (Section 3.1 of paper)
 */
struct CompactState {
    int site;
    int a[3];  // Cell index offset from initial
    int b[3];  // Proxy for displacement
    
    bool operator==(const CompactState& other) const {
        if (site != other.site) return false;
        for (int i = 0; i < 3; i++) {
            if (a[i] != other.a[i] || b[i] != other.b[i]) return false;
        }
        return true;
    }
};

/**
 * @brief Hash function for CompactState
 */
struct CompactStateHasher {
    std::size_t operator()(const CompactState& state) const {
        std::size_t hash = std::hash<int>()(state.site);
        for (int i = 0; i < 3; i++) {
            hash ^= std::hash<int>()(state.a[i]) << ((i * 2 + 1) * 4);
            hash ^= std::hash<int>()(state.b[i]) << ((i * 2 + 2) * 4);
        }
        return hash;
    }
};

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
     * @brief Solve triple-row optimization problem
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
    std::vector<DPNode*> all_nodes;
    std::unordered_map<CompactState, DPNode*, CompactStateHasher> node_lookup;
    
    // State tracking
    int initial_s[3];  // Initial cell indices for compact encoding
    std::vector<Staple> prev_staples;  // Staples from previous triple-rows
    int row_start;  // Starting row of current triple-row problem
    
    // Statistics
    DPNode* best_final_node;
    int best_final_case;
    size_t nodes_created;
    
    // Core functions
    void generateExtensions(DPNode* current, 
                           const std::vector<std::vector<Cell*>>& cells_in_rows,
                           std::queue<DPNode*>& Q);
    
    bool isLegalExtension(DPNode* current, bool place_cell[3],
                         const std::vector<std::vector<Cell*>>& cells_in_rows,
                         int new_s[3], int new_l[3]);
    
    DPNode* createOrFindNode(int site, int s1, int l1, int s2, int l2, int s3, int l3);
    
    CompactState getCompactState(DPNode* node);
    
    void updateNodeCases(DPNode* source, DPNode* target,
                        const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    bool isValidCaseTransition(DPNode* source, int src_case,
                              DPNode* target, int tgt_case,
                              const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    std::vector<Staple> getStaplesForCase(int case_type, DPNode* node);
    
    bool hasStaggeringViolation(const std::vector<Staple>& new_staples);
    
    bool hasOverlapWithPins(const Staple& staple, DPNode* node,
                           const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    int calculateCellPosition(DPNode* node, int row_idx, int cell_idx,
                             const std::vector<std::vector<Cell*>>& cells_in_rows);
    
    int calculateCaseBenefit(int case_type, DPNode* node, int current_vdd, int current_vss);
    
    std::pair<int, int> getStapleTypesForCase(int case_type, int site);
    
    bool determineStapleType(int row_boundary);
    
    std::vector<Staple> extractBestSolution(
        const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start);
    
    void pruneOldNodes(int keep_after_site);
    
    void cleanup();
};

#endif // DP_SOLVER_HPP