/**
 * @file dp_solver.hpp
 * @brief Enhanced dynamic programming solver for manufacturing-aware power staple insertion
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
#include <iomanip>

/**
 * @brief Enhanced staple case definitions following paper Figure 7
 */
enum EnhancedStapleCase {
    CASE_1_NO_STAPLE = 0,        // No staple inserted
    CASE_2_R1_R2_ONLY = 1,       // Only R1-R2 staple
    CASE_3_R2_R3_ONLY = 2,       // Only R2-R3 staple  
    CASE_4_BOTH_ALIGNED = 3,     // Both R1-R2 and R2-R3 staples, aligned
    CASE_5_BOTH_STAGGERED = 4    // Both staples, but staggered (special manufacturing case)
};

/**
 * @brief Enhanced compact state with proper displacement tracking
 */
struct EnhancedCompactState {
    int site;
    int a[3];           // Cell offset from initial placement
    int b[3];           // Actual displacement of cells (not l values)
    bool flipped[3];    // Flipping status of last placed cells
    
    bool operator==(const EnhancedCompactState& other) const {
        if (site != other.site) return false;
        for (int i = 0; i < 3; i++) {
            if (a[i] != other.a[i] || b[i] != other.b[i] || flipped[i] != other.flipped[i]) 
                return false;
        }
        return true;
    }
};

/**
 * @brief Hash function for EnhancedCompactState
 */
struct EnhancedCompactStateHasher {
    std::size_t operator()(const EnhancedCompactState& state) const {
        std::size_t hash = std::hash<int>()(state.site);
        for (int i = 0; i < 3; i++) {
            hash ^= std::hash<int>()(state.a[i]) << ((i * 3 + 1) * 4);
            hash ^= std::hash<int>()(state.b[i]) << ((i * 3 + 2) * 4);
            hash ^= std::hash<bool>()(state.flipped[i]) << ((i * 3 + 3) * 4);
        }
        return hash;
    }
};

/**
 * @brief Enhanced DP node with proper displacement tracking
 */
struct EnhancedDPNode {
    int site;           // Current site i
    int s[3];          // Last placed cell index in each row
    int l[3];          // Distance from site to cell's left boundary
    int displacement[3]; // Actual displacement of last placed cells
    bool is_flipped[3]; // Flipping status of last placed cells
    
    // Five staple insertion cases as per Figure 7 in the paper
    static const int NUM_CASES = 5;
    
    // For each case, track:
    int benefit[NUM_CASES];         // Maximum accumulated staple benefit
    EnhancedDPNode* prev_node[NUM_CASES];   // Previous node pointer
    int prev_case[NUM_CASES];       // Case in previous node
    
    // Staple balance tracking
    int vdd_staples[NUM_CASES];
    int vss_staples[NUM_CASES];
    
    // Detailed staple information for consistency
    std::vector<Staple> staples_this_step[NUM_CASES];
    
    // Constructor
    EnhancedDPNode(int i, int s1, int l1, int s2, int l2, int s3, int l3) 
        : site(i) {
        s[0] = s1; s[1] = s2; s[2] = s3;
        l[0] = l1; l[1] = l2; l[2] = l3;
        
        // Initialize displacement and flipping
        for (int j = 0; j < 3; j++) {
            displacement[j] = 0;
            is_flipped[j] = false;
        }
        
        for (int c = 0; c < NUM_CASES; c++) {
            benefit[c] = (c == 0) ? 0 : -1000000; // Case 1 (no staple) starts at 0
            prev_node[c] = nullptr;
            prev_case[c] = -1;
            vdd_staples[c] = 0;
            vss_staples[c] = 0;
        }
    }
};

/**
 * @brief Staple insertion information for validation
 */
struct StapleInsertionInfo {
    int x, y;
    bool is_vdd;
    int between_rows;  // 0: R1-R2, 1: R2-R3
    EnhancedStapleCase case_type;
};

/**
 * @brief Enhanced dynamic programming solver for triple-row optimization
 */
class DPSolver {
public:
    DPSolver(const ChipInfo& chip_info, 
             const std::vector<CellType>& cell_types,
             const AlgorithmParams& params = AlgorithmParams());
    
    ~DPSolver();
    
    /**
     * @brief Solve triple-row optimization problem using enhanced MATRO algorithm
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
    
    // Enhanced DAG structure
    std::queue<EnhancedDPNode*> node_queue;
    std::unordered_map<EnhancedCompactState, EnhancedDPNode*, EnhancedCompactStateHasher> node_lookup;
    std::vector<EnhancedDPNode*> all_nodes;
    
    // Tracking
    EnhancedDPNode* best_final_node;
    int best_final_case;
    size_t nodes_created;
    
    // Initial cell positions for compact encoding
    int initial_s[3];
    int initial_cell_x[3][1000];  // Store initial x positions for displacement calculation
    
    // Enhanced pruning strategy
    static const int BEAM_WIDTH = 800;        // Keep top nodes per site
    static const int PRUNING_FREQUENCY = 5;   // Prune every 5 sites
    static const int MAX_TOTAL_NODES = 100000; // Global node limit
    
    /**
     * @brief Enhanced legality checking with staggering validation
     */
    bool isLegalExtension(EnhancedDPNode* source, int target_site,
                         int new_s[3], int new_l[3], int new_disp[3], bool new_flipped[3],
                         bool placement_mask[3],
                         const std::vector<std::vector<Cell*>>& cells_in_rows,
                         const std::vector<Staple>& prev_staples, int row_start);
    
    /**
     * @brief Enhanced benefit update with early staggering check
     */
    void updateBenefitEnhanced(EnhancedDPNode* source, EnhancedDPNode* target,
                              const std::vector<std::vector<Cell*>>& cells_in_rows,
                              const std::vector<Staple>& prev_staples, int row_start);
    
    /**
     * @brief Calculate staples for a specific case
     */
    std::vector<Staple> calculateStaplesForCase(EnhancedStapleCase case_type, 
                                               EnhancedDPNode* node, int row_start);
    
    /**
     * @brief Enhanced staggering violation check
     */
    bool hasStaggeringViolationEnhanced(const std::vector<Staple>& new_staples,
                                       const std::vector<Staple>& prev_staples,
                                       int current_site, int row_start);
    
    /**
     * @brief Detailed pin and staple overlap checking
     */
    bool hasOverlapWithPinsOrStaples(const Staple& new_staple,
                                    const std::vector<std::vector<Cell*>>& cells_in_rows,
                                    const std::vector<Staple>& existing_staples,
                                    int row_start);
    
    /**
     * @brief Check if a staple overlaps with cell pins
     */
    bool overlapWithPins(const Staple& staple, const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start);
    
    /**
     * @brief Check if a staple overlaps with other staples
     */
    bool overlapWithStaples(const Staple& staple, const std::vector<Staple>& existing_staples);
    
    /**
     * @brief Intelligent node pruning strategy
     */
    void pruneNodesIntelligently(int current_site);
    
    /**
     * @brief Create or find enhanced node
     */
    EnhancedDPNode* createOrFindEnhancedNode(int site, int s1, int l1, int s2, int l2, int s3, int l3,
                                           int disp1, int disp2, int disp3,
                                           bool flip1, bool flip2, bool flip3);
    
    /**
     * @brief Enhanced solution extraction
     */
    std::vector<Staple> extractBestSolutionEnhanced(
        const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start);
    
    /**
     * @brief Calculate proper case benefit based on paper
     */
    int calculateEnhancedCaseBenefit(EnhancedStapleCase case_type, int current_vdd, int current_vss);
    
    /**
     * @brief Determine staple type with advanced balancing
     */
    bool determineStapleTypeAdvanced(int between_rows, int current_vdd, int current_vss, int site);
    
    /**
     * @brief Check if row is VDD or VSS
     */
    bool isVDDRow(int row_idx) const;
    
    /**
     * @brief Generate fallback solution
     */
    std::vector<Staple> generateEnhancedFallback(
        const std::vector<std::vector<Cell*>>& cells_in_rows, int row_start);
    
    /**
     * @brief Clean up allocated memory
     */
    void cleanup();
};

#endif // DP_SOLVER_HPP