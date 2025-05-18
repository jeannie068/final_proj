/**
 * @file data_structure.hpp
 * @brief Data structures for power staple insertion optimization
 * 
 * This file contains all the data structures needed for the power staple insertion
 * optimization algorithm based on the DAG-based dynamic programming approach
 * described in the paper "Manufacturing-Aware Power Staple Insertion Optimization
 * by Enhanced Multi-Row Detailed Placement Refinement".
 */

#ifndef DATA_STRUCTURE_HPP
#define DATA_STRUCTURE_HPP

#include <vector>
#include <unordered_map>
#include <queue>
#include <memory>
#include <functional>

/**
 * @brief Basic chip layout information
 */
struct ChipInfo {
    int left_x;        // X-coordinate of the chip left boundary
    int bottom_y;      // Y-coordinate of the chip bottom boundary
    int right_x;       // X-coordinate of the chip right boundary
    int top_y;         // Y-coordinate of the chip top boundary
    int num_rows;      // Number of placement rows
    int row_height;    // Height of a placement row
    int site_width;    // Width of a placement site
    
    // Derived information
    int total_sites;   // Total number of sites in a row
    
    // Constructor
    ChipInfo() : left_x(0), bottom_y(0), right_x(0), top_y(0),
                num_rows(0), row_height(0), site_width(0), total_sites(0) {}
    
    // Calculate sites per row
    void calculateTotalSites() {
        total_sites = (right_x - left_x) / site_width;
    }
    
    // Get Y coordinate for a specific row
    int getRowY(int row_index) const {
        return bottom_y + row_index * row_height;
    }
};

/**
 * @brief Cell type definition
 */
struct CellType {
    int type_index;             // Cell type index
    int width;                  // Width of the cell type
    int height;                 // Height of the cell type
    int cellSiteWidth;         // Width in site units
    std::vector<int> pin_sites; // Sites with pins (relative to cell's left boundary)
    
    // Default constructor
    CellType() : type_index(0), width(0), height(0), cellSiteWidth(0) {}
    
    // Constructor
    CellType(int idx, int w, int h, int site_width) : type_index(idx), width(w), height(h), cellSiteWidth(0) {
        cellSiteWidth = w / site_width; // Calculate width in site units
    }
    
    // Check if a specific site has a pin
    bool hasPinAt(int site_idx, bool is_flipped) const {
        if (is_flipped) {
            // Adjust index for flipped cell
            int flipped_idx = cellSiteWidth - 1 - site_idx;
            for (int pin : pin_sites) {
                if (pin == flipped_idx) return true;
            }
            return false;
        } else {
            for (int pin : pin_sites) {
                if (pin == site_idx) return true;
            }
            return false;
        }
    }
};

/**
 * @brief Cell information
 */
struct Cell {
    int cell_index;         // Cell index
    int type_index;         // Cell type index
    int initial_x;          // Initial X-coordinate of the cell left boundary
    int initial_y;          // Initial Y-coordinate of the cell bottom boundary
    int max_displacement;   // Maximum allowed displacement in site units
    
    // Refined placement information
    int current_x;          // Current X-coordinate in the refined placement
    int current_y;          // Current Y-coordinate in the refined placement
    bool is_flipped;        // Whether the cell is horizontally flipped
    
    // Constructor
    Cell(int idx, int type, int x, int y, int max_disp) :
        cell_index(idx), type_index(type), initial_x(x), initial_y(y),
        max_displacement(max_disp), current_x(x), current_y(y), is_flipped(false) {}
    
    // Get row index
    int getRowIndex(int row_height) const {
        return initial_y / row_height;
    }
    
    // Get initial site index
    int getInitialSite(int site_width) const {
        return initial_x / site_width;
    }
    
    // Check if a position is within maximum displacement
    bool isWithinDisplacement(int x, int site_width) const {
        int displacement = std::abs(x - initial_x) / site_width;
        return displacement <= max_displacement;
    }
};

/**
 * @brief Staple information
 */
struct Staple {
    int x;          // X-coordinate of the staple left boundary
    int y;          // Y-coordinate of the staple bottom boundary
    bool is_vdd;    // true for VDD staple, false for VSS staple
    
    // Constructor
    Staple(int x_pos, int y_pos, bool vdd) : x(x_pos), y(y_pos), is_vdd(vdd) {}
    
    // Equal operator for comparing staples
    bool operator==(const Staple& other) const {
        return x == other.x && y == other.y;
    }
};

/**
 * @brief Represents the state of a triple-row optimization problem
 */
struct TripleRowState {
    // Row information for triple-row problem
    int row_start;          // Start row index
    int row_end;            // End row index
    
    // Cell information for each row
    std::vector<Cell*> cells_in_row[3];   // Cells in each row
    
    // Staple constraints
    std::vector<Staple> prev_staples;     // Staples inserted in previous rounds
    
    // Result
    std::vector<Staple> inserted_staples;  // Staples inserted in this round
    
    // Constructor
    TripleRowState(int start, int end) : row_start(start), row_end(end) {}
};

/**
 * @brief Compact state encoding for memory efficiency
 * 
 * This uses a relative encoding to reduce state space and memory usage
 */
struct CompactState {
    int site;               // Current site
    
    // For each row (3 rows)
    int cell_offset[3];     // Relative cell index from initial placement
    int displacement[3];    // Displacement from initial position
    bool is_flipped[3];     // Cell flipping status
    
    // Constructor
    CompactState() : site(0) {
        for (int i = 0; i < 3; i++) {
            cell_offset[i] = 0;
            displacement[i] = 0;
            is_flipped[i] = false;
        }
    }
    
    // Compare operator for map/set
    bool operator==(const CompactState& other) const {
        if (site != other.site) return false;
        
        for (int i = 0; i < 3; i++) {
            if (cell_offset[i] != other.cell_offset[i] ||
                displacement[i] != other.displacement[i] ||
                is_flipped[i] != other.is_flipped[i])
                return false;
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
            hash ^= std::hash<int>()(state.cell_offset[i]) << (i * 10);
            hash ^= std::hash<int>()(state.displacement[i]) << (i * 10 + 3);
            hash ^= std::hash<bool>()(state.is_flipped[i]) << (i * 10 + 6);
        }
        
        return hash;
    }
};

/**
 * @brief Staple insertion cases as defined in the paper
 */
enum StapleCase {
    NO_STAPLE = 0,      // Case 1: No staple inserted
    R1_R2_STAPLE = 1,   // Case 2: Staple between R1 and R2
    R2_R3_STAPLE = 2,   // Case 3: Staple between R2 and R3
    BOTH_STAPLES = 3,   // Case 4: Both staples (R1-R2 and R2-R3)
};

/**
 * @brief DAG node for dynamic programming
 */
struct DPNode {
    CompactState state;     // Current state
    
    // For each staple insertion case (4 cases simplify than paper)
    int benefit[4];         // Max accumulated staple benefit
    DPNode* prev_node[4];   // Previous node pointer
    int case_from_prev[4];  // Case in previous node
    
    // For staple balance constraint
    int vdd_staples[4];     // VDD staple count
    int vss_staples[4];     // VSS staple count
    
    // Constructor
    DPNode(const CompactState& s) : state(s) {
        for (int i = 0; i < 4; i++) {
            benefit[i] = (i == 0) ? 0 : -1000000; // Initialize case 0 with 0, others with large negative
            prev_node[i] = nullptr;
            case_from_prev[i] = -1;
            vdd_staples[i] = 0;
            vss_staples[i] = 0;
        }
    }
};

/**
 * @brief Extension object representing a way to extend a partial solution
 */
struct Extension {
    int site;                   // Target site
    int row_idx;                // Row being extended
    int new_cell_idx[3];        // New cell index for each row
    int new_displacement[3];    // New displacement for each row
    bool new_is_flipped[3];     // New flipping status for each row
    
    // Calculated staple benefits
    int staple_benefit;         // Staple benefit of this extension
    int vdd_count;              // VDD staples added
    int vss_count;              // VSS staples added
    
    // Constructor
    Extension(int s, int r) : site(s), row_idx(r), staple_benefit(0), vdd_count(0), vss_count(0) {
        for (int i = 0; i < 3; i++) {
            new_cell_idx[i] = -1;
            new_displacement[i] = 0;
            new_is_flipped[i] = false;
        }
    }
};

/**
 * @brief Solution structure containing refined placement and inserted staples
 */
struct Solution {
    std::vector<Cell> refined_cells;     // Refined cell placement
    std::vector<Staple> inserted_staples; // Inserted staples
    
    // Statistics
    int total_staples;
    int vdd_staples;
    int vss_staples;
    double staple_ratio;
    
    // Constructor
    Solution() : total_staples(0), vdd_staples(0), vss_staples(0), staple_ratio(0.0) {}
    
    // Update statistics
    void updateStats() {
        vdd_staples = 0;
        vss_staples = 0;
        
        for (const Staple& staple : inserted_staples) {
            if (staple.is_vdd) vdd_staples++;
            else vss_staples++;
        }
        
        total_staples = vdd_staples + vss_staples;
        
        if (vdd_staples > 0 && vss_staples > 0) {
            double max_count = std::max(vdd_staples, vss_staples);
            double min_count = std::min(vdd_staples, vss_staples);
            staple_ratio = max_count / min_count;
        } else {
            staple_ratio = 0.0; // Cannot calculate ratio if either is zero
        }
    }
    
    // Check if staple balance constraint is satisfied
    bool isBalanceConstraintSatisfied() const {
        if (vdd_staples == 0 || vss_staples == 0) return true; // No staples of one type
        return staple_ratio <= 1.1;
    }
};

/**
 * @brief Global algorithm parameters
 */
struct AlgorithmParams {
    double balance_factor;   // Beta factor for staple balance encouragement
    int verbose_level;       // Verbosity level for debugging
    
    // Constructor with default values
    AlgorithmParams() : balance_factor(0.4), verbose_level(0) {}
};

#endif // DATA_STRUCTURE_HPP