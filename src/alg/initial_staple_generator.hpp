/**
 * @file initial_staple_generator.hpp
 * @brief Generate initial staples without cell refinement
 */

#ifndef INITIAL_STAPLE_GENERATOR_HPP
#define INITIAL_STAPLE_GENERATOR_HPP

#include "../data_structure/data_structure.hpp"
#include "../Logger.hpp"
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

class InitialStapleGenerator {
private:
    const ChipInfo& chip_info;
    const std::vector<CellType>& cell_types;
    const std::vector<Cell>& cells;
    
    // Track occupied sites by pins
    std::vector<std::vector<bool>> pin_occupied; // [row][site]
    
    // Track existing staples for staggering check
    std::set<std::pair<int, int>> staple_positions; // (site, row)
    
public:
    InitialStapleGenerator(const ChipInfo& chip_info,
                          const std::vector<CellType>& cell_types,
                          const std::vector<Cell>& cells)
        : chip_info(chip_info), cell_types(cell_types), cells(cells) {
        
        // Initialize pin occupation map
        pin_occupied.resize(chip_info.num_rows);
        for (int row = 0; row < chip_info.num_rows; row++) {
            pin_occupied[row].resize(chip_info.total_sites, false);
        }
        
        // Mark sites occupied by pins
        markPinOccupiedSites();
    }
    
    /**
     * @brief Generate initial staples without moving any cells
     * @return Vector of valid staples
     */
    std::vector<Staple> generateInitialStaples() {
        std::vector<Staple> staples;
        int vdd_count = 0;
        int vss_count = 0;
        
        // Process all sites
        for (int site = 0; site < chip_info.total_sites; site++) {
            // Try to place staples at this site
            std::vector<Staple> candidates = generateCandidatesAtSite(site, vdd_count, vss_count);
            
            // Add valid candidates
            for (const Staple& candidate : candidates) {
                if (isValidStaplePosition(candidate, staples)) {
                    staples.push_back(candidate);
                    staple_positions.insert({candidate.x / chip_info.site_width, 
                                           (candidate.y - chip_info.bottom_y) / chip_info.row_height});
                    
                    if (candidate.is_vdd) vdd_count++;
                    else vss_count++;
                    
                    // Early termination if balance is getting bad
                    if (vdd_count > 0 && vss_count > 0) {
                        double ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                                      static_cast<double>(std::min(vdd_count, vss_count));
                        if (ratio > 1.08) { // Leave some margin
                            // Skip this site to maintain balance
                            break;
                        }
                    }
                }
            }
        }
        
        std::cout << "Initial staple generation completed:" << std::endl;
        std::cout << "  Total staples: " << staples.size() << std::endl;
        std::cout << "  VDD staples: " << vdd_count << std::endl;
        std::cout << "  VSS staples: " << vss_count << std::endl;
        if (vdd_count > 0 && vss_count > 0) {
            double ratio = static_cast<double>(std::max(vdd_count, vss_count)) / 
                          static_cast<double>(std::min(vdd_count, vss_count));
            std::cout << "  Balance ratio: " << ratio << std::endl;
        }
        
        // Final validation
        std::cout << "Performing final validation..." << std::endl;
        for (const Staple& s : staples) {
            int s_site = s.x / chip_info.site_width;
            int s_row = (s.y - chip_info.bottom_y) / chip_info.row_height;
            
            // Check bounds
            if (s_row < 0 || s_row >= chip_info.num_rows) {
                std::cerr << "ERROR: Staple at (" << s.x << ", " << s.y << ") has invalid row boundary" << std::endl;
            }
            
            // Check pins
            if (pin_occupied[s_row][s_site] || pin_occupied[s_row+1][s_site]) {
                std::cerr << "ERROR: Staple at (" << s.x << ", " << s.y << ") overlaps with pins!" << std::endl;
            }
        }
        
        return staples;
    }
    
private:
    /**
     * @brief Mark all sites occupied by cell pins
     */
    void markPinOccupiedSites() {
        for (const Cell& cell : cells) {
            int row = cell.getRowIndex(chip_info.row_height);
            if (row < 0 || row >= chip_info.num_rows) continue;
            
            const CellType& cell_type = cell_types[cell.type_index];
            int cell_left_site = cell.initial_x / chip_info.site_width;
            
            // Debug: Check if cell type is valid
            if (cell.type_index < 0 || cell.type_index >= static_cast<int>(cell_types.size())) {
                std::cerr << "Warning: Invalid cell type index " << cell.type_index << std::endl;
                continue;
            }
            
            // Mark all pin sites as occupied
            for (int pin_site : cell_type.pin_sites) {
                int absolute_site = cell_left_site + pin_site;
                if (absolute_site >= 0 && absolute_site < chip_info.total_sites) {
                    pin_occupied[row][absolute_site] = true;
                }
            }
        }
        
        // Debug: Print some statistics
        int total_pins = 0;
        for (int row = 0; row < chip_info.num_rows; row++) {
            for (int site = 0; site < chip_info.total_sites; site++) {
                if (pin_occupied[row][site]) total_pins++;
            }
        }
        std::cout << "Marked " << total_pins << " sites as occupied by pins" << std::endl;
    }
    
    /**
     * @brief Generate candidate staples at a specific site
     */
    std::vector<Staple> generateCandidatesAtSite(int site, int current_vdd, int current_vss) {
        std::vector<Staple> candidates;
        int x = site * chip_info.site_width;
        
        // Check all possible staple positions (between adjacent rows)
        // A staple at row boundary 'row' connects rows (row) and (row+1)
        for (int row = 0; row < chip_info.num_rows-1; row++) {
            // Check if both adjacent rows are free at this site
            if (!pin_occupied[row][site] && !pin_occupied[row+1][site]) {
                int y = chip_info.getRowY(row);
                
                // Determine staple type based on balance
                bool should_be_vdd = determineStapleType(row, current_vdd, current_vss);
                
                candidates.push_back(Staple(x, y, should_be_vdd));
            }
        }
        
        return candidates;
    }
    
    /**
     * @brief Determine if a staple should be VDD or VSS based on balance
     */
    bool determineStapleType(int row, int current_vdd, int current_vss) {
        // First priority: ensure we have at least one of each type
        if (current_vdd == 0) return true;
        if (current_vss == 0) return false;
        
        // Second priority: maintain balance
        if (current_vdd < current_vss) return true;
        if (current_vss < current_vdd) return false;
        
        // If equal, alternate based on row for spatial distribution
        return (row % 2 == 0);
    }
    
    /**
     * @brief Check if a staple position is valid (no overlaps and no staggering violations)
     */
    bool isValidStaplePosition(const Staple& candidate, const std::vector<Staple>& existing_staples) {
        int cand_site = candidate.x / chip_info.site_width;
        int cand_row = (candidate.y - chip_info.bottom_y) / chip_info.row_height;
        
        // Check for overlaps with existing staples (same position)
        for (const Staple& existing : existing_staples) {
            if ((existing.x == candidate.x && existing.y == candidate.y) || 
                (existing.x == candidate.x && existing.y+chip_info.row_height == candidate.y)) {
                return false; // Exact position overlap
            }
        }
        
        // Check for staggering with existing staples
        for (const Staple& existing : existing_staples) {
            int exist_site = existing.x / chip_info.site_width;
            int exist_row = (existing.y - chip_info.bottom_y) / chip_info.row_height;
            
            if (exist_row == cand_row) {
                continue; // Same row, no staggering violation
            }

            // if (exist_row<cand_row) exist_row++;
            // else if (exist_row>cand_row) cand_row++;
            
            if (hasStaggeringViolation(cand_site, cand_row, exist_site, exist_row)) {
                return false;
            }
        }
        
        return true;
    }
    
    /**
     * @brief Check if two staples create a staggering violation
     * Based on TA's definition: diagonal staples without blocking staples
     */
    bool hasStaggeringViolation(int site1, int row1, int site2, int row2) {
        // Check if they are diagonal
        if (std::abs(site1 - site2) == 1 && std::abs(row1 - row2) == 2) {
            // Identify lower and upper staples
            int lower_site, lower_row, upper_site, upper_row;
            if (row1 < row2) {
                lower_site = site1; lower_row = row1;
                upper_site = site2; upper_row = row2;
            } else {
                lower_site = site2; lower_row = row2;
                upper_site = site1; upper_row = row1;
            }
            
            // Check for blocking staples
            bool has_blocker_above_lower = staple_positions.count({lower_site, lower_row + 2}) > 0;
            bool has_blocker_below_upper = staple_positions.count({upper_site, upper_row - 2}) > 0;

            // Logger::log(Logger::INFO, "Checking staggering between (" + 
            //             std::to_string(lower_site) + ", " + std::to_string(lower_row) + 
            //             ") and (" + std::to_string(upper_site) + ", " + std::to_string(upper_row) + 
            //             "): Blocker above lower: " + std::to_string(has_blocker_above_lower) + 
            //             ", Blocker below upper: " + std::to_string(has_blocker_below_upper));
            
            // Staggering occurs only if both blockers are missing
            return !has_blocker_above_lower && !has_blocker_below_upper;
        }
        else {
            return false; // Not diagonal, no staggering
        }
        
        
    }
};

/**
 * @brief Simple function to generate initial staples for testing
 * Can be called from Optimizer::run() instead of the complex DP
 */
std::vector<Staple> generateSimpleInitialSolution(const ChipInfo& chip_info,
                                                  const std::vector<CellType>& cell_types,
                                                  const std::vector<Cell>& cells) {
    InitialStapleGenerator generator(chip_info, cell_types, cells);
    return generator.generateInitialStaples();
}

#endif // INITIAL_STAPLE_GENERATOR_HPP