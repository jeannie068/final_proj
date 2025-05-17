/**
 * @file parser.hpp
 * @brief Input/output parsing functions for power staple insertion optimization
 * 
 * This file contains functions for parsing input files and writing output files
 * according to the specified formats for the power staple insertion optimization project.
 */

#ifndef PARSER_HPP
#define PARSER_HPP

#include "../data_structure/data_structure.hpp"
#include <string>
#include <fstream>
#include <sstream>
#include <limits>
#include <vector>
#include <iostream>
#include <algorithm>

/**
 * @brief Main parser class for input/output handling
 */
class Parser {
public:
    /**
     * @brief Parse input file and fill data structures
     * @param filename Input file path
     * @param chip_info ChipInfo structure to fill
     * @param cell_types Vector of CellType to fill
     * @param cells Vector of Cell to fill
     * @return true if parsing was successful, false otherwise
     */
    static bool parseInputFile(const std::string& filename, 
                              ChipInfo& chip_info,
                              std::vector<CellType>& cell_types,
                              std::vector<Cell>& cells);
    
    /**
     * @brief Write solution to output file
     * @param filename Output file path
     * @param solution Solution structure with refined placement and staples
     * @return true if writing was successful, false otherwise
     */
    static bool writeOutputFile(const std::string& filename,
                               const Solution& solution);
    
    /**
     * @brief Print summary of parsed data for verification
     * @param chip_info ChipInfo structure
     * @param cell_types Vector of CellType
     * @param cells Vector of Cell
     */
    static void printSummary(const ChipInfo& chip_info,
                            const std::vector<CellType>& cell_types,
                            const std::vector<Cell>& cells);
    
    /**
     * @brief Print solution summary
     * @param solution Solution structure with refined placement and staples
     */
    static void printSolutionSummary(const Solution& solution);
};

/**
 * @brief Parse input file and fill data structures
 */
bool Parser::parseInputFile(const std::string& filename, 
                           ChipInfo& chip_info,
                           std::vector<CellType>& cell_types,
                           std::vector<Cell>& cells) {
    std::ifstream infile(filename);
    if (!infile.is_open()) {
        std::cerr << "Error: Could not open input file " << filename << std::endl;
        return false;
    }
    
    // Parse first part: chip info (4 lines)
    
    // Line 1: chip boundaries
    if (!(infile >> chip_info.left_x >> chip_info.bottom_y >> 
          chip_info.right_x >> chip_info.top_y)) {
        std::cerr << "Error: Failed to parse chip boundaries" << std::endl;
        return false;
    }
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    // Line 2: row info and site width
    if (!(infile >> chip_info.num_rows >> chip_info.row_height >> 
          chip_info.site_width)) {
        std::cerr << "Error: Failed to parse row info" << std::endl;
        return false;
    }
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    // Line 3: number of cell types
    int num_cell_types;
    if (!(infile >> num_cell_types)) {
        std::cerr << "Error: Failed to parse number of cell types" << std::endl;
        return false;
    }
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    // Line 4: number of cells
    int num_cells;
    if (!(infile >> num_cells)) {
        std::cerr << "Error: Failed to parse number of cells" << std::endl;
        return false;
    }
    infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    
    // Calculate total sites
    chip_info.calculateTotalSites();
    
    // Parse second part: cell types
    cell_types.clear();
    cell_types.resize(num_cell_types);
    
    for (int i = 0; i < num_cell_types; i++) {
        std::string line;
        if (!std::getline(infile, line)) {
            std::cerr << "Error: Failed to read cell type line " << i << std::endl;
            return false;
        }
        
        std::istringstream iss(line);
        int type_idx, width, height;
        if (!(iss >> type_idx >> width >> height)) {
            std::cerr << "Error: Failed to parse cell type info for type " << i << std::endl;
            return false;
        }
        
        // Create cell type with index, width, and height
        cell_types[type_idx] = CellType(type_idx, width, height, chip_info.site_width);
        
        // Read pin sites
        int pin_site;
        while (iss >> pin_site) {
            cell_types[type_idx].pin_sites.push_back(pin_site);
        }
    }
    
    // Parse third part: cells
    cells.clear();
    cells.reserve(num_cells);
    
    for (int i = 0; i < num_cells; i++) {
        int cell_idx, type_idx, x, y, max_disp;
        if (!(infile >> cell_idx >> type_idx >> x >> y >> max_disp)) {
            std::cerr << "Error: Failed to parse cell info for cell " << i << std::endl;
            return false;
        }
        infile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        
        // Verify that type_idx is valid
        if (type_idx < 0 || type_idx >= static_cast<int>(cell_types.size())) {
            std::cerr << "Error: Cell " << cell_idx << " has invalid type index " << type_idx << std::endl;
            return false;
        }
        
        // Create cell with the parsed information
        cells.push_back(Cell(cell_idx, type_idx, x, y, max_disp));
    }
    
    infile.close();
    return true;
}

/**
 * @brief Write solution to output file
 */
bool Parser::writeOutputFile(const std::string& filename,
                            const Solution& solution) {
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Error: Could not open output file " << filename << std::endl;
        return false;
    }
    
    // Write first part: refined cell placement
    for (const Cell& cell : solution.refined_cells) {
        outfile << cell.cell_index << " " 
                << cell.current_x << " " 
                << cell.current_y << " "
                << (cell.is_flipped ? 1 : 0) << std::endl;
    }
    
    // Write second part: inserted staples
    for (const Staple& staple : solution.inserted_staples) {
        outfile << staple.x << " " << staple.y << std::endl;
    }
    
    outfile.close();
    return true;
}

/**
 * @brief Print summary of parsed data for verification
 */
void Parser::printSummary(const ChipInfo& chip_info,
                         const std::vector<CellType>& cell_types,
                         const std::vector<Cell>& cells) {
    std::cout << "===== Input Summary =====" << std::endl;
    
    // Chip info
    std::cout << "Chip boundaries: (" << chip_info.left_x << ", " << chip_info.bottom_y << ") to ("
              << chip_info.right_x << ", " << chip_info.top_y << ")" << std::endl;
    std::cout << "Rows: " << chip_info.num_rows << ", Row height: " << chip_info.row_height
              << ", Site width: " << chip_info.site_width << std::endl;
    std::cout << "Total sites per row: " << chip_info.total_sites << std::endl;
    
    // Cell types
    std::cout << "Cell types: " << cell_types.size() << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), cell_types.size()); i++) {
        const CellType& type = cell_types[i];
        std::cout << "  Type " << type.type_index << ": " << type.width << " x " << type.height
                  << ", Pins at sites:";
        for (int pin : type.pin_sites) {
            std::cout << " " << pin;
        }
        std::cout << std::endl;
    }
    if (cell_types.size() > 5) {
        std::cout << "  ... and " << (cell_types.size() - 5) << " more types" << std::endl;
    }
    
    // Cells
    std::cout << "Cells: " << cells.size() << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), cells.size()); i++) {
        const Cell& cell = cells[i];
        std::cout << "  Cell " << cell.cell_index << ": Type " << cell.type_index
                  << ", Position (" << cell.initial_x << ", " << cell.initial_y
                  << "), Max displacement: " << cell.max_displacement << std::endl;
    }
    if (cells.size() > 5) {
        std::cout << "  ... and " << (cells.size() - 5) << " more cells" << std::endl;
    }
    
    std::cout << "=========================" << std::endl;
}

/**
 * @brief Print solution summary
 */
void Parser::printSolutionSummary(const Solution& solution) {
    std::cout << "===== Solution Summary =====" << std::endl;
    
    // Refined cells
    std::cout << "Refined cells: " << solution.refined_cells.size() << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), solution.refined_cells.size()); i++) {
        const Cell& cell = solution.refined_cells[i];
        std::cout << "  Cell " << cell.cell_index << ": Position ("
                  << cell.current_x << ", " << cell.current_y
                  << "), Flipped: " << (cell.is_flipped ? "Yes" : "No") << std::endl;
    }
    if (solution.refined_cells.size() > 5) {
        std::cout << "  ... and " << (solution.refined_cells.size() - 5) << " more cells" << std::endl;
    }
    
    // Staples
    std::cout << "Inserted staples: " << solution.inserted_staples.size() << std::endl;
    for (size_t i = 0; i < std::min(size_t(5), solution.inserted_staples.size()); i++) {
        const Staple& staple = solution.inserted_staples[i];
        std::cout << "  Staple at (" << staple.x << ", " << staple.y << "), Type: "
                  << (staple.is_vdd ? "VDD" : "VSS") << std::endl;
    }
    if (solution.inserted_staples.size() > 5) {
        std::cout << "  ... and " << (solution.inserted_staples.size() - 5) << " more staples" << std::endl;
    }
    
    // Statistics
    std::cout << "Total staples: " << solution.total_staples << std::endl;
    std::cout << "VDD staples: " << solution.vdd_staples << std::endl;
    std::cout << "VSS staples: " << solution.vss_staples << std::endl;
    std::cout << "Staple ratio: " << solution.staple_ratio << std::endl;
    std::cout << "Balance constraint satisfied: " << (solution.isBalanceConstraintSatisfied() ? "Yes" : "No") << std::endl;
    
    std::cout << "============================" << std::endl;
}

#endif // PARSER_HPP