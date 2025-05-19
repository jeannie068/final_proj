/**
 * @file main.cpp
 * @brief Main program for power staple insertion optimization
 * 
 * This file contains the main function that ties together the parser and optimizer
 * to solve the power staple insertion optimization problem.
 */

#include "data_structure/data_structure.hpp"
#include "parser/parser.hpp"
#include "alg/optimizer.hpp"
#include <iostream>
#include <string>
#include <chrono>
#include <cstdlib>
#include <memory>

/**
 * @brief Main function
 * @param argc Number of command line arguments
 * @param argv Command line arguments
 * @return Exit code
 */
int main(int argc, char* argv[]) {
    // Check command line arguments
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <input file> <output file>" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::string input_file = argv[1];
    std::string output_file = argv[2];
    
    // Start timer
    auto start_time = std::chrono::high_resolution_clock::now();
    
    // Parse input file
    ChipInfo chip_info;
    std::vector<CellType> cell_types;
    std::vector<Cell> cells;
    
    std::cout << "Parsing input file: " << input_file << std::endl;
    if (!Parser::parseInputFile(input_file, chip_info, cell_types, cells)) {
        std::cerr << "Failed to parse input file" << std::endl;
        return EXIT_FAILURE;
    }
    
    // Print summary of parsed data
    Parser::printSummary(chip_info, cell_types, cells);
    
    // Create algorithm parameters
    AlgorithmParams params;
    params.balance_factor = 0.4;  // Beta factor for staple balance
    params.verbose_level = 1;     // Verbosity level
    
    // Create optimizer
    std::cout << "Running optimization..." << std::endl;
    Optimizer optimizer(chip_info, cell_types, cells, params);
    
    // Run optimization
    Solution solution = optimizer.run();
    
    // Print solution summary
    Parser::printSolutionSummary(solution);
    
    // Check if solution is valid
    if (!solution.isBalanceConstraintSatisfied()) {
        std::cerr << "Warning: Staple balance constraint not satisfied" << std::endl;
    }
    
    // Write output file
    std::cout << "Writing output file: " << output_file << std::endl;
    if (!Parser::writeOutputFile(output_file, solution)) {
        std::cerr << "Failed to write output file" << std::endl;
        return EXIT_FAILURE;
    }
    
    // End timer
    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
    
    std::cout << "Optimization completed in " << duration.count() << " seconds" << std::endl;
    std::cout << "Total staples inserted: " << solution.total_staples << std::endl;
    std::cout << "VDD staples: " << solution.vdd_staples << ", VSS staples: " << solution.vss_staples << std::endl;
    std::cout << "Staple ratio: " << solution.staple_ratio << std::endl;
    
    return EXIT_SUCCESS;
}