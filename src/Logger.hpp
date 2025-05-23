// Modified Logger.hpp

#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <chrono>
#include <ctime>
#include <iomanip>

/**
 * @brief A simple logger class for debugging
 */
class Logger {
private:
    static std::ofstream logFile;
    static bool initialized;
    static int indent;
    static int logLevel;  // Log level: 0=ERROR, 1=WARNING, 2=INFO, 3=DEBUG, 4=TRACE
    
public:
    // Log level enumeration
    enum LogLevel {
        ERROR = 0,
        WARNING = 1,
        INFO = 2,
        DEBUG = 3,
        TRACE = 4
    };

    static void init(const std::string& filename = "debug_log.txt", int level = INFO) {
        if (!initialized) {
            logFile.open(filename, std::ios::out | std::ios::trunc);
            initialized = true;
            indent = 0;
            logLevel = level;
            log(INFO, "Logger initialized with level " + std::to_string(logLevel));
        }
    }
    
    static void close() {
        if (initialized) {
            logFile.close();
            initialized = false;
        }
    }
    
    static void increaseIndent() {
        indent += 2;
    }
    
    static void decreaseIndent() {
        indent = std::max(0, indent - 2);
    }
    
    static void setLogLevel(int level) {
        logLevel = level;
        log(INFO, "Log level set to " + std::to_string(logLevel));
    }
    
    // New method with log level parameter
    template<typename T>
    static void log(LogLevel level, const T& message) {
        if (!initialized) {
            init();
        }
        
        // Only log if message level is less than or equal to current log level
        if (level <= logLevel) {
            // Get current time
            auto now = std::chrono::system_clock::now();
            auto time = std::chrono::system_clock::to_time_t(now);
            
            // Format timestamp
            std::stringstream timestamp;
            timestamp << std::put_time(std::localtime(&time), "%H:%M:%S");
            
            // Get level string
            std::string levelStr;
            switch(level) {
                case ERROR: levelStr = "ERROR"; break;
                case WARNING: levelStr = "WARN "; break;
                case INFO: levelStr = "INFO "; break;
                case DEBUG: levelStr = "DEBUG"; break;
                case TRACE: levelStr = "TRACE"; break;
                default: levelStr = "?????"; break;
            }
            
            // Create indentation
            std::string indentation(indent, ' ');
            
            // Write to log file
            logFile << timestamp.str() << " | " << levelStr << " | " 
                    << indentation << message << std::endl;
            
            // Also print critical messages to console
            // if (level <= WARNING) {
            //     std::cout << "LOG: " << levelStr << " | " << message << std::endl;
            // }
        }
    }
    
    // Legacy method for backward compatibility
    template<typename T>
    static void log(const T& message) {
        log(TRACE, message);
    }
};

