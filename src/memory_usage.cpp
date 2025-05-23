#include "memory_usage.hpp"

size_t getCurrentMemoryUsage() {
#ifdef _WIN32
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), (PROCESS_MEMORY_COUNTERS*)&pmc, sizeof(pmc))) {
        return pmc.WorkingSetSize / (1024 * 1024); // Convert to MB
    }
    return 0;
#elif __linux__
    std::ifstream file("/proc/self/status");
    std::string line;
    while (std::getline(file, line)) {
        if (line.substr(0, 6) == "VmRSS:") {
            std::string memory_str = line.substr(6);
            size_t memory_kb = std::stoul(memory_str);
            return memory_kb / 1024; // Convert to MB
        }
    }
    return 0;
#elif __APPLE__
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    
    if (task_info(mach_task_self(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count) == KERN_SUCCESS) {
        return t_info.resident_size / (1024 * 1024); // Convert to MB
    }
    return 0;
#else
    return 0; // Fallback for other platforms
#endif
}