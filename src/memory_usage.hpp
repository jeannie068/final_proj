#pragma once

#ifdef _WIN32
#include <windows.h>
#include <psapi.h>
#elif __linux__
#include <unistd.h>
#include <fstream>
#include <string>
#elif __APPLE__
#include <mach/mach.h>
#endif

// 跨平台記憶體使用監控
size_t getCurrentMemoryUsage() ;