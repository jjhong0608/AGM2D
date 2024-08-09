//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_UTIL_H
#define AGM_UTIL_H

#include <omp.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "indicators/block_progress_bar.hpp"
#include "indicators/cursor_control.hpp"

constexpr auto UNITVALUE{1.0000000000000000E0};
constexpr auto HALFVALUE{5.0000000000000000E-1};
constexpr auto ZEROVALUE{0.0000000000000000E0};
constexpr auto NEARZERO{1.0000000000000000E-10};

namespace AGM {
auto isclose(double x, double y, double eps = NEARZERO) -> bool;

void printError(const std::string &function_name);

void printError(const char *function_name, const char *fmt, ...);

auto iszero(double x, double eps = NEARZERO) -> bool;

auto sgn(double d) -> double;

auto ispositive(double d) -> bool;

auto isnegative(double d) -> bool;

auto ispositivezero(double d) -> bool;

auto isnegativezero(double d) -> bool;

[[nodiscard]] auto min(double d0, double d1) -> double;

[[nodiscard]] auto max(double d0, double d1) -> double;

enum EWNS { E,
            W,
            N,
            S,
            EN,
            ES,
            WN,
            WS,
            NE,
            NW,
            SE,
            SW };

static std::unordered_map<std::string, int> valueMap{{"sol", 0}, {"phi", 1}, {"rhs", 2}, {"bdv", 3}, {"dx", 4}, {"dy", 5}, {"dxx", 6}, {"dyy", 7}};
}// namespace AGM

#endif//AGM_UTIL_H