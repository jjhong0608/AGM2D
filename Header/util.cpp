//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include "util.h"

auto AGM::isclose(double x, double y, double eps) -> bool {
    return std::fabs(x - y) < eps;
}

void AGM::printError(const std::string &function_name) {
    std::cout << std::endl << function_name << std::endl;
    exit(1);
}

void AGM::printError(const char *function_name, const char *fmt, ...) {
    char buf[256] = {0,};
    va_list ap;

    printf("\n");
    printf("Fatal error has occur in %s\n", function_name);
    sprintf(buf, "Massage: ");

    va_start(ap, fmt);
    vsprintf(buf + strlen(buf), fmt, ap);
    va_end(ap);

    puts(buf);
    exit(1);
}

auto AGM::iszero(double x, double eps) -> bool {
    return std::fabs(x) < eps;
}

auto AGM::sgn(double d) -> double {
    return (d > ZEROVALUE) - (d < ZEROVALUE);
}

auto AGM::ispositive(double d) -> bool {
    return d > NEARZERO;
}

auto AGM::isnegative(double d) -> bool {
    return d < -NEARZERO;
}

auto AGM::ispositivezero(double d) -> bool {
    return ispositive(d) || iszero(d);
}

auto AGM::isnegativezero(double d) -> bool {
    return isnegative(d) || iszero(d);
}

auto AGM::min(double d0, double d1) -> double {
    if (d0 < d1) {
        return d0;
    } else {
        return d1;
    }
}

auto AGM::max(double d0, double d1) -> double {
    if (d0 < d1) {
        return d1;
    } else {
        return d0;
    }
}