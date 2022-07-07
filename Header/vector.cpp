//
// Created by NIMS-JUNHONG on 2022/07/04.
//

#include "vector.h"

AGM::vector::vector() : std::vector<double>() {

}

AGM::vector::vector(const std::vector<double> &x) : std::vector<double>(x) {

}

AGM::vector::~vector() = default;

double AGM::vector::norm() {
    double sum{};
    for (const auto &i: *this) {
        sum += i * i;
    }
    return sqrt(sum);
}

double AGM::vector::dot(const AGM::vector &src) {
    if (size() != src.size()) printError("AGM::vector::dot", "The sizes of the two vectors do not match.");
    double sum{};
    for (int i = 0; i < size(); ++i) {
        sum += at(i) * src.at(i);
    }
    return sum;
}

AGM::vector AGM::vector::cross(const AGM::vector &src) {
    return AGM::vector(std::vector<double>{at(1) * src.at(2) - at(2) * src.at(1), at(2) * src.at(0) - at(0) * src.at(2),
                                           at(0) * src.at(1) - at(1) * src.at(0)});
}

AGM::vector AGM::vector::unitVector() {
    AGM::vector vec = vector{};
    double n = norm();
    for (const auto &i: *this) {
        vec.emplace_back(i / n);
    }
    return vec;
}

AGM::vector AGM::vector::operator+(const vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator+",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    AGM::vector vec = vector{};
    for (int i = 0; i < size(); ++i) {
        vec.emplace_back(at(i) + src.at(i));
    }
    return vec;
}

AGM::vector AGM::vector::operator-(const vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator-",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    AGM::vector vec = vector{};
    for (int i = 0; i < size(); ++i) {
        vec.emplace_back(at(i) - src.at(i));
    }
    return vec;
}

AGM::vector AGM::vector::operator*(double d) {
    AGM::vector vec = vector{};
    for (const auto &i: *this) {
        vec.emplace_back(i * d);
    }
    return vec;
}

AGM::vector AGM::vector::operator/(double d) {
    AGM::vector vec = vector{};
    for (const auto &i: *this) {
        vec.emplace_back(i / d);
    }
    return vec;
}

AGM::vector AGM::vector::operator+(const vector &src) const {
    if (size() != src.size())
        printError("AGM::src::operator+",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    AGM::vector vec = vector{};
    for (int i = 0; i < size(); ++i) {
        vec.emplace_back(at(i) + src.at(i));
    }
    return vec;
}

AGM::vector AGM::vector::operator-(const vector &src) const {
    if (size() != src.size())
        printError("AGM::src::operator+",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    AGM::vector vec = vector{};
    for (int i = 0; i < size(); ++i) {
        vec.emplace_back(at(i) - src.at(i));
    }
    return vec;
}

AGM::vector AGM::vector::operator*(double d) const {
    AGM::vector vec = vector{};
    for (const auto &i: *this) {
        vec.emplace_back(i * d);
    }
    return vec;
}

AGM::vector AGM::vector::operator/(double d) const {
    AGM::vector vec = vector{};
    for (const auto &i: *this) {
        vec.emplace_back(i / d);
    }
    return vec;
}

double AGM::vector::operator*(const vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator*",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    double sum{};
    for (int i = 0; i < size(); ++i) {
        sum += at(i) * src.at(i);
    }
    return sum;
}

AGM::vector &AGM::vector::operator+=(const vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator+=",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    for (int i = 0; i < size(); ++i) {
        at(i) += src.at(i);
    }
    return *this;
}

AGM::vector &AGM::vector::operator-=(const vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator-=",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    for (int i = 0; i < size(); ++i) {
        at(i) -= src.at(i);
    }
    return *this;
}

AGM::vector &AGM::vector::operator*=(double d) {
    for (auto &i: *this) {
        i *= d;
    }
    return *this;
}

bool AGM::vector::operator<(const AGM::vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator<",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    for (int i = 0; i < size(); ++i) {
        if (!isclose(at(i), src.at(i))) {
            return at(i) < src.at(i);
        }
    }
    return false;
}

bool AGM::vector::operator>(const AGM::vector &src) {
    if (size() != src.size())
        printError("AGM::src::operator>",
                   "size of this vector (which is %d) is not equal to vector src (which is %d)", size(),
                   src.size());
    for (int i = 0; i < size(); ++i) {
        if (!isclose(at(i), src.at(i))) {
            return at(i) > src.at(i);
        }
    }
    return false;
}
