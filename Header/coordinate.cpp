//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include "coordinate.h"

AGM::coordinate::coordinate() : std::array<double, 2>{} {

}

AGM::coordinate::coordinate(double x, double y) : std::array<double, 2>{x, y} {

}

double AGM::coordinate::norm() const {
    return std::sqrt(at(0) * at(0) + at(1) * at(1));
}

AGM::coordinate AGM::coordinate::operator+(const AGM::coordinate &src) const {
    return {at(0) + src.at(0), at(1) + src.at(1)};
}

AGM::coordinate AGM::coordinate::operator-(const AGM::coordinate &src) const {
    return {at(0) - src.at(0), at(1) - src.at(1)};
}

AGM::coordinate AGM::coordinate::operator*(double d) const {
    return {at(0) * d, at(1) * d};
}

bool AGM::coordinate::operator==(const AGM::coordinate &src) const {
    return isclose(at(0), src.at(0)) && isclose(at(1), src.at(1));
}

bool AGM::coordinate::operator!=(const AGM::coordinate &src) const {
    return !(*this == src);
}

AGM::coordinate::~coordinate() = default;