//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include "coordinate.h"

AGM::coordinate::coordinate() : std::array<double, 2>{} {
}

AGM::coordinate::coordinate(double x, double y) : std::array<double, 2>{x, y} {
}

auto AGM::coordinate::norm() const -> double {
  return std::sqrt(at(0) * at(0) + at(1) * at(1));
}

auto AGM::coordinate::operator+(const AGM::coordinate &src) const -> AGM::coordinate {
  return {at(0) + src.at(0), at(1) + src.at(1)};
}

auto AGM::coordinate::operator-(const AGM::coordinate &src) const -> AGM::coordinate {
  return {at(0) - src.at(0), at(1) - src.at(1)};
}

auto AGM::coordinate::operator*(double d) const -> AGM::coordinate {
  return {at(0) * d, at(1) * d};
}

auto AGM::coordinate::operator==(const AGM::coordinate &src) const -> bool {
  return isclose(at(0), src.at(0)) && isclose(at(1), src.at(1));
}

auto AGM::coordinate::operator!=(const AGM::coordinate &src) const -> bool {
  return !(*this == src);
}

AGM::coordinate::~coordinate() = default;