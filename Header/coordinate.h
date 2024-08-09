//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_COORDINATE_H
#define AGM_COORDINATE_H

#include "matrixRow.h"

namespace AGM {
class coordinate : public std::array<double, 2> {
 public:
  coordinate();

  virtual ~coordinate();

  coordinate(double x, double y);

  [[nodiscard]] auto norm() const -> double;

  auto operator+(const coordinate &src) const -> coordinate;

  auto operator-(const coordinate &src) const -> coordinate;

  auto operator*(double d) const -> coordinate;

  auto operator==(const coordinate &src) const -> bool;

  auto operator!=(const coordinate &rhs) const -> bool;
};
}// namespace AGM

#endif//AGM_COORDINATE_H
