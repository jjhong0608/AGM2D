//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_MATRIXROW_H
#define AGM_MATRIXROW_H

#include "boundaryLine.h"

namespace AGM {
struct matrixElement {
  int idx{};
  double value{};
};

class matrixRow : public std::vector<matrixElement> {
 public:
  void remove(int i);

  auto operator[](int i) -> double &;

  auto operator+(const matrixRow &src) const -> matrixRow;

  auto operator-(const matrixRow &src) const -> matrixRow;

  auto operator*(double d) const -> matrixRow;

  auto operator+=(const matrixRow &src) -> matrixRow;

  auto operator-=(const matrixRow &src) -> matrixRow;
};
}// namespace AGM

#endif//AGM_MATRIXROW_H
