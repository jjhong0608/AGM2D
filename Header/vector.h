//
// Created by NIMS-JUNHONG on 2022/07/04.
//

#ifndef AGM2D_VECTOR_H
#define AGM2D_VECTOR_H

#include "util.h"

namespace AGM {
class vector : public std::vector<double> {
 public:
  vector();

  explicit vector(const std::vector<double> &x);

  virtual ~vector();

  auto norm() -> double;

  auto dot(const vector &src) -> double;

  auto cross(const vector &src) -> vector;

  auto unitVector() -> vector;

  auto operator+(const vector &src) -> vector;

  auto operator-(const vector &src) -> vector;

  auto operator*(double d) -> vector;

  auto operator/(double d) -> vector;

  auto operator+(const vector &src) const -> vector;

  auto operator-(const vector &src) const -> vector;

  auto operator*(double d) const -> vector;

  auto operator/(double d) const -> vector;

  auto operator*(const vector &src) -> double;

  auto operator+=(const vector &src) -> vector &;

  auto operator-=(const vector &src) -> vector &;

  auto operator*=(double d) -> vector &;

  auto operator<(const vector &src) -> bool;

  auto operator>(const vector &src) -> bool;
};

}// namespace AGM

#endif//AGM2D_VECTOR_H
