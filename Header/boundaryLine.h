//
// Created by NIMS-JUNHONG on 2022/07/04.
//

#ifndef AGM2D_BOUNDARYLINE_H
#define AGM2D_BOUNDARYLINE_H

#include "axialLine.h"

namespace AGM {
class denseMatrix : public std::vector<std::vector<double>> {
 private:
  int row{}, col{};

 public:
  denseMatrix(int row, int col);

  virtual ~denseMatrix();

  auto operator+(const denseMatrix &src) -> denseMatrix;

  auto operator-(const denseMatrix &src) -> denseMatrix;

  auto operator*(double d) -> denseMatrix;

  auto operator*(const denseMatrix &src) -> denseMatrix;

  auto operator*(const AGM::vector &src) -> AGM::vector;

  auto operator=(const denseMatrix &src) -> denseMatrix &;
};

class line2D {
 private:
  vector start{}, end{}, normal{}, tangent{};
  double length{};

 public:
  line2D();

  line2D(const vector &start, const vector &anEnd);

  line2D(double s0, double s1, double e0, double e1);

  void calcProperties();

  [[nodiscard]] auto getStart() -> vector &;

  [[nodiscard]] auto getStart() const -> const vector &;

  void setStart(const vector &vector);

  [[nodiscard]] auto getAnEnd() -> vector &;

  [[nodiscard]] auto getAnEnd() const -> const vector &;

  void setAnEnd(const vector &vector);

  [[nodiscard]] auto getNormal() -> vector &;

  [[nodiscard]] auto getNormal() const -> const vector &;

  [[nodiscard]] auto getTangent() -> vector &;

  [[nodiscard]] auto getTangent() const -> const vector &;

  [[nodiscard]] auto getLength() const -> double;

  auto iscross(const line2D &src, AGM::vector &vec) -> bool;

  virtual ~line2D();
};

class boundaryLine2D : public line2D {
 private:
  char condition{};
  double boundary_value{};

 public:
  boundaryLine2D();

  boundaryLine2D(const vector &start, const vector &anEnd, char condition, double boundaryValue);

  boundaryLine2D(double s0, double s1, double e0, double e1, char condition, double boundaryValue);

  [[nodiscard]] auto getCondition() const -> char;

  void setCondition(char i);

  [[nodiscard]] auto getBoundaryValue() const -> double;

  void setBoundaryValue(double boundaryValue);
};
}// namespace AGM

#endif//AGM2D_BOUNDARYLINE_H
