//
// Created by 조준홍 on 10/27/23.
//

#ifndef AGM_POINTSTOKES_H
#define AGM_POINTSTOKES_H

#include "pointHeat.h"

namespace AGM {
class pointStokes : public point {
 protected:
  std::array<matrixRow, 2> pMatrixRow{}, prhsMatrixRow{};
  std::array<double, 2> srb{};

 public:
  virtual void findStencil(const axialElement *axialElement1, std::vector<pointStokes> *vector);

  void makeRepresentationFormulaStokes();

  void makeRepresentationFormulaStokesCross();

  void makeRepresentationFormulaStokesBoundary();

  void makeRepresentationFormulaStokesInterface();

  auto calculateRepresentationFormulaStokesBoundaryAxial(char axis, int axisInt) -> AGM::matrixRow;

  auto calculateRepresentationFormulaStokes(const std::vector<pointStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp, int comp) -> double;

  void calculateDerivatives(const std::vector<pointStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp);

  void approximateNaNDerivatives(std::vector<pointStokes> *points);

  void approximateNaNPressure(std::vector<pointStokes> *points);

  [[nodiscard]] auto getPMatrixRow() const -> const std::array<matrixRow, 2> &;

  void setPMatrixRow(const std::array<matrixRow, 2> &matrix_row);

  [[nodiscard]] auto getPartMatrixRow() const -> const std::array<matrixRow, 2> &;

  void makeRepresentationFormulaStokesFull();

  void makeRepresentationFormulaStokesCrossFull();

  void makeRepresentationFormulaStokesBoundaryFull();

  void makeRepresentationFormulaStokesInterfaceFull();

  void updateRightHandSideFull(const std::function<double(int)> &f, const std::function<double(int)> &g);

  void updateRightHandSideCrossFull(const std::function<double(int)> &f, const std::function<double(int)> &g);

  void updateRightHandSideInterfaceFull(const std::function<double(int)> &f, const std::function<double(int)> &g);

  [[nodiscard]] auto getSrb() const -> const std::array<double, 2> &;

  void setSrb(const std::array<double, 2> &array);
};
}// namespace AGM

#endif//AGM_POINTSTOKES_H
