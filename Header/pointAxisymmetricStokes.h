//
// Created by 조준홍 on 2/19/24.
//

#ifndef AGM_POINTAXISYMMETRICSTOKES_H
#define AGM_POINTAXISYMMETRICSTOKES_H

#include "pointStokes.h"

namespace AGM {
class pointAxisymmetricStokes : public pointStokes {
 private:
  bool is_on_axis{false};

 public:
  [[nodiscard]] auto isOnAxis() const -> bool;

  void setIsOnAxis(bool isOnAxis);

  void findStencil(const axialElement *axialElement1, std::vector<pointAxisymmetricStokes> *vector);

  void checkOnAxis();

  void EquationOnAxis();

  void calculateRepresentationFormula(char rz, int order);

  void calculateRepresentationFormulaCross(char rz);

  auto calculateRepresentationFormulaCrossSymmetricElliptic() -> matrixRow;

  auto calculateRepresentationFormulaCrossSymmetricReactionDiffusion() -> matrixRow;

  auto calculateRepresentationFormulaCrossSymmetricNearAxis() -> matrixRow;

  auto calculateRepresentationFormulaCrossNonSymmetricR() -> matrixRow;

  auto calculateRepresentationFormulaCrossNonSymmetricZ() -> matrixRow;

  void calculateRepresentationFormulaNeumann(int order, char rz);

  auto calculateRepresentationFormulaNeumannOnAxial(char rz, char axis, int axisInt) -> matrixRow;

  auto calculateRepresentationFormulaNeumannOnAxialSymmetricElliptic() -> matrixRow;

  auto calculateRepresentationFormulaNeumannOnAxialNonSymmetricElliptic() -> matrixRow;

  auto calculateRepresentationFormulaNeumannOffAxial(char rz, char axis, int axisInt) -> matrixRow;

  auto calculateRepresentationFormulaNeumannOffAxialSymmetricElliptic() -> matrixRow;

  auto calculateRepresentationFormulaNeumannOffAxialNonSymmetricElliptic() -> matrixRow;

  void calculateRepresentationFormulaInterface(char rz);

  auto calculateRepresentationFormulaInterfaceSymmetricElliptic() -> matrixRow;

  auto calculateRepresentationFormulaInterfaceSymmetricReactionDiffusion() -> matrixRow;

  auto calculateRepresentationFormulaInterfaceSymmetricNearAxis() -> matrixRow;

  auto calculateRepresentationFormulaInterfaceNonSymmetricR() -> matrixRow;

  auto calculateRepresentationFormulaInterfaceNonSymmetricZ() -> matrixRow;

  auto calculateRepresentationFormulaStokes(const std::vector<pointAxisymmetricStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp, int comp) -> double;

  void approximateNaNPressure(std::vector<pointAxisymmetricStokes> *points);

  void constantExtensionPressure(std::vector<pointAxisymmetricStokes> *points);

  void makeDerivatives(char rz);

  void makeDerivativesCross(char rz);

  void makeDerivativesCrossSymmetricElliptic();

  void makeDerivativesCrossSymmetricReactionDiffusion();

  void makeDerivativesCrossSymmetricNearAxis();

  void makeDerivativesCrossNonSymmetricR();

  void makeDerivativesCrossNonSymmetricZ();

  void makeDerivativesBoundary(char rz);

  void makeDerivativesInterface(char rz);

  void makeDerivativesInterfaceSymmetricElliptic();

  void makeDerivativesInterfaceSymmetricReactionDiffusion();

  void makeDerivativesInterfaceSymmetricNearAxis();

  void makeDerivativesInterfaceNonSymmetricR();

  void makeDerivativesInterfaceNonSymmetricZ();

  void calculateDerivatives(const std::vector<pointAxisymmetricStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp);

  void approximateNaNDerivatives(std::vector<pointAxisymmetricStokes> *points);

  void makeRepresentationFormulaStokes(char rz);

  void makeRepresentationFormulaStokesCross(char rz);

  void makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion0();

  void makeRepresentationFormulaStokesCrossNonSymmetricZ0();

  void makeRepresentationFormulaStokesBoundary(char rz);

  void makeRepresentationFormulaStokesInterface(char rz);

  void makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion0();

  void makeRepresentationFormulaStokesInterfaceNonSymmetricZ0();

  void makeRepresentationFormulaStokesFull(char rz);

  void makeRepresentationFormulaStokesCrossFull(char rz);

  void makeRepresentationFormulaStokesCrossSymmetricElliptic();

  void makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion();

  void makeRepresentationFormulaStokesCrossSymmetricNearAxis();

  void makeRepresentationFormulaStokesCrossNonSymmetricR();

  void makeRepresentationFormulaStokesCrossNonSymmetricZ();

  void makeRepresentationFormulaStokesInterfaceFull(char rz);

  void makeRepresentationFormulaStokesInterfaceSymmetricElliptic();

  void makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion();

  void makeRepresentationFormulaStokesInterfaceSymmetricNearAxis();

  void makeRepresentationFormulaStokesInterfaceNonSymmetricR();

  void makeRepresentationFormulaStokesInterfaceNonSymmetricZ();
};

}// namespace AGM

#endif//AGM_POINTAXISYMMETRICSTOKES_H
