//
// Created by 조준홍 on 2022/02/13.
//

#ifndef AGM_POINTHEAT_H
#define AGM_POINTHEAT_H

#include "pointAxisymmetric.h"

namespace AGM {
class pointHeat : public point {
 protected:
  static double time, delta;

 public:
  static auto getTime() -> double;

  static void setTime(double d);

  static auto getDelta() -> double;

  static void setDelta(double d);

  void findStencil(const axialElement *axialElement1, std::vector<pointHeat> *vector);

  void calculateRepresentationFormulaCross() override;

  auto calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> matrixRow override;

  auto calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> matrixRow override;

  void calculateRepresentationFormulaInterface() override;

  void makePhiCoefficient(std::vector<pointHeat> *vector);

  void updateRightHandSidePhiPressure(
      const std::function<double(int)> &f,
      const std::function<double(int)> &g,
      std::vector<pointHeat> *points);

  void updateRightHandSidePhiPressureCross(
      const std::function<double(int)> &f,
      const std::function<double(int)> &g,
      std::vector<pointHeat> *points);

  void updateRightHandSidePhiPressureDirichlet(
      const std::function<double(int)> &f,
      const std::function<double(int)> &g,
      std::vector<pointHeat> *points);

  void updateRightHandSidePhiPressureInterface(
      const std::function<double(int)> &f,
      const std::function<double(int)> &g,
      std::vector<pointHeat> *points);

  void makeDerivativesCross() override;

  void calculateDerivatives(const std::vector<pointHeat> *points, const std::function<double(int)> &f,
                            const std::function<double(int)> &g, const std::function<double(int)> &fp,
                            const std::function<double(int)> &gp);

  void approximateNaNDerivatives(std::vector<pointHeat> *points);

  void approximateDiff(std::vector<pointHeat> *points);

  void makeDerivativesInterface() override;
};

}// namespace AGM

#endif//AGM_POINTHEAT_H
