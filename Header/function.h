//
// Created by 조준홍 on 2022/02/13.
//

#ifndef AGM_FUNCTION_H
#define AGM_FUNCTION_H

#include "structure.h"

namespace AGM {
class ellipticFunction {
 public:
  ellipticFunction();

  virtual ~ellipticFunction();

  auto u(const point &pt) -> double;

  auto phi(const point &pt) -> double;

  auto f(const point &pt) -> double;

  auto ux(const point &pt) -> double;

  auto uy(const point &pt) -> double;

  auto isAssignBoundaryValue() -> bool;

  void assignBoundaryValue(point &pt);

  auto resultPath() -> std::string;

  auto FlowDataPath() -> std::string;
};

class heatFunction {
 public:
  heatFunction();

  virtual ~heatFunction();

  static auto initialTime() -> double;

  static auto terminalTime() -> double;

  static auto deltaTime() -> double;

  auto u(double t, const point &pt) -> double;

  auto phi(double t, const point &pt) -> double;

  auto f(double t, const point &pt) -> double;

  auto ux(double t, const point &pt) -> double;

  auto uy(double t, const point &pt) -> double;

  auto isLoadPreviousValue() -> bool;

  auto isAssignBoundaryValue() -> bool;

  void assignPreviousValue(AGM::value &value, point &pt);

  void assignBoundaryValue(point &pt);

  auto resultPath() -> std::string;
};

class NavierStokesFunction {
 public:
  NavierStokesFunction();

  virtual ~NavierStokesFunction();

  static auto initialTime() -> double;

  static auto terminalTime() -> double;

  static auto deltaTime() -> double;

  static auto writeTime() -> double;

  auto u(double t, const point &pt) -> double;

  auto v(double t, const point &pt) -> double;

  auto p(double t, const point &pt) -> double;

  auto phi(double t, const point &pt) -> double;

  auto psi(double t, const point &pt) -> double;

  auto ux(double t, const point &pt) -> double;

  auto uy(double t, const point &pt) -> double;

  auto vx(double t, const point &pt) -> double;

  auto vy(double t, const point &pt) -> double;

  auto px(double t, const point &pt) -> double;

  auto py(double t, const point &pt) -> double;

  auto f1(double t, const point &pt) -> double;

  auto f2(double t, const point &pt) -> double;

  void assignPreviousValue(AGM::value &pu, AGM::value &pv, AGM::value &pp, point &uvel, point &vvel, point &pres);

  void assignBoundaryValue(AGM::point &uvel, AGM::point &vvel, int presentIter);

  auto isLoadPreviousFile() -> bool;

  auto isAssignBoundaryValue() -> bool;

  void loadPreviousValue(std::vector<value> *pu, std::vector<value> *pv, std::vector<value> *pp);

  auto isNormalEq() -> bool;

  auto findFixedPointIndex(std::vector<AGM::point> *pts) -> int;

  auto tolerance() -> double;

  auto resultPath() -> std::string;
};
}// namespace AGM

#endif//AGM_FUNCTION_H
