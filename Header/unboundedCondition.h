//
// Created by NIMS-JUNHONG on 2022/08/12.
//

#ifndef AGM2D_UNBOUNDEDCONDITION_H
#define AGM2D_UNBOUNDEDCONDITION_H

#include "GreenfunctionAxisymmetricStokes.h"

namespace AGM {
struct unbounded_f_param {
  double tm{}, tau{}, tp{}, tau0{};
  std::function<double(double, double)> *function{};
};

auto unwrap(double x, void *p) -> double;

class asymptoticBehavior {
 private:
  std::function<double(double)> u0{}, u1{}, phi0{}, phi1{}, f0{}, f1{};

 public:
  asymptoticBehavior();

  asymptoticBehavior(const std::function<double(double)> &u0,
                     const std::function<double(double)> &u1,
                     const std::function<double(double)> &phi0,
                     const std::function<double(double)> &phi1,
                     const std::function<double(double)> &f0,
                     const std::function<double(double)> &f1);

  auto getU0() const -> const std::function<double(double)> &;

  void setU0(const std::function<double(double)> &function);

  auto getU1() const -> const std::function<double(double)> &;

  void setU1(const std::function<double(double)> &function);

  auto getPhi0() const -> const std::function<double(double)> &;

  void setPhi0(const std::function<double(double)> &function);

  auto getPhi1() const -> const std::function<double(double)> &;

  void setPhi1(const std::function<double(double)> &function);

  auto getF0() const -> const std::function<double(double)> &;

  void setF0(const std::function<double(double)> &function);

  auto getF1() const -> const std::function<double(double)> &;

  void setF1(const std::function<double(double)> &function);

  virtual ~asymptoticBehavior();
};

class unboundedCondition {
 private:
  double tm{}, tau{}, tp{}, tau0{}, gamma{6e0};
  asymptoticBehavior ab{};

 public:
  unboundedCondition(double tm, double tau, double tp, double tau0, const asymptoticBehavior &ab);

  auto getAb() const -> const asymptoticBehavior &;

  void setAb(const asymptoticBehavior &asymptoticBehavior);

  auto rho_p(double t) const -> double;

  auto rho_m(double t) const -> double;

  auto A_DF(const std::function<double(double, double)> &f) -> double;

  auto A_FD(const std::function<double(double, double)> &f) -> double;

  virtual ~unboundedCondition();
};

}// namespace AGM

#endif//AGM2D_UNBOUNDEDCONDITION_H
