//
// Created by 조준홍 on 12/3/23.
//

#ifndef AGM_GREENFUNCTIONAXISYMMETRICSTOKES_H
#define AGM_GREENFUNCTIONAXISYMMETRICSTOKES_H

#include "GreenfunctionReactionDiffusion.h"

namespace AGM {

class GreenfunctionAxisymmetricStokes : public Greenfunction {
 private:
  double a1{}, a2{};
  double b1{}, b2{};
  double c1{}, c2{};
  double d1{}, d2{};
  double e1{}, e2{};

 public:
  GreenfunctionAxisymmetricStokes(double tm, double tau, double tp, double mpl, double mpr);

  [[nodiscard]] auto green_function(double t) const -> double override;

  [[nodiscard]] auto green_function_t(double t) const -> double override;

  [[nodiscard]] auto green_function_tau(double t) const -> double override;

  [[nodiscard]] auto green_function_ttau(double t) const -> double override;

  [[nodiscard]] auto integrate_square(char i) const -> double override;

  [[nodiscard]] auto integrate_linear(char i) const -> double override;

  [[nodiscard]] auto integrate_const(char i) const -> double override;

  [[nodiscard]] auto integrate_square_t(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_t(char i) const -> double override;

  [[nodiscard]] auto integrate_const_t(char i) const -> double override;

  [[nodiscard]] auto integrate_square_tau(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_tau(char i) const -> double override;

  [[nodiscard]] auto integrate_const_tau(char i) const -> double override;

  [[nodiscard]] auto integrate_square_ttau(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_ttau(char i) const -> double override;

  [[nodiscard]] auto integrate_const_ttau(char i) const -> double override;
};

}// namespace AGM

#endif//AGM_GREENFUNCTIONAXISYMMETRICSTOKES_H
