//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_GREENFUNCTIONREACTIONDIFFUSION_H
#define AGM_GREENFUNCTIONREACTIONDIFFUSION_H

#include "GreenfunctionAxisymmetric.h"

namespace AGM {
class GreenfunctionReactionDiffusion : public Greenfunction {
 protected:
  double c{};
  double alpha{};

 public:
  GreenfunctionReactionDiffusion(double tm, double tau, double tp, double mpl, double mpr, double c);

  [[nodiscard]] auto E(double a, double b) const -> double;

  [[nodiscard]] auto F(double a, double b) const -> double;

  [[nodiscard]] auto l2p(double s) const -> double;

  [[nodiscard]] auto l2m(double s) const -> double;

  [[nodiscard]] auto l1p(double s) const -> double;

  [[nodiscard]] auto l1m(double s) const -> double;

  [[nodiscard]] auto l0(double a, double b) const -> double;

  [[nodiscard]] auto green_function(double t) const -> double override;

  [[nodiscard]] auto green_function_t(double t) const -> double override;

  [[nodiscard]] auto green_function_tau(double t) const -> double override;

  [[nodiscard]] auto green_function_ttau(double t) const -> double override;

  [[nodiscard]] auto green_function_ND(double t) const -> double override;

  [[nodiscard]] auto green_function_t_ND(double t) const -> double override;

  [[nodiscard]] auto green_function_tau_ND(double t) const -> double override;

  [[nodiscard]] auto green_function_ttau_ND(double t) const -> double override;

  [[nodiscard]] auto green_function_DN(double t) const -> double override;

  [[nodiscard]] auto green_function_t_DN(double t) const -> double override;

  [[nodiscard]] auto green_function_tau_DN(double t) const -> double override;

  [[nodiscard]] auto green_function_ttau_DN(double t) const -> double override;

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

  [[nodiscard]] auto integrate_square_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_const_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_square_t_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_t_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_const_t_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_square_tau_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_tau_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_const_tau_ND(char i) const -> double override;

  [[nodiscard]] auto integrate_square_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_const_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_square_t_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_t_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_const_t_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_square_tau_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_linear_tau_DN(char i) const -> double override;

  [[nodiscard]] auto integrate_const_tau_DN(char i) const -> double override;
};
}// namespace AGM

#endif//AGM_GREENFUNCTIONREACTIONDIFFUSION_H
