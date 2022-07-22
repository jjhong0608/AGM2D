//
// Created by NIMS-JUNHONG on 2022/07/19.
//

#ifndef AGM2D_GREENFUNCTIONAXISYMMETRIC_H
#define AGM2D_GREENFUNCTIONAXISYMMETRIC_H

#include "Greenfunction.h"

namespace AGM {
    class GreenfunctionAxisymmetric : public Greenfunction {
    public:
        GreenfunctionAxisymmetric(double tm, double tau, double tp, double mpl, double mpr);

        [[nodiscard]] static double L(double d) ;

        [[nodiscard]] double green_function(double t) const override;

        [[nodiscard]] double green_function_t(double t) const override;

        [[nodiscard]] double green_function_tau(double t) const override;

        [[nodiscard]] double green_function_ttau(double t) const override;

        [[nodiscard]] double integrate_square(char i) const override;

        [[nodiscard]] double integrate_linear(char i) const override;

        [[nodiscard]] double integrate_const(char i) const override;

        [[nodiscard]] double integrate_square_t(char i) const override;

        [[nodiscard]] double integrate_linear_t(char i) const override;

        [[nodiscard]] double integrate_const_t(char i) const override;

        [[nodiscard]] double integrate_square_tau(char i) const override;

        [[nodiscard]] double integrate_linear_tau(char i) const override;

        [[nodiscard]] double integrate_const_tau(char i) const override;

        [[nodiscard]] double integrate_square_ttau(char i) const override;

        [[nodiscard]] double integrate_linear_ttau(char i) const override;

        [[nodiscard]] double integrate_const_ttau(char i) const override;

        [[nodiscard]] double green_function_ND(double t) const override;

        [[nodiscard]] double green_function_t_ND(double t) const override;

        [[nodiscard]] double green_function_tau_ND(double t) const override;

        [[nodiscard]] double green_function_ttau_ND(double t) const override;

        [[nodiscard]] double integrate_square_ND(char pos) const override;

        [[nodiscard]] double integrate_linear_ND(char pos) const override;

        [[nodiscard]] double integrate_const_ND(char pos) const override;

        [[nodiscard]] double integrate_square_t_ND(char pos) const override;

        [[nodiscard]] double integrate_linear_t_ND(char pos) const override;

        [[nodiscard]] double integrate_const_t_ND(char pos) const override;

        [[nodiscard]] double integrate_square_tau_ND(char pos) const override;

        [[nodiscard]] double integrate_linear_tau_ND(char pos) const override;

        [[nodiscard]] double integrate_const_tau_ND(char pos) const override;

        [[nodiscard]] double integrate_square_ttau_ND(char pos) const override;

        [[nodiscard]] double integrate_linear_ttau_ND(char pos) const override;

        [[nodiscard]] double integrate_const_ttau_ND(char pos) const override;

        [[nodiscard]] double green_function_DN(double t) const override;

        [[nodiscard]] double green_function_t_DN(double t) const override;

        [[nodiscard]] double green_function_tau_DN(double t) const override;

        [[nodiscard]] double green_function_ttau_DN(double t) const override;

        [[nodiscard]] double integrate_square_DN(char pos) const override;

        [[nodiscard]] double integrate_linear_DN(char pos) const override;

        [[nodiscard]] double integrate_const_DN(char pos) const override;

        [[nodiscard]] double integrate_square_t_DN(char pos) const override;

        [[nodiscard]] double integrate_linear_t_DN(char pos) const override;

        [[nodiscard]] double integrate_const_t_DN(char pos) const override;

        [[nodiscard]] double integrate_square_tau_DN(char pos) const override;

        [[nodiscard]] double integrate_linear_tau_DN(char pos) const override;

        [[nodiscard]] double integrate_const_tau_DN(char pos) const override;

        [[nodiscard]] double integrate_square_ttau_DN(char pos) const override;

        [[nodiscard]] double integrate_linear_ttau_DN(char pos) const override;

        [[nodiscard]] double integrate_const_ttau_DN(char pos) const override;

    };


}


#endif //AGM2D_GREENFUNCTIONAXISYMMETRIC_H
