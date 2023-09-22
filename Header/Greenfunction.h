//
// Created by NIMS-JUNHONG on 2022/01/21.
//


#ifndef AGM_GREENFUNCTION_H
#define AGM_GREENFUNCTION_H

#include "vector.h"

namespace AGM {
    class Greenfunction {
    protected:
        double tm{}, tau{}, tp{}, mpl{}, mpr{};

    public:
        Greenfunction(double tm, double tau, double tp, double mpl, double mpr);

        virtual ~Greenfunction();

        [[nodiscard]] virtual auto integrate_square(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_t(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_t(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ttau(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ttau(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_ttau(char pos) const -> double;

        [[nodiscard]] virtual auto green_function(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau(double t) const -> double;

        [[nodiscard]] auto green_integral(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_t(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_tau(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_ttau(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_t_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_tau_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_square(char pos) const -> double;

        [[nodiscard]] auto green_integral_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_t_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_tau_linear(char pos) const -> double;

        [[nodiscard]] auto green_integral_ttau_linear(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_t_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_t_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ttau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ttau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_ttau_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_function_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau_ND(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau_ND(double t) const -> double;

        [[nodiscard]] auto green_integral_ND(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_t_ND(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_tau_ND(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_ttau_ND(char pos, int order = 2) const -> double;

        [[nodiscard]] virtual auto green_integral_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_square_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_linear_ND(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_t_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_t_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ttau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ttau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_ttau_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_function_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau_DN(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau_DN(double t) const -> double;

        [[nodiscard]] auto green_integral_DN(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_t_DN(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_tau_DN(char pos, int order = 2) const -> double;

        [[nodiscard]] auto green_integral_ttau_DN(char pos, int order = 2) const -> double;

        [[nodiscard]] virtual auto green_integral_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_square_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_linear_DN(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_t_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_t_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_t_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_tau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_tau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_tau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_square_ttau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_linear_ttau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto integrate_const_ttau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_function_DF(double t) const -> double;

        [[nodiscard]] virtual auto green_function_t_DF(double t) const -> double;

        [[nodiscard]] virtual auto green_function_tau_DF(double t) const -> double;

        [[nodiscard]] virtual auto green_function_ttau_DF(double t) const -> double;

        [[nodiscard]] virtual auto green_integral_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_square_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_square_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_square_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_square_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_linear_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_t_linear_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_tau_linear_DF(char pos) const -> double;

        [[nodiscard]] virtual auto green_integral_ttau_linear_DF(char pos) const -> double;

    };
}

#endif //AGM_GREENFUNCTION_H
