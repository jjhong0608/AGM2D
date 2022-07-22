//
// Created by NIMS-JUNHONG on 2022/07/19.
//

#include "GreenfunctionAxisymmetric.h"

AGM::GreenfunctionAxisymmetric::GreenfunctionAxisymmetric(double tm, double tau, double tp, double mpl, double mpr)
        : Greenfunction(tm, tau, tp, mpl, mpr) {}


double AGM::GreenfunctionAxisymmetric::L(double d) {
    return std::log(std::abs(d));
}

double AGM::GreenfunctionAxisymmetric::green_function(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (L(t) - L(tm)) * (L(tau) - L(tp)) / (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)));
    } else {
        return (L(t) - L(tp)) * (L(tau) - L(tm)) / (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)));
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_t(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (L(tau) - L(tp)) / (t * (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm))));
    } else {
        return (L(tau) - L(tm)) / (t * (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm))));
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_tau(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (L(t) - L(tm)) * (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)) - (mpl - mpr) * (L(tau) - L(tp))) /
               (tau * pow(mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)), 2));
    } else {
        return (L(t) - L(tp)) * (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)) - (mpl - mpr) * (L(tau) - L(tm))) /
               (tau * pow(mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)), 2));
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_ttau(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)) - (mpl - mpr) * (L(tau) - L(tp))) /
               (t * tau * pow(mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)), 2));
    } else {
        return (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)) - (mpl - mpr) * (L(tau) - L(tm))) /
               (t * tau * pow(mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)), 2));
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 9.0) * (3 * pow(tau, 3) * pow(L(tau), 2) - 3 * pow(tau, 3) * L(tau) * L(tm) -
                                  3 * pow(tau, 3) * L(tau) * L(tp) - pow(tau, 3) * L(tau) +
                                  3 * pow(tau, 3) * L(tm) * L(tp) + pow(tau, 3) * L(tp) + pow(tm, 3) * L(tau) -
                                  pow(tm, 3) * L(tp)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        case 'r':
            return (1.0 / 9.0) * (-3 * pow(tau, 3) * pow(L(tau), 2) + 3 * pow(tau, 3) * L(tau) * L(tm) +
                                  3 * pow(tau, 3) * L(tau) * L(tp) + pow(tau, 3) * L(tau) -
                                  3 * pow(tau, 3) * L(tm) * L(tp) - pow(tau, 3) * L(tm) - pow(tp, 3) * L(tau) +
                                  pow(tp, 3) * L(tm)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 4.0) * (2 * pow(tau, 2) * pow(L(tau), 2) - 2 * pow(tau, 2) * L(tau) * L(tm) -
                                  2 * pow(tau, 2) * L(tau) * L(tp) - pow(tau, 2) * L(tau) +
                                  2 * pow(tau, 2) * L(tm) * L(tp) + pow(tau, 2) * L(tp) + pow(tm, 2) * L(tau) -
                                  pow(tm, 2) * L(tp)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        case 'r':
            return (1.0 / 4.0) * (-2 * pow(tau, 2) * pow(L(tau), 2) + 2 * pow(tau, 2) * L(tau) * L(tm) +
                                  2 * pow(tau, 2) * L(tau) * L(tp) + pow(tau, 2) * L(tau) -
                                  2 * pow(tau, 2) * L(tm) * L(tp) - pow(tau, 2) * L(tm) - pow(tp, 2) * L(tau) +
                                  pow(tp, 2) * L(tm)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const(char i) const {
    switch (i) {
        case 'l':
            return (tau * pow(L(tau), 2) - tau * L(tau) * L(tm) - tau * L(tau) * L(tp) - tau * L(tau) +
                    tau * L(tm) * L(tp) + tau * L(tp) + tm * L(tau) - tm * L(tp)) /
                   (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        case 'r':
            return (-tau * pow(L(tau), 2) + tau * L(tau) * L(tm) + tau * L(tau) * L(tp) + tau * L(tau) -
                    tau * L(tm) * L(tp) - tau * L(tm) - tp * L(tau) + tp * L(tm)) /
                   (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_t(char i) const {
    switch (i) {
        case 'l':
            return ((1.0 / 2.0) * pow(tau, 2) * (L(tau) - L(tp)) + (1.0 / 2.0) * pow(tm, 2) * (-L(tau) + L(tp))) /
                   (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        case 'r':
            return ((1.0 / 2.0) * pow(tau, 2) * (-L(tau) + L(tm)) + (1.0 / 2.0) * pow(tp, 2) * (L(tau) - L(tm))) /
                   (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_t", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_t(char i) const {
    switch (i) {
        case 'l':
            return (tau - tm) * (L(tau) - L(tp)) / (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)));
        case 'r':
            return (-tau + tp) * (L(tau) - L(tm)) / (mpl * (L(tau) - L(tp)) + mpr * (-L(tau) + L(tm)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_t", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_t(char i) const {
    switch (i) {
        case 'l':
            return (L(tau) - L(tm)) * (L(tau) - L(tp)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        case 'r':
            return (-L(tau) + L(tp)) * (L(tau) - L(tm)) / (mpl * L(tau) - mpl * L(tp) - mpr * L(tau) + mpr * L(tm));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_t", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_tau(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 9.0) * mpr * (pow(tau, 3) * (-3 * pow(L(tm), 2) + 3 * L(tm) * L(tp) - L(tm) + L(tp)) +
                                        pow(tm, 3) * (3 * pow(L(tm), 2) - 3 * L(tm) * L(tp) + L(tm) - L(tp)) +
                                        3 * (pow(tau, 3) * L(tau) - pow(tm, 3) * L(tm)) * (L(tm) - L(tp))) / (tau *
                                                                                                              (pow(mpl,
                                                                                                                   2) *
                                                                                                               pow(L(tau),
                                                                                                                   2) -
                                                                                                               2 *
                                                                                                               pow(mpl,
                                                                                                                   2) *
                                                                                                               L(tau) *
                                                                                                               L(tp) +
                                                                                                               pow(mpl,
                                                                                                                   2) *
                                                                                                               pow(L(tp),
                                                                                                                   2) -
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               pow(L(tau),
                                                                                                                   2) +
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tau) *
                                                                                                               L(tm) +
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tau) *
                                                                                                               L(tp) -
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tm) *
                                                                                                               L(tp) +
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               pow(L(tau),
                                                                                                                   2) -
                                                                                                               2 *
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               L(tau) *
                                                                                                               L(tm) +
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               pow(L(tm),
                                                                                                                   2)));
        case 'r':
            return (1.0 / 9.0) * mpl * (pow(tau, 3) * (3 * L(tm) * L(tp) + L(tm) - 3 * pow(L(tp), 2) - L(tp)) +
                                        pow(tp, 3) * (-3 * L(tm) * L(tp) - L(tm) + 3 * pow(L(tp), 2) + L(tp)) +
                                        3 * (-pow(tau, 3) * L(tau) + pow(tp, 3) * L(tp)) * (L(tm) - L(tp))) / (tau *
                                                                                                               (pow(mpl,
                                                                                                                    2) *
                                                                                                                pow(L(tau),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                pow(mpl,
                                                                                                                    2) *
                                                                                                                L(tau) *
                                                                                                                L(tp) +
                                                                                                                pow(mpl,
                                                                                                                    2) *
                                                                                                                pow(L(tp),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                pow(L(tau),
                                                                                                                    2) +
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tau) *
                                                                                                                L(tm) +
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tau) *
                                                                                                                L(tp) -
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tm) *
                                                                                                                L(tp) +
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                pow(L(tau),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                L(tau) *
                                                                                                                L(tm) +
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                pow(L(tm),
                                                                                                                    2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_tau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_tau(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 4.0) * mpr * (pow(tau, 2) * (-2 * pow(L(tm), 2) + 2 * L(tm) * L(tp) - L(tm) + L(tp)) +
                                        pow(tm, 2) * (2 * pow(L(tm), 2) - 2 * L(tm) * L(tp) + L(tm) - L(tp)) +
                                        2 * (pow(tau, 2) * L(tau) - pow(tm, 2) * L(tm)) * (L(tm) - L(tp))) / (tau *
                                                                                                              (pow(mpl,
                                                                                                                   2) *
                                                                                                               pow(L(tau),
                                                                                                                   2) -
                                                                                                               2 *
                                                                                                               pow(mpl,
                                                                                                                   2) *
                                                                                                               L(tau) *
                                                                                                               L(tp) +
                                                                                                               pow(mpl,
                                                                                                                   2) *
                                                                                                               pow(L(tp),
                                                                                                                   2) -
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               pow(L(tau),
                                                                                                                   2) +
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tau) *
                                                                                                               L(tm) +
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tau) *
                                                                                                               L(tp) -
                                                                                                               2 * mpl *
                                                                                                               mpr *
                                                                                                               L(tm) *
                                                                                                               L(tp) +
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               pow(L(tau),
                                                                                                                   2) -
                                                                                                               2 *
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               L(tau) *
                                                                                                               L(tm) +
                                                                                                               pow(mpr,
                                                                                                                   2) *
                                                                                                               pow(L(tm),
                                                                                                                   2)));
        case 'r':
            return (1.0 / 4.0) * mpl * (pow(tau, 2) * (2 * L(tm) * L(tp) + L(tm) - 2 * pow(L(tp), 2) - L(tp)) +
                                        pow(tp, 2) * (-2 * L(tm) * L(tp) - L(tm) + 2 * pow(L(tp), 2) + L(tp)) +
                                        2 * (-pow(tau, 2) * L(tau) + pow(tp, 2) * L(tp)) * (L(tm) - L(tp))) / (tau *
                                                                                                               (pow(mpl,
                                                                                                                    2) *
                                                                                                                pow(L(tau),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                pow(mpl,
                                                                                                                    2) *
                                                                                                                L(tau) *
                                                                                                                L(tp) +
                                                                                                                pow(mpl,
                                                                                                                    2) *
                                                                                                                pow(L(tp),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                pow(L(tau),
                                                                                                                    2) +
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tau) *
                                                                                                                L(tm) +
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tau) *
                                                                                                                L(tp) -
                                                                                                                2 *
                                                                                                                mpl *
                                                                                                                mpr *
                                                                                                                L(tm) *
                                                                                                                L(tp) +
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                pow(L(tau),
                                                                                                                    2) -
                                                                                                                2 *
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                L(tau) *
                                                                                                                L(tm) +
                                                                                                                pow(mpr,
                                                                                                                    2) *
                                                                                                                pow(L(tm),
                                                                                                                    2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_tau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_tau(char i) const {
    switch (i) {
        case 'l':
            return mpr * (tau * L(tau) * L(tm) - tau * L(tau) * L(tp) - tau * pow(L(tm), 2) + tau * L(tm) * L(tp) -
                          tau * L(tm) + tau * L(tp) + tm * L(tm) - tm * L(tp)) / (tau * (pow(mpl, 2) * pow(L(tau), 2) -
                                                                                         2 * pow(mpl, 2) * L(tau) *
                                                                                         L(tp) +
                                                                                         pow(mpl, 2) * pow(L(tp), 2) -
                                                                                         2 * mpl * mpr *
                                                                                         pow(L(tau), 2) +
                                                                                         2 * mpl * mpr * L(tau) *
                                                                                         L(tm) +
                                                                                         2 * mpl * mpr * L(tau) *
                                                                                         L(tp) -
                                                                                         2 * mpl * mpr * L(tm) * L(tp) +
                                                                                         pow(mpr, 2) * pow(L(tau), 2) -
                                                                                         2 * pow(mpr, 2) * L(tau) *
                                                                                         L(tm) +
                                                                                         pow(mpr, 2) * pow(L(tm), 2)));
        case 'r':
            return mpl * (-tau * L(tau) * L(tm) + tau * L(tau) * L(tp) + tau * L(tm) * L(tp) + tau * L(tm) -
                          tau * pow(L(tp), 2) - tau * L(tp) - tp * L(tm) + tp * L(tp)) / (tau * (pow(mpl, 2) *
                                                                                                 pow(L(tau), 2) -
                                                                                                 2 * pow(mpl, 2) *
                                                                                                 L(tau) * L(tp) +
                                                                                                 pow(mpl, 2) *
                                                                                                 pow(L(tp), 2) -
                                                                                                 2 * mpl * mpr *
                                                                                                 pow(L(tau), 2) +
                                                                                                 2 * mpl * mpr *
                                                                                                 L(tau) * L(tm) +
                                                                                                 2 * mpl * mpr *
                                                                                                 L(tau) * L(tp) -
                                                                                                 2 * mpl * mpr * L(tm) *
                                                                                                 L(tp) + pow(mpr, 2) *
                                                                                                         pow(L(tau),
                                                                                                             2) -
                                                                                                 2 * pow(mpr, 2) *
                                                                                                 L(tau) * L(tm) +
                                                                                                 pow(mpr, 2) *
                                                                                                 pow(L(tm), 2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_tau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_ttau(char i) const {
    switch (i) {
        case 'l':
            return (1.0 / 2.0) * mpr * (pow(tau, 2) * (L(tm) - L(tp)) + pow(tm, 2) * (-L(tm) + L(tp))) / (tau *
                                                                                                          (pow(mpl, 2) *
                                                                                                           pow(L(tau),
                                                                                                               2) - 2 *
                                                                                                                    pow(mpl,
                                                                                                                        2) *
                                                                                                                    L(tau) *
                                                                                                                    L(tp) +
                                                                                                           pow(mpl, 2) *
                                                                                                           pow(L(tp),
                                                                                                               2) -
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           pow(L(tau),
                                                                                                               2) +
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           L(tau) *
                                                                                                           L(tm) +
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           L(tau) *
                                                                                                           L(tp) -
                                                                                                           2 * mpl *
                                                                                                           mpr * L(tm) *
                                                                                                           L(tp) +
                                                                                                           pow(mpr, 2) *
                                                                                                           pow(L(tau),
                                                                                                               2) - 2 *
                                                                                                                    pow(mpr,
                                                                                                                        2) *
                                                                                                                    L(tau) *
                                                                                                                    L(tm) +
                                                                                                           pow(mpr, 2) *
                                                                                                           pow(L(tm),
                                                                                                               2)));
        case 'r':
            return (1.0 / 2.0) * mpl * (pow(tau, 2) * (-L(tm) + L(tp)) + pow(tp, 2) * (L(tm) - L(tp))) / (tau *
                                                                                                          (pow(mpl, 2) *
                                                                                                           pow(L(tau),
                                                                                                               2) - 2 *
                                                                                                                    pow(mpl,
                                                                                                                        2) *
                                                                                                                    L(tau) *
                                                                                                                    L(tp) +
                                                                                                           pow(mpl, 2) *
                                                                                                           pow(L(tp),
                                                                                                               2) -
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           pow(L(tau),
                                                                                                               2) +
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           L(tau) *
                                                                                                           L(tm) +
                                                                                                           2 * mpl *
                                                                                                           mpr *
                                                                                                           L(tau) *
                                                                                                           L(tp) -
                                                                                                           2 * mpl *
                                                                                                           mpr * L(tm) *
                                                                                                           L(tp) +
                                                                                                           pow(mpr, 2) *
                                                                                                           pow(L(tau),
                                                                                                               2) - 2 *
                                                                                                                    pow(mpr,
                                                                                                                        2) *
                                                                                                                    L(tau) *
                                                                                                                    L(tm) +
                                                                                                           pow(mpr, 2) *
                                                                                                           pow(L(tm),
                                                                                                               2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_ttau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_ttau(char i) const {
    switch (i) {
        case 'l':
            return mpr * (tau * L(tm) - tau * L(tp) - tm * L(tm) + tm * L(tp)) / (tau * (pow(mpl, 2) * pow(L(tau), 2) -
                                                                                         2 * pow(mpl, 2) * L(tau) *
                                                                                         L(tp) +
                                                                                         pow(mpl, 2) * pow(L(tp), 2) -
                                                                                         2 * mpl * mpr *
                                                                                         pow(L(tau), 2) +
                                                                                         2 * mpl * mpr * L(tau) *
                                                                                         L(tm) +
                                                                                         2 * mpl * mpr * L(tau) *
                                                                                         L(tp) -
                                                                                         2 * mpl * mpr * L(tm) * L(tp) +
                                                                                         pow(mpr, 2) * pow(L(tau), 2) -
                                                                                         2 * pow(mpr, 2) * L(tau) *
                                                                                         L(tm) +
                                                                                         pow(mpr, 2) * pow(L(tm), 2)));
        case 'r':
            return mpl * (-tau * L(tm) + tau * L(tp) + tp * L(tm) - tp * L(tp)) / (tau * (pow(mpl, 2) * pow(L(tau), 2) -
                                                                                          2 * pow(mpl, 2) * L(tau) *
                                                                                          L(tp) +
                                                                                          pow(mpl, 2) * pow(L(tp), 2) -
                                                                                          2 * mpl * mpr *
                                                                                          pow(L(tau), 2) +
                                                                                          2 * mpl * mpr * L(tau) *
                                                                                          L(tm) +
                                                                                          2 * mpl * mpr * L(tau) *
                                                                                          L(tp) -
                                                                                          2 * mpl * mpr * L(tm) *
                                                                                          L(tp) +
                                                                                          pow(mpr, 2) * pow(L(tau), 2) -
                                                                                          2 * pow(mpr, 2) * L(tau) *
                                                                                          L(tm) +
                                                                                          pow(mpr, 2) * pow(L(tm), 2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_ttau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_ttau(char i) const {
    switch (i) {
        case 'l':
            return mpr * (L(tau) * L(tm) - L(tau) * L(tp) - pow(L(tm), 2) + L(tm) * L(tp)) / (tau * (pow(mpl, 2) *
                                                                                                     pow(L(tau), 2) -
                                                                                                     2 * pow(mpl, 2) *
                                                                                                     L(tau) * L(tp) +
                                                                                                     pow(mpl, 2) *
                                                                                                     pow(L(tp), 2) -
                                                                                                     2 * mpl * mpr *
                                                                                                     pow(L(tau), 2) +
                                                                                                     2 * mpl * mpr *
                                                                                                     L(tau) * L(tm) +
                                                                                                     2 * mpl * mpr *
                                                                                                     L(tau) * L(tp) -
                                                                                                     2 * mpl * mpr *
                                                                                                     L(tm) * L(tp) +
                                                                                                     pow(mpr, 2) *
                                                                                                     pow(L(tau), 2) -
                                                                                                     2 * pow(mpr, 2) *
                                                                                                     L(tau) * L(tm) +
                                                                                                     pow(mpr, 2) *
                                                                                                     pow(L(tm), 2)));
        case 'r':
            return mpl * (-L(tau) * L(tm) + L(tau) * L(tp) + L(tm) * L(tp) - pow(L(tp), 2)) / (tau * (pow(mpl, 2) *
                                                                                                      pow(L(tau), 2) -
                                                                                                      2 * pow(mpl, 2) *
                                                                                                      L(tau) * L(tp) +
                                                                                                      pow(mpl, 2) *
                                                                                                      pow(L(tp), 2) -
                                                                                                      2 * mpl * mpr *
                                                                                                      pow(L(tau), 2) +
                                                                                                      2 * mpl * mpr *
                                                                                                      L(tau) * L(tm) +
                                                                                                      2 * mpl * mpr *
                                                                                                      L(tau) * L(tp) -
                                                                                                      2 * mpl * mpr *
                                                                                                      L(tm) * L(tp) +
                                                                                                      pow(mpr, 2) *
                                                                                                      pow(L(tau), 2) -
                                                                                                      2 * pow(mpr, 2) *
                                                                                                      L(tau) * L(tm) +
                                                                                                      pow(mpr, 2) *
                                                                                                      pow(L(tm), 2)));
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_ttau", "char i (which is %c) wrong", i);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (-L(tau) + L(tp)) / mpr;
    } else {
        return (-L(t) + L(tp)) / mpr;
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_t_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return -UNITVALUE / (mpr * t);
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_tau_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return -UNITVALUE / (mpr * tau);
    } else {
        return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::green_function_ttau_ND(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 3.0) * (pow(tau, 3) * (-L(tau) + L(tp)) + pow(tm, 3) * (L(tau) - L(tp))) / mpr;
        case 'r':
            return (1.0 / 9.0) * (3 * pow(tau, 3) * L(tau) - 3 * pow(tau, 3) * L(tp) - pow(tau, 3) + pow(tp, 3)) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 2.0) * (pow(tau, 2) * (-L(tau) + L(tp)) + pow(tm, 2) * (L(tau) - L(tp))) / mpr;
        case 'r':
            return (1.0 / 4.0) * (2 * pow(tau, 2) * L(tau) - 2 * pow(tau, 2) * L(tp) - pow(tau, 2) + pow(tp, 2)) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (-tau + tm) * (L(tau) - L(tp)) / mpr;
        case 'r':
            return (tau * L(tau) - tau * L(tp) - tau + tp) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_t_ND(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 2.0) * (pow(tau, 2) - pow(tp, 2)) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_t_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_t_ND(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (tau - tp) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_t_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_t_ND(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (L(tau) - L(tp)) / mpr;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_t_ND", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_tau_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 3.0) * (-pow(tau, 3) + pow(tm, 3)) / (mpr * tau);
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_tau_ND", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_tau_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 2.0) * (-pow(tau, 2) + pow(tm, 2)) / (mpr * tau);
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_tau_ND", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_tau_ND(char pos) const {
    switch (pos) {
        case 'l':
            return (-tau + tm) / (mpr * tau);
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_tau_ND", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_ttau_ND(char pos) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_ttau_ND(char pos) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionAxisymmetric::integrate_const_ttau_ND(char pos) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionAxisymmetric::green_function_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return (L(t) - L(tm)) / mpl;
    } else {
        return (L(tau) - L(tm)) / mpl;
    }

}

double AGM::GreenfunctionAxisymmetric::green_function_t_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return 1 / (mpl * t);
    } else {
        return ZEROVALUE;
    }

}

double AGM::GreenfunctionAxisymmetric::green_function_tau_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return 1 / (mpl * tau);
    }

}

double AGM::GreenfunctionAxisymmetric::green_function_ttau_DN(double t) const {
    if (t < tau || isclose(tm, t)) {
        return ZEROVALUE;
    } else {
        return ZEROVALUE;
    }

}

double AGM::GreenfunctionAxisymmetric::integrate_square_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 9.0) * (3 * pow(tau, 3) * L(tau) - 3 * pow(tau, 3) * L(tm) - pow(tau, 3) + pow(tm, 3)) / mpl;
        case 'r':
            return (1.0 / 3.0) * (pow(tau, 3) * (-L(tau) + L(tm)) + pow(tp, 3) * (L(tau) - L(tm))) / mpl;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 4.0) * (2 * pow(tau, 2) * L(tau) - 2 * pow(tau, 2) * L(tm) - pow(tau, 2) + pow(tm, 2)) / mpl;
        case 'r':
            return (1.0 / 2.0) * (pow(tau, 2) * (-L(tau) + L(tm)) + pow(tp, 2) * (L(tau) - L(tm))) / mpl;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (tau * L(tau) - tau * L(tm) - tau + tm) / mpl;
        case 'r':
            return (-tau + tp) * (L(tau) - L(tm)) / mpl;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_t_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (1.0 / 2.0) * (pow(tau, 2) - pow(tm, 2)) / mpl;
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_t_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_t_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (tau - tm) / mpl;
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_t_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_t_DN(char pos) const {
    switch (pos) {
        case 'l':
            return (L(tau) - L(tm)) / mpl;
        case 'r':
            return ZEROVALUE;
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_t_DN", "char pos (which is %c) is wrong", pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_tau_DN(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 3.0) * (-pow(tau, 3) + pow(tp, 3)) / (mpl * tau);
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_square_tau_DN", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_tau_DN(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (1.0 / 2.0) * (-pow(tau, 2) + pow(tp, 2)) / (mpl * tau);
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_linear_tau_DN", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_const_tau_DN(char pos) const {
    switch (pos) {
        case 'l':
            return ZEROVALUE;
        case 'r':
            return (-tau + tp) / (mpl * tau);
        default:
            printError("AGM::GreenfunctionAxisymmetric::integrate_const_tau_DN", "char pos (which is %c) is wrong",
                       pos);
            return ZEROVALUE;
    }
}

double AGM::GreenfunctionAxisymmetric::integrate_square_ttau_DN(char pos) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionAxisymmetric::integrate_linear_ttau_DN(char pos) const {
    return ZEROVALUE;
}

double AGM::GreenfunctionAxisymmetric::integrate_const_ttau_DN(char pos) const {
    return ZEROVALUE;
}
