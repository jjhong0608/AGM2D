//
// Created by NIMS-JUNHONG on 2022/08/12.
//

#include "unboundedCondition.h"
#include "gsl/gsl_integration.h"
#include <gsl/gsl_errno.h>

auto AGM::unwrap(double x, void *p) -> double {
    auto params{(struct unbounded_f_param *) p};
    auto fp{static_cast<std::function<double(double, double)> *>(params->function)};
    return (*fp)(x, params->tau0);
}

AGM::asymptoticBehavior::asymptoticBehavior() = default;

AGM::asymptoticBehavior::asymptoticBehavior(const std::function<double(double)> &u0,
                                            const std::function<double(double)> &u1,
                                            const std::function<double(double)> &phi0,
                                            const std::function<double(double)> &phi1,
                                            const std::function<double(double)> &f0,
                                            const std::function<double(double)> &f1) {

}

auto AGM::asymptoticBehavior::getU0() const -> const std::function<double(double)> & {
    return u0;
}

void AGM::asymptoticBehavior::setU0(const std::function<double(double)> &function) {
    asymptoticBehavior::u0 = function;
}

auto AGM::asymptoticBehavior::getU1() const -> const std::function<double(double)> & {
    return u1;
}

void AGM::asymptoticBehavior::setU1(const std::function<double(double)> &function) {
    asymptoticBehavior::u1 = function;
}

auto AGM::asymptoticBehavior::getPhi0() const -> const std::function<double(double)> & {
    return phi0;
}

void AGM::asymptoticBehavior::setPhi0(const std::function<double(double)> &function) {
    asymptoticBehavior::phi0 = function;
}

auto AGM::asymptoticBehavior::getPhi1() const -> const std::function<double(double)> & {
    return phi1;
}

void AGM::asymptoticBehavior::setPhi1(const std::function<double(double)> &function) {
    asymptoticBehavior::phi1 = function;
}

auto AGM::asymptoticBehavior::getF0() const -> const std::function<double(double)> & {
    return f0;
}

void AGM::asymptoticBehavior::setF0(const std::function<double(double)> &function) {
    asymptoticBehavior::f0 = function;
}

auto AGM::asymptoticBehavior::getF1() const -> const std::function<double(double)> & {
    return f1;
}

void AGM::asymptoticBehavior::setF1(const std::function<double(double)> &function) {
    asymptoticBehavior::f1 = function;
}

AGM::asymptoticBehavior::~asymptoticBehavior() = default;

AGM::unboundedCondition::unboundedCondition(double tm, double tau, double tp, double tau0,
                                            const AGM::asymptoticBehavior &ab) : tm(tm), tau(tau), tp(tp), tau0(tau0),
                                                                                 ab(ab) {}

auto AGM::unboundedCondition::getAb() const -> const AGM::asymptoticBehavior & {
    return ab;
}

void AGM::unboundedCondition::setAb(const AGM::asymptoticBehavior &asymptoticBehavior) {
    unboundedCondition::ab = asymptoticBehavior;
}

auto AGM::unboundedCondition::rho_p(double t) const -> double {
    return std::pow(t - tm + UNITVALUE, -gamma);
}

auto AGM::unboundedCondition::rho_m(double t) const -> double {
    return std::pow(tp + UNITVALUE - t, -gamma);
}

auto AGM::unboundedCondition::A_DF(const std::function<double(double, double)> &f) -> double {
    auto f0 = std::function<double(double, double)>{[this](double x, double y) -> double {
        return rho_p(x);
    }};
    auto f1 = std::function<double(double, double)>{[this, &f](double x, double y) -> double {
        return f(x, y) * rho_p(x);
    }};
    auto params0{unbounded_f_param{tm, tau, tp, tau0, &f0}};
    auto params1{unbounded_f_param{tm, tau, tp, tau0, &f1}};
    auto w0{gsl_integration_workspace_alloc(1000)};
    auto w1{gsl_integration_workspace_alloc(1000)};
    double result0{}, result1{}, error0{}, error1{};

    auto F0{gsl_function{unwrap, &params0}};
    auto F1{gsl_function{unwrap, &params1}};

    int status0{gsl_integration_qagiu(&F0, tm, 0, 1e-10, 1000, w0, &result0, &error0)};
    int status1{gsl_integration_qagiu(&F1, tm, 0, 1e-10, 1000, w1, &result1, &error1)};

    if (status0) {
        printf("error: %s\n", gsl_strerror(status0));
        printError("AGM::unboundedCondition::A_DF", "status0 error");
    }

    if (status1) {
        printf("error: %s\n", gsl_strerror(status1));
        printError("AGM::unboundedCondition::A_DF", "status1 error");
    }

    gsl_integration_workspace_free(w0);
    gsl_integration_workspace_free(w1);

    return result1 / result0;
}

auto AGM::unboundedCondition::A_FD(const std::function<double(double, double)> &f) -> double {
    auto f0 = std::function<double(double, double)>{[this](double x, double y) -> double {
        return rho_m(x);
    }};
    auto f1 = std::function<double(double, double)>{[this, &f](double x, double y) -> double {
        return f(x, y) * rho_m(x);
    }};
    auto params0{unbounded_f_param{tm, tau, tp, tau0, &f0}};
    auto params1{unbounded_f_param{tm, tau, tp, tau0, &f1}};
    auto w0{gsl_integration_workspace_alloc(1000)};
    auto w1{gsl_integration_workspace_alloc(1000)};
    double result0{}, result1{}, error0{}, error1{};

    auto F0{gsl_function{unwrap, &params0}};
    auto F1{gsl_function{unwrap, &params1}};

    int status0{gsl_integration_qagil(&F0, tp, 0, 1e-10, 1000, w0, &result0, &error0)};
    int status1{gsl_integration_qagil(&F1, tp, 0, 1e-10, 1000, w1, &result1, &error1)};

    if (status0) {
        printf("error: %s\n", gsl_strerror(status0));
        printError("AGM::unboundedCondition::A_FD", "status0 error");
    }

    if (status1) {
        printf("error: %s\n", gsl_strerror(status1));
        printError("AGM::unboundedCondition::A_FD", "status1 error");
    }

    gsl_integration_workspace_free(w0);
    gsl_integration_workspace_free(w1);

    return result1 / result0;
}

AGM::unboundedCondition::~unboundedCondition() = default;
