//
// Created by 조준홍 on 2022/02/13.
//

#include "function.h"

AGM::ellipticFunction::ellipticFunction() = default;

double AGM::ellipticFunction::u(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return sin(x) + cos(y) + y;
}

double AGM::ellipticFunction::phi(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return sin(x) - HALFVALUE * f(pt);
}

double AGM::ellipticFunction::f(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
//    return ZEROVALUE;
    return sin(x) + cos(y);
}

double AGM::ellipticFunction::ux(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return cos(x);
}

double AGM::ellipticFunction::uy(const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return -sin(y) + UNITVALUE;
}

void AGM::ellipticFunction::assignBoundaryValue(AGM::point &pt) {
//    return;
    if (pt.getCondition() == 'D') {
        pt["bdv"] = u(pt);
    } else if (pt.getCondition() == 'N') {
        pt["bdv"] = ux(pt) * pt.getNormal()[0] + uy(pt) * pt.getNormal()[1];
    }
    pt["rhs"] = f(pt);
}

AGM::ellipticFunction::~ellipticFunction() = default;

AGM::heatFunction::heatFunction() = default;

double AGM::heatFunction::initialTime() {
    return UNITVALUE;
}

double AGM::heatFunction::terminalTime() {
    return 1.25;
}

double AGM::heatFunction::deltaTime() {
    return 0.01;
}

double AGM::heatFunction::u(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    double r{std::sqrt(std::pow(x, 2) + std::pow(y, 2))};
    return std::exp(-std::pow(r, 2) / (4 * t)) / (4 * M_PI * t);
}

double AGM::heatFunction::phi(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return (1.0 / 32.0) * (-pow(x, 2) + pow(y, 2)) * exp(-1.0 / 4.0 * (pow(x, 2) + pow(y, 2)) / t) / (M_PI * pow(t, 3));
}

double AGM::heatFunction::f(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

double AGM::heatFunction::ux(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return -1.0 / 8.0 * x * exp(-1.0 / 4.0 * (pow(x, 2) + pow(y, 2)) / t) / (M_PI * pow(t, 2));
}

double AGM::heatFunction::uy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return -1.0 / 8.0 * y * exp(-1.0 / 4.0 * (pow(x, 2) + pow(y, 2)) / t) / (M_PI * pow(t, 2));
}

void AGM::heatFunction::assignPreviousValue(AGM::value &value, AGM::point &pt) {
    value["sol"] = u(pointHeat::getTime(), pt);
    value["phi"] = phi(pointHeat::getTime(), pt);
    value["rhs"] = f(pointHeat::getTime(), pt);
    value["dx"] = ux(pointHeat::getTime(), pt);
    value["dy"] = uy(pointHeat::getTime(), pt);
}

void AGM::heatFunction::assignBoundaryValue(AGM::point &pt) {
    if (pt.getCondition() == 'D') {
        pt["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), pt);
    } else if (pt.getCondition() == 'N') {
        pt["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[0] +
                    uy(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[1];
    }
    pt["rhs"] = f(pointHeat::getTime() + pointHeat::getDelta(), pt);
}

AGM::heatFunction::~heatFunction() = default;

AGM::NavierStokesFunction::NavierStokesFunction() = default;

double AGM::NavierStokesFunction::initialTime() {
    return ZEROVALUE;
}

double AGM::NavierStokesFunction::terminalTime() {
    return 1.5e2;
}

double AGM::NavierStokesFunction::deltaTime() {
    return 1e-2;
}

double AGM::NavierStokesFunction::u(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    // Kalman vortex
    double a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
        return ZEROVALUE;
    }
    return UNITVALUE - pow(a, 2) / (pow(x, 2) + pow(y, 2)) + 2 * pow(a * y, 2) / pow(pow(x, 2) + pow(y, 2), 2);
}

double AGM::NavierStokesFunction::v(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    // Kalman vortex
    double a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
        return ZEROVALUE;
    }
    return -2 * pow(a, 2) * x * y / pow(pow(x, 2) + pow(y, 2), 2);
}

double AGM::NavierStokesFunction::p(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::phi(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::psi(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::ux(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::uy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::vx(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::vy(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::px(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::py(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::f1(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

double AGM::NavierStokesFunction::f2(double t, const AGM::point &pt) {
    double x{pt[0]}, y{pt[1]};
    return 0;
}

void
AGM::NavierStokesFunction::assignPreviousValue(AGM::value &pu, AGM::value &pv, AGM::value &pp, point &uvel, point &vvel,
                                               point &pres) {
    pu["sol"] = u(pointHeat::getTime(), uvel);
    pv["sol"] = v(pointHeat::getTime(), vvel);
    pp["sol"] = p(pointHeat::getTime(), pres);
    pu["phi"] = phi(pointHeat::getTime(), uvel);
    pv["phi"] = psi(pointHeat::getTime(), vvel);
    pu["rhs"] = f1(pointHeat::getTime(), uvel);
    pv["rhs"] = f2(pointHeat::getTime(), vvel);
    pu["dx"] = ux(pointHeat::getTime(), uvel);
    pv["dx"] = vx(pointHeat::getTime(), vvel);
    pp["dx"] = px(pointHeat::getTime(), pres);
    pu["dy"] = uy(pointHeat::getTime(), uvel);
    pv["dy"] = vy(pointHeat::getTime(), vvel);
    pp["dy"] = py(pointHeat::getTime(), pres);
}

void AGM::NavierStokesFunction::assignBoundaryValue(AGM::point &uvel, AGM::point &vvel) {
    if (uvel.getCondition() == 'D') {
        uvel["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    } else if (uvel.getCondition() == 'N') {
        uvel["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] +
                      uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1];
    }
    uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    if (vvel.getCondition() == 'D') {
        vvel["bdv"] = v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else if (vvel.getCondition() == 'N') {
        vvel["bdv"] = vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] +
                      vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1];
    }
    vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
}

AGM::NavierStokesFunction::~NavierStokesFunction() = default;
