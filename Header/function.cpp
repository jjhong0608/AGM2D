//
// Created by 조준홍 on 2022/02/13.
//

#include "function.h"

AGM::ellipticFunction::ellipticFunction() = default;

auto AGM::ellipticFunction::u(const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};

    // test problem
    return std::pow(x, 3.) - 2. * std::pow(x, 2.) - std::pow(y, 2.) + 3. * y + 2.;

    // for Lid-driven cavity flow
    return ZEROVALUE;

    // for BFS flow
    if (isclose(x, 3e1)) {
        return (1.5e0 - y) * y * y;
    } else if (isclose(x, -1e1) && y > HALFVALUE) {
        return -2e0 * y * (4e0 * y * y - 9e0 * y + 6e0) + 2.5e0;
    } else if (isclose(y, UNITVALUE)) {
        return HALFVALUE;
    }
    return ZEROVALUE;
//    return UNITVALUE / std::sqrt(x * x + y * y);
}

auto AGM::ellipticFunction::phi(const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::ellipticFunction::f(const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    // test problem
    return -6. * x + 6.;
    return ZEROVALUE;
}

auto AGM::ellipticFunction::ux(const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    // test problem
    return 3. * std::pow(x, 2.) - 4. * x;
    return ZEROVALUE;
}

auto AGM::ellipticFunction::uy(const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    // test problem
    return -2. * y + 3.;
    return ZEROVALUE;
}

void AGM::ellipticFunction::assignBoundaryValue(AGM::point &pt) {
//    return;
    double x{pt[0]}, y{pt[1]};
//    if (iszero(x)) {
//        pt.setCondition('N');
//    } else if (isclose(y, 130.0)) {
//        pt.setCondition('N');
//        if (isclose(x, 31.5)) {
//            pt.setCondition('n');
//        } else if (isclose(x, 42.0)) {
//            pt.setCondition('n');
//        }
//    } else if (isclose(y, -194.738)) {
//        pt.setCondition('N');
//    } else if (x < 1e1 && y < 1e-6) {
//        pt["bdv"] = 1e2;
//    }

    if (pt.getCondition() == 'D') {
        pt["bdv"] = u(pt);
    } else if (pt.getCondition() == 'N') {
        pt["bdv"] = ux(pt) * pt.getNormal()[0] + uy(pt) * pt.getNormal()[1];
    }
    pt["rhs"] = f(pt);
}

AGM::ellipticFunction::~ellipticFunction() = default;

AGM::heatFunction::heatFunction() = default;

auto AGM::heatFunction::initialTime() -> double {
    return UNITVALUE;
}

auto AGM::heatFunction::terminalTime() -> double {
    return 1.25;
}

auto AGM::heatFunction::deltaTime() -> double {
    return 0.01;
}

auto AGM::heatFunction::u(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    double r{std::sqrt(std::pow(x, 2) + std::pow(y, 2))};
    return std::exp(-std::pow(r, 2) / (4 * t)) / (4 * M_PI * t);
}

auto AGM::heatFunction::phi(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return (1.0 / 32.0) * (-pow(x, 2) + pow(y, 2)) * exp(-1.0 / 4.0 * (pow(x, 2) + pow(y, 2)) / t) / (M_PI * pow(t, 3));
}

auto AGM::heatFunction::f(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::heatFunction::ux(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return -1.0 / 8.0 * x * exp(-1.0 / 4.0 * (pow(x, 2) + pow(y, 2)) / t) / (M_PI * pow(t, 2));
}

auto AGM::heatFunction::uy(double t, const AGM::point &pt) -> double {
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

auto AGM::NavierStokesFunction::initialTime() -> double {
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::terminalTime() -> double {
//    return 50.;
//    return 2.;
    return 1000.;
}

auto AGM::NavierStokesFunction::deltaTime() -> double {
    return 1e-3;
//    return 1e-5;
//    return (M_PI / 40.) * (M_PI / 40.) * 2.;
}

auto AGM::NavierStokesFunction::writeTime() -> double {
    return 1e-1;
    return 1e-2;
}

auto AGM::NavierStokesFunction::u(double t, const AGM::point &pt) -> double {
    auto x{pt[0]}, y{pt[1]};
    // Tesla valve (T45_R, N4, Forward)
//    auto ymin{0.725044}, ymax{0.825044};
//    return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);

    // Tesla valve (T45_R, N4, Reverse)
//    auto ymin{-0.398711}, ymax{-0.298711};
//    return 6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);


    // Tesla valve (T45_R, N2, Forward)
//    auto ymin{0.163044}, ymax{0.263044};
//    return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);

    // Tesla valve (T45_R, N2, Reverse)
//    auto ymin{-0.398711}, ymax{-0.298711};
//    return 6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);

    // Tesla valve (T45_R, Forward)
//    auto ymin{ZEROVALUE}, ymax{0.1};
//    return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);

    // Tesla valve (T45_R, Reverse)
//    auto xmin{-1.3413800000000002e+00}, xmax{-1.2706700000000002e+00};
//    auto ymin{-8.4823000000000004e-01}, ymax{-7.7751999999999999e-01};
//    auto d0{std::sqrt((x - xmin) * (x - xmin) + (y - ymax) * (y - ymax))};
//    auto d1{std::sqrt((x - xmax) * (x - xmax) + (y - ymin) * (y - ymin))};
//    auto dm{std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin))};
//    return 6. * d0 * d1 / std::pow(dm, 3.) / std::sqrt(2);

    // Tesla valve (left to right)
//    auto ymin{-1.9173022499999999e+00}, ymax{-1.8522147499999999e+00};
//    if (isclose(x, 1.2203700000000040e-02)) {
//        return 3. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
//    } else {
//        return ZEROVALUE;
//    }

    // Tesla valve (right to left)
//    auto ymin{-1.9157432600000002e+00}, ymax{-1.8470397900000002e+00};
//    if (isclose(x, 2.9891414400000000e+00)) {
//        return -3. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
//    } else {
//        return ZEROVALUE;
//    }

    // Tesla valve cylinder
//    auto ymin{-1.9173022499999999e+00}, ymax{-1.8522147499999999e+00};
//    auto cx{0.5}, cy{HALFVALUE * (ymax + ymin)}, r{HALFVALUE * HALFVALUE * (ymax - ymin) + 1e-3};
//    auto xmin{0.49 - NEARZERO}, ymin{-1.89 - NEARZERO};
//    auto xmax{0.51 + NEARZERO}, ymax{-1.88 + NEARZERO};
//    if ((x - cx) * (x - cx) + (y - cy) * (y - cy) < r * r) {
//    if ((xmin < x) && (x < xmax) && (ymin < y) && (y < ymax)) {
//        return 10. * std::sin(1e4 * t * M_PI);
//    } else if ((xmin < x - 2.) && (x - 2. < xmax) && (ymin < y) && (y < ymax)) {
//        return 10. * std::sin(1e4 * t * M_PI);
//    } else {
//        return ZEROVALUE;
//    }

// two-square cylinders
//    if (isclose(x, -7.5) | isclose(y, -7.5) | isclose(y, 7.5)) {
//        return UNITVALUE;
//    } else {
//        return ZEROVALUE;
//    }

// Kin and Moin
//    return -std::cos(x) * std::sin(y) * std::exp(-2. * t);

// FSI
//    return isclose(y, UNITVALUE) ? UNITVALUE
//                                 : isclose(y, -UNITVALUE) ? -UNITVALUE
//                                                          : ZEROVALUE;

auto Re{50.};
// Taylor-Green vortex
//    return -std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y) * std::exp(-8 * std::pow(M_PI, 2) * t / Re);

// Lid-driven cavity
//    if (isclose(y, UNITVALUE)) {
//        return UNITVALUE;
//    }
//    return ZEROVALUE;

// BFS flow
    if (isclose(x, -1e1) && y > HALFVALUE) {
        return 2.4e1 * (y - HALFVALUE) * (UNITVALUE - y);
    }
    return ZEROVALUE;

// Air Foil
//    if (isclose(x, -7)) {
//        return UNITVALUE;
//    } else if (isclose(y, -7)) {
//        return UNITVALUE;
//    } else if (isclose(y, 6.5)) {
//        return UNITVALUE;
//    }
//    return ZEROVALUE;

// Saccular aneurysm
//    double c{-3. / 32};
//    if (isclose(x, -12.)) {
//        return c * y * (y - 4.);
//    }
//    return ZEROVALUE;
// Kalman vortex
double a{HALFVALUE};
if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
    return ZEROVALUE;
}
return UNITVALUE - pow(a, 2) / (pow(x, 2) + pow(y, 2)) + 2 * pow(a * y, 2) / pow(pow(x, 2) + pow(y, 2), 2);
}

auto AGM::NavierStokesFunction::v(double t, const AGM::point &pt) -> double {
    auto x{pt[0]}, y{pt[1]};
    auto Re{50.};
    // two-square cylinders
    return ZEROVALUE;

    // Tesla valve (T45_R, Reverse)
//    auto xmin{-1.3413800000000002e+00}, xmax{-1.2706700000000002e+00};
//    auto ymin{-8.4823000000000004e-01}, ymax{-7.7751999999999999e-01};
//    auto d0{std::sqrt((x - xmin) * (x - xmin) + (y - ymax) * (y - ymax))};
//    auto d1{std::sqrt((x - xmax) * (x - xmax) + (y - ymin) * (y - ymin))};
//    auto dm{std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin))};
//    return 6. * d0 * d1 / std::pow(dm, 3.) / std::sqrt(2);

    // Kim and Moin
//    return std::sin(x) * std::cos(y) * std::exp(-2. * t);

    // Taylor-Green vortex
//    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::exp(-8 * std::pow(M_PI, 2) * t / Re);

    return ZEROVALUE;
    // Kalman vortex
    double a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
        return ZEROVALUE;
    }
    return -2 * pow(a, 2) * x * y / pow(pow(x, 2) + pow(y, 2), 2);
}

auto AGM::NavierStokesFunction::p(double t, const AGM::point &pt) -> double {
    auto x{pt[0]}, y{pt[1]};
    auto Re{50.};
    // two-sqaure cylinders
    return ZEROVALUE;

    // Kim and Moin
//    return -(std::cos(2. * x) + std::cos(2. * y)) * std::exp(-4. * t) / 4.;

    // Taylor-Green vortex
//    return -(std::cos(4 * M_PI * x) + std::cos(4 * M_PI * y)) / 4 * std::exp(-16 * std::pow(M_PI, 2) * t / Re);

    //
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::phi(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::psi(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::ux(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::uy(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::vx(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::vy(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::px(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::py(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::f1(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

auto AGM::NavierStokesFunction::f2(double t, const AGM::point &pt) -> double {
    double x{pt[0]}, y{pt[1]};
    return ZEROVALUE;
}

void
AGM::NavierStokesFunction::assignPreviousValue(
        AGM::value &pu,
        AGM::value &pv,
        AGM::value &pp,
        point &uvel,
        point &vvel,
        point &pres
) {
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

void AGM::NavierStokesFunction::assignBoundaryValue(
        AGM::point &uvel,
        AGM::point &vvel,
        int presentIter
) {
    if (presentIter == 0) {
        if (uvel.getCondition() == 'D') {
            uvel["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
        } else if (uvel.getCondition() == 'N') {
            uvel["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] +
                          uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1];
        }
//        uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
        if (vvel.getCondition() == 'D') {
            vvel["bdv"] = v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
        } else if (vvel.getCondition() == 'N') {
            vvel["bdv"] = vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] +
                          vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1];
        }
//        vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else {
        if (uvel.getCondition() == 'D') {
            uvel["bdv"] = u(pointHeat::getTime(), uvel) + u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
        } else if (uvel.getCondition() == 'N') {
            uvel["bdv"] = ux(pointHeat::getTime(), uvel) * uvel.getNormal()[0] +
                          uy(pointHeat::getTime(), uvel) * uvel.getNormal()[1] +
                          ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] +
                          uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1];
        }
//        uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
        if (vvel.getCondition() == 'D') {
            vvel["bdv"] = v(pointHeat::getTime(), vvel) + v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
        } else if (vvel.getCondition() == 'N') {
            vvel["bdv"] = vx(pointHeat::getTime(), vvel) * vvel.getNormal()[0] +
                          vy(pointHeat::getTime(), vvel) * vvel.getNormal()[1] +
                          vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] +
                          vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1];
        }
//        vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    }
}

void AGM::NavierStokesFunction::loadPreviousValue(
        const std::string &filename,
        std::vector<AGM::value> *pu,
        std::vector<AGM::value> *pv,
        std::vector<AGM::value> *pp
) {
    int idx{}, bc{};
    double x{}, y{};
    std::ifstream f(filename);
    if (f.fail()) {
        printError("AGM::NavierStokesFunction::loadPreviousValue", "file %s is not opened", filename.c_str());
    }
    std::cout << "previous file: " << filename << " is opened\n";
    for (int i = 0; i < point::getNPts(); ++i) {
        f >> idx >> x >> y;
        if (idx > pu->size()) {
            printError("AGM::NavierStokesFunction::loadPreviousValue",
                       "idx (which is %d) is greater(or equal) then size of the point (which is %d)", idx, pu->size());
        }
        f >> pu->at(idx)["sol"]; // u
        f >> pv->at(idx)["sol"]; // v
        f >> pp->at(idx)["sol"]; // p
        f >> pu->at(idx)["dx"]; // ux
        f >> pu->at(idx)["dy"]; // uy
        f >> pv->at(idx)["dx"]; // vx
        f >> pv->at(idx)["dy"]; // vy
//        f >> pp->at(idx)["dx"]; // px
//        f >> pp->at(idx)["dy"]; // py
        f >> pu->at(idx)["phi"]; // phi
        f >> pv->at(idx)["phi"]; // psi
//        f >> pp->at(idx)["phi"]; // phi of p
        f >> bc;
    }
    f.close();
}

AGM::NavierStokesFunction::~NavierStokesFunction() = default;
