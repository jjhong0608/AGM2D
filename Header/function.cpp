//
// Created by 조준홍 on 2022/02/13.
//

#include "function.h"

AGM::ellipticFunction::ellipticFunction() = default;

auto AGM::ellipticFunction::u(const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};

  auto test_problem = [&]() -> double {
    return std::sin(std::log(x + 2) * std::cos(y)) + 2;
  };

  auto lid_driven_cavity_flow = [&]() -> double {
    return ZEROVALUE;
  };

  auto bfs_flow = [&]() -> double {
    if (isclose(x, 150.)) {
      return (1.5e0 - y) * y * y;
    } else if (isclose(x, -1e1) && y > HALFVALUE) {
      return -2e0 * y * (4e0 * y * y - 9e0 * y + 6e0) + 2.5e0;
    } else if (isclose(y, UNITVALUE)) {
      return HALFVALUE;
    }
    return ZEROVALUE;
  };

  auto external_flow = [&]() -> double {
    return y * (1. - 0.25 / (x * x + y * y));
  };

  return test_problem();
}

auto AGM::ellipticFunction::phi(const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};
  return ZEROVALUE;
}

auto AGM::ellipticFunction::f(const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};
  return UNITVALUE;
}

auto AGM::ellipticFunction::ux(const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};
  return ZEROVALUE;
}

auto AGM::ellipticFunction::uy(const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};
  return ZEROVALUE;
}

auto AGM::ellipticFunction::isAssignBoundaryValue() -> bool {
  return false;
}

void AGM::ellipticFunction::assignBoundaryValue(AGM::point &pt) {
  if (!isAssignBoundaryValue()) {
    return;
  }
  auto x{pt[0]}, y{pt[1]};

  if (pt.getCondition() == 'D') {
    pt["bdv"] = u(pt);
  } else if (pt.getCondition() == 'N') {
    pt["bdv"] = ux(pt) * pt.getNormal()[0] + uy(pt) * pt.getNormal()[1];
  }
  pt["rhs"] = f(pt);
}

auto AGM::ellipticFunction::resultPath() -> std::string {
  return std::string{"/home/jhjo/AGM_AT_NIMS/Results/Elliptic_L_shape.result"};
}

auto AGM::ellipticFunction::FlowDataPath() -> std::string {
  return std::string{"/home/jhjo/extHD2/BFS/extension/re1600/AGM_Results_282.400000000"};
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

auto AGM::heatFunction::isLoadPreviousValue() -> bool {
  return false;
}

auto AGM::heatFunction::isAssignBoundaryValue() -> bool {
  return false;
}

void AGM::heatFunction::assignPreviousValue(AGM::value &value, AGM::point &pt) {
  if (!isLoadPreviousValue()) {
    return;
  }
  value["sol"] = u(pointHeat::getTime(), pt);
  value["phi"] = phi(pointHeat::getTime(), pt);
  value["rhs"] = f(pointHeat::getTime(), pt);
  value["dx"] = ux(pointHeat::getTime(), pt);
  value["dy"] = uy(pointHeat::getTime(), pt);
}

void AGM::heatFunction::assignBoundaryValue(AGM::point &pt) {
  if (!isAssignBoundaryValue()) {
    return;
  }
  if (pt.getCondition() == 'D') {
    pt["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), pt);
  } else if (pt.getCondition() == 'N') {
    pt["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[0] + uy(pointHeat::getTime() + pointHeat::getDelta(), pt) * pt.getNormal()[1];
  }
  pt["rhs"] = f(pointHeat::getTime() + pointHeat::getDelta(), pt);
}

auto AGM::heatFunction::resultPath() -> std::string {
  return std::string{"/home/jhjo/extHD2/BFS/extension/AGM_Result_stream_re1600"};
}

AGM::heatFunction::~heatFunction() = default;

AGM::NavierStokesFunction::NavierStokesFunction() = default;

auto AGM::NavierStokesFunction::initialTime() -> double {
  return ZEROVALUE;
}

auto AGM::NavierStokesFunction::terminalTime() -> double {
  return 50. + deltaTime();
}

auto AGM::NavierStokesFunction::deltaTime() -> double {
  return 0.01;
}

auto AGM::NavierStokesFunction::writeTime() -> double {
  return 10.0;
}

auto AGM::NavierStokesFunction::u(double t, const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};

  auto tesla_valve = [&](std::string &direction, std::string &stage) -> double {
    if (stage == "N4") {
      if (direction == "forward") {
        auto ymin{0.725044}, ymax{0.825044};
        return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
      } else if (direction == "reverse") {
        auto ymin{-0.398711}, ymax{-0.298711};
        return 6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
      }
    } else if (stage == "N2") {
      if (direction == "forward") {
        auto ymin{0.163044}, ymax{0.263044};
        return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
      } else if (direction == "reverse") {
        auto ymin{-0.398711}, ymax{-0.298711};
        return 6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
      }
    } else if (stage == "N1") {
      if (direction == "forward") {
        auto ymin{ZEROVALUE}, ymax{0.1};
        return -6. * (y - ymin) * (ymax - y) / std::pow(ymax - ymin, 3.);
      } else if (direction == "reverse") {
        auto xmin{-1.3413800000000002e+00}, xmax{-1.2706700000000002e+00};
        auto ymin{-8.4823000000000004e-01}, ymax{-7.7751999999999999e-01};
        auto d0{std::sqrt((x - xmin) * (x - xmin) + (y - ymax) * (y - ymax))};
        auto d1{std::sqrt((x - xmax) * (x - xmax) + (y - ymin) * (y - ymin))};
        auto dm{std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin))};
        return 6. * d0 * d1 / std::pow(dm, 3.) / std::sqrt(2);
      }
    }
    return ZEROVALUE;
  };

  auto kim_and_moin = [&]() -> double {
    return -std::cos(x) * std::sin(y) * std::exp(-2. * t);
  };

  auto taylor_green_vortex = [&]() -> double {
    auto Re{1000.};
    return -std::cos(2 * M_PI * x) * std::sin(2 * M_PI * y) * std::exp(-8 * std::pow(M_PI, 2) * t / Re);
  };

  auto galton_board = [&]() -> double {
    return -6. * (y - HALFVALUE) * (y + HALFVALUE);
  };

  auto lid_driven_cavity_flow = [&]() -> double {
    if (isclose(y, UNITVALUE)) {
      return UNITVALUE;
    }
    return ZEROVALUE;
  };

  auto bfs_flow = [&]() -> double {
    if (isclose(x, -1e1) && y > HALFVALUE) {
      return 2.4e1 * (y - HALFVALUE) * (UNITVALUE - y);
    }
    return ZEROVALUE;
  };

  auto kalman_vortex = [&]() -> double {
    auto a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
      return ZEROVALUE;
    }
    return UNITVALUE - pow(a, 2) / (pow(x, 2) + pow(y, 2)) + 2 * pow(a * y, 2) / pow(pow(x, 2) + pow(y, 2), 2);
  };

  auto axisymmetric_poiseuille = [&]() -> double {
    return ZEROVALUE;
  };

  auto axisymmetric_externalflow = [&]() -> double {
    auto U{UNITVALUE}, R{HALFVALUE};
    auto r{std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0))}, theta{std::atan2(x, y)};
    auto u_r = [&](double r, double theta) -> double {
      return -U * std::cos(theta) * (UNITVALUE - 3.0 * R / 2.0 / r + std::pow(R, 3.0) / 2.0 / std::pow(r, 3.0));
    };
    auto u_theta = [&](double r, double theta) -> double {
      return U * std::sin(theta) * (UNITVALUE - 3.0 * R / 4.0 / r - std::pow(R, 3.0) / 4.0 / std::pow(r, 3.0));
    };
    return u_r(r, theta) * std::sin(theta) + u_theta(r, theta) * std::cos(theta);
  };

  auto heart_depiction = [&]() -> double {
    if (iszero(x)) {
      if (ispositivezero(y)) {
        return 5.0 * y * (4.0 - y);
      } else {
        return 5.0 * y * (4.0 + y);
      }
    } else {
      return ZEROVALUE;
    }
  };

  return lid_driven_cavity_flow();
}

auto AGM::NavierStokesFunction::v(double t, const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};
  auto Re{1000.};

  auto tesla_valve = [&]() -> double {
    auto xmin{-1.3413800000000002e+00}, xmax{-1.2706700000000002e+00};
    auto ymin{-8.4823000000000004e-01}, ymax{-7.7751999999999999e-01};
    auto d0{std::sqrt((x - xmin) * (x - xmin) + (y - ymax) * (y - ymax))};
    auto d1{std::sqrt((x - xmax) * (x - xmax) + (y - ymin) * (y - ymin))};
    auto dm{std::sqrt((xmax - xmin) * (xmax - xmin) + (ymax - ymin) * (ymax - ymin))};
    return 6. * d0 * d1 / std::pow(dm, 3.) / std::sqrt(2);
  };

  auto kim_and_moin = [&]() -> double {
    return std::sin(x) * std::cos(y) * std::exp(-2. * t);
  };

  auto taylor_green_vortex = [&]() -> double {
    return std::sin(2 * M_PI * x) * std::cos(2 * M_PI * y) * std::exp(-8 * std::pow(M_PI, 2) * t / Re);
  };

  auto kalman_vortex = [&]() -> double {
    double a{HALFVALUE};
    if (pow(x, 2) + pow(y, 2) < pow(0.51, 2)) {
      return ZEROVALUE;
    }
    return -2 * pow(a, 2) * x * y / pow(pow(x, 2) + pow(y, 2), 2);
  };

  auto zero = [&]() -> double {
    return ZEROVALUE;
  };

  auto axisymmetric_poiseuille = [&]() -> double {
    return (HALFVALUE + x) * (HALFVALUE - x);
  };
  auto axisymmetric_externalflow = [&]() -> double {
    auto U{UNITVALUE}, R{HALFVALUE};
    auto r{std::sqrt(std::pow(x, 2.0) + std::pow(y, 2.0))}, theta{std::atan2(x, y)};
    auto u_r = [&](double r, double theta) -> double {
      return -U * std::cos(theta) * (UNITVALUE - 3.0 * R / 2.0 / r + std::pow(R, 3.0) / 2.0 / std::pow(r, 3.0));
    };
    auto u_theta = [&](double r, double theta) -> double {
      return U * std::sin(theta) * (UNITVALUE - 3.0 * R / 4.0 / r - std::pow(R, 3.0) / 4.0 / std::pow(r, 3.0));
    };
    return u_r(r, theta) * std::cos(theta) - u_theta(r, theta) * std::sin(theta);
  };

  auto heart_depiction = [&]() -> double {
    return ZEROVALUE;
  };

  return zero();
}

auto AGM::NavierStokesFunction::p(double t, const AGM::point &pt) -> double {
  auto x{pt[0]}, y{pt[1]};

  auto kim_and_moin = [&]() -> double {
    return -(std::cos(2. * x) + std::cos(2. * y)) * std::exp(-4. * t) / 4.;
  };

  auto taylor_green_vortex = [&]() -> double {
    auto Re{1000.};
    return -(std::cos(4 * M_PI * x) + std::cos(4 * M_PI * y)) / 4 * std::exp(-16 * std::pow(M_PI, 2) * t / Re);
  };

  auto zero = [&]() -> double {
    return ZEROVALUE;
  };
  return zero();
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

void AGM::NavierStokesFunction::assignPreviousValue(AGM::value &pu, AGM::value &pv, AGM::value &pp, point &uvel, point &vvel, point &pres) {
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

auto AGM::NavierStokesFunction::isAssignBoundaryValue() -> bool {
  return true;
}

void AGM::NavierStokesFunction::assignBoundaryValue(AGM::point &uvel, AGM::point &vvel, int presentIter) {
  if (!isAssignBoundaryValue()) {
    return;
  }
  if (presentIter == 0) {
    if (uvel.getCondition() == 'D') {
      uvel["bdv"] = u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    } else if (uvel.getCondition() == 'N') {
      uvel["bdv"] = ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] + uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1];
    }
    //        uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    if (vvel.getCondition() == 'D') {
      vvel["bdv"] = v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else if (vvel.getCondition() == 'N') {
      vvel["bdv"] = vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] + vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1];
    }
    //        vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
  } else {
    if (uvel.getCondition() == 'D') {
      uvel["bdv"] = u(pointHeat::getTime(), uvel) + u(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    } else if (uvel.getCondition() == 'N') {
      uvel["bdv"] = ux(pointHeat::getTime(), uvel) * uvel.getNormal()[0] + uy(pointHeat::getTime(), uvel) * uvel.getNormal()[1] + ux(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[0] + uy(pointHeat::getTime() + pointHeat::getDelta(), uvel) * uvel.getNormal()[1];
    }
    //        uvel["rhs"] = f1(pointHeat::getTime() + pointHeat::getDelta(), uvel);
    if (vvel.getCondition() == 'D') {
      vvel["bdv"] = v(pointHeat::getTime(), vvel) + v(pointHeat::getTime() + pointHeat::getDelta(), vvel);
    } else if (vvel.getCondition() == 'N') {
      vvel["bdv"] = vx(pointHeat::getTime(), vvel) * vvel.getNormal()[0] + vy(pointHeat::getTime(), vvel) * vvel.getNormal()[1] + vx(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[0] + vy(pointHeat::getTime() + pointHeat::getDelta(), vvel) * vvel.getNormal()[1];
    }
    //        vvel["rhs"] = f2(pointHeat::getTime() + pointHeat::getDelta(), vvel);
  }
}

auto AGM::NavierStokesFunction::isLoadPreviousFile() -> bool {
  return false;
}

void AGM::NavierStokesFunction::loadPreviousValue(std::vector<AGM::value> *pu, std::vector<AGM::value> *pv, std::vector<AGM::value> *pp) {
  if (!isLoadPreviousFile()) {
    return;
  }
  auto idx{0u}, bc{0u};
  auto x{ZEROVALUE}, y{ZEROVALUE};
  auto const &filename{std::string("/home/jhjo/extHD2/heat_depiction/test/AGM_Results_50.00000000")};
  std::ifstream f(filename);
  if (f.fail()) {
    printError("AGM::NavierStokesFunction::loadPreviousValue", "file %s is not opened", filename.c_str());
  }
  printf("Previous file: %s is opened\n", filename.c_str());
  for (int i = 0; i < point::getNPts(); ++i) {
    f >> idx >> x >> y;
    if (idx > pu->size()) {
      printError("AGM::NavierStokesFunction::loadPreviousValue", "idx (which is %d) is greater(or equal) than size of the point (which is %d)", idx, pu->size());
    }
    f >> pu->at(idx)["sol"];// u
    f >> pv->at(idx)["sol"];// v
    f >> pp->at(idx)["phi"];// p
    f >> pu->at(idx)["dx"]; // ux
    f >> pu->at(idx)["dy"]; // uy
    f >> pv->at(idx)["dx"]; // vx
    f >> pv->at(idx)["dy"]; // vy
    f >> pp->at(idx)["dx"]; // px
    f >> pp->at(idx)["dy"]; // py
    f >> pu->at(idx)["phi"];// phi
    f >> pv->at(idx)["phi"];// psi
    f >> pp->at(idx)["sol"];// phi of p
    f >> bc;
  }
  f.close();
}

auto AGM::NavierStokesFunction::isNormalEq() -> bool {
  return true;
}

auto AGM::NavierStokesFunction::findFixedPointIndex(std::vector<AGM::point> *pts) -> int {
  if (isNormalEq()) {
    for (const auto &item : *pts) {
      auto x{item[0]}, y{item[1]}
      if (isclose(x, 0.5) && isclose(y, 0.5)) {
        printf("Find Fixed Point at (%f, %f)\n", item[0], item[1]);
        return item.getIdx();
      }
    }
    printError("AGM::NavierStokesFunction::findFixedPointIndex", "Not Found Fixed Point Index");
  }
  return -1;
}

auto AGM::NavierStokesFunction::tolerance() -> double {
  return 1e-6;
}

auto AGM::NavierStokesFunction::resultPath() -> std::string {
  return {"/home/jhjo/AGM_AT_NIMS/Results/Lid-driven_cavity_re100.result"};
}

AGM::NavierStokesFunction::~NavierStokesFunction() = default;
