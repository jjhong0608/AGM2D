//
// Created by 조준홍 on 2022/02/19.
//

#include "solver.h"

#include "util.h"

AGM::solver::solver(std::vector<point> *pts) : pts(pts) {}

AGM::solver::~solver() = default;

auto AGM::solver::getPts() const -> std::vector<AGM::point> * {
  return pts;
}

void AGM::solver::setPts(std::vector<point> *vector) {
  solver::pts = vector;
}

void AGM::solver::ellipticSolver() {
  auto f{AGM::ellipticFunction()};
  auto rhsX = [&](int i) -> double {
    return HALFVALUE * f.f(pts->at(i));
  };
  auto rhsY = [&](int i) -> double {
    return HALFVALUE * f.f(pts->at(i));
  };
  auto rhsXp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto rhsYp = [&](int i) -> double {
    return ZEROVALUE;
  };
#pragma omp parallel for
  for (auto item = pts->begin(); item != pts->end(); ++item) {
    //        if (item->getCondition() == 'D') {
    //            item->setCondition('N');
    //        } else if (item->getCondition() == 'N') {
    //            item->setCondition('D');
    //        }
    f.assignBoundaryValue(*item);
    item->setOrd(1);
  }
  point::setNPts(int(pts->size()));

  //    auto row_array{std::array<matrixRow, 2>{}};
  //    auto row{matrixRow()};
  //    auto rb_array{std::array<double, 2>{}};

#pragma omp parallel for
  for (auto item = pts->begin(); item != pts->end(); ++item) {
    item->calculateRepresentationFormula(1);
    item->makeDerivatives();
    item->updateRightHandSide(rhsX, rhsY);
    item->updateRightHandSidePart(rhsXp, rhsXp);

    //        if (item->getCondition() == 'D') {
    //            row_array = item->getSolMatrixRow();
    //            row = matrixRow();
    //            row[item->getIdx() + point::getNPts()] = UNITVALUE;
    //            row_array[1] = row;
    //            rb_array = item->getRb();
    //            rb_array[1] = f.phi(*item);
    //            item->setSolMatrixRow(row_array);
    //            item->setRb(rb_array);
    //        }
  }

  auto matrix = AGM::matrix<point>(pts);
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();
  matrix.releaseMatrix();

#pragma omp parallel for
  for (auto item = pts->begin(); item != pts->end(); ++item) {
    item->calculateDerivatives(pts, rhsX, rhsY, rhsXp, rhsYp);
  }
#pragma omp parallel for
  for (auto item = pts->begin(); item != pts->end(); ++item) {
    item->approximateNaNDerivatives(pts);
  }

  auto wf{AGM::writeFile<point>(pts)};
  //    wf.writeResult("/home/jhjo/extHD1/tesla_valve/test/AGM_Result");
  std::cout << "Relative Error of the solution = " << wf.calculateError("sol") << "\n";
}

void AGM::solver::streamSolver() {
  auto f{AGM::ellipticFunction()};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<value>(point::getNPts())};
  auto vvel{std::vector<value>(point::getNPts())};
  auto rightbound{std::vector<double>()};
  auto righty{std::vector<double>()};

  auto rhsX = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto rhsY = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto rhsXp = [&](int i) -> double {
    return -vvel.at(i)["sol"];
  };
  auto rhsYp = [&](int i) -> double {
    return uvel.at(i)["sol"];
  };

  auto NS{std::ifstream(
      "/home/jhjo/extHD2/BFS/extension/re1600/AGM_Results_282.40000000")};
  if (NS.fail()) {
    printError("file is not opened");
  }
  std::cout << "file is opened\n";
  int idx{}, bc{}, bdnum{};
  double x{}, y{}, p{}, px{}, py{}, pphi{};
  auto assignBoundary = [&](int i) -> double {
    auto sum{ZEROVALUE};
    for (int j = 1; j < i + 1; ++j) {
      sum += HALFVALUE * (rightbound.at(j) + rightbound.at(j - 1)) * (righty.at(j) - righty.at(j - 1));
    }
    return sum;
  };

  rightbound.emplace_back(ZEROVALUE);
  righty.emplace_back(ZEROVALUE);

  std::for_each(pts->begin(), pts->end(), [&](point &pt) -> void {
    NS >> idx >> x >> y >> uvel.at(idx)["sol"] >> vvel.at(idx)["sol"] >> p >> uvel.at(idx)["dx"] >> uvel.at(idx)["dy"] >> vvel.at(idx)["dx"] >> vvel.at(idx)["dy"] >> px >> py >> uvel.at(idx)["phi"] >> vvel.at(idx)["phi"] >> pphi >> bc;
    if (isclose(x, 150.)) {
      rightbound.emplace_back(uvel.at(idx)["sol"]);
      righty.emplace_back(y);
    }
    if (pt.getCondition() == 'D' || pt.getCondition() == 'd') {
      pt["bdv"] = f.u(pt);
    } else if (pt.getCondition() == 'N') {
      pt.setCondition('D');
      pt["bdv"] = assignBoundary(++bdnum);
      //            pt["bdv"] = f.u(pt);
    }
    //        pt["bdv"] = ZEROVALUE;
    pt.setMp(UNITVALUE);
    idx++;
  });

  NS.close();

#pragma omp parallel for
  for (auto item = pts->begin(); item != pts->end(); ++item) {
    item->calculateRepresentationFormula(2);
    item->makeDerivatives();
    item->updateRightHandSide(rhsX, rhsY);
    item->updateRightHandSidePart(rhsXp, rhsYp);
  }

  auto matrix = AGM::matrix<point>(pts);
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();
  matrix.releaseMatrix();

  //    #pragma omp parallel for
  //    for (auto item = pts->begin(); item != pts->end(); ++item) {
  //        item->calculateDerivatives(pts, rhsX, rhsY, rhsXp, rhsYp);
  //    }
  //    #pragma omp parallel for
  //    for (auto item = pts->begin(); item != pts->end(); ++item) {
  //        item->approximateNaNDerivatives(pts);
  //    }

  auto wf{AGM::writeFile<point>(pts)};
  wf.writeResult(
      "/home/jhjo/extHD2/BFS/extension/AGM_Result_stream_re1600");
}

void AGM::solver::axisymmetricEllipticSolver() {
  auto f{AGM::ellipticFunction()};
  point::setNPts(int(pts->size()));
  auto ptsAxis{std::vector<pointAxisymmetric>(point::getNPts())};
  auto rhsX = [&](int i) -> double {
    return HALFVALUE * ptsAxis.at(i)[0] * f.f(ptsAxis.at(i));
  };
  auto rhsY = [&](int i) -> double {
    return HALFVALUE * ptsAxis.at(i)[0] * f.f(ptsAxis.at(i));
  };
  auto rhsXp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto rhsYp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto copyPointInformation = [this, &ptsAxis]() -> void {
    for (int j = 0; j < point::getNPts(); ++j) {
      ptsAxis.at(j).point::operator=(pts->at(j));
      ptsAxis.at(j).findStencil(&(pts->at(j).getElement()), &ptsAxis);
    }
  };
  auto assignBoundaryValue = [&f, &ptsAxis]() -> void {
#pragma omp parallel for
    for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
      f.assignBoundaryValue(*item);
    }
  };
  auto makeMatrix = [&rhsX, &rhsY, &rhsXp, &rhsYp, &ptsAxis]() -> void {
#pragma omp parallel for
    for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
      item->checkOnAxis();
      item->calculateRepresentationFormula(1);
      item->makeDerivatives();
      item->updateRightHandSide(rhsX, rhsY);
      item->updateRightHandSidePart(rhsXp, rhsYp);
    }
#pragma omp parallel for
    for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
      item->EquationOnAxis();
    }
  };
  auto calculateDifferentiation = [&rhsX, &rhsY, &rhsXp, &rhsYp, &ptsAxis]() -> void {
#pragma omp parallel for
    for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
      item->calculateDerivatives(&ptsAxis, rhsX, rhsY, rhsXp, rhsYp);
    }
#pragma omp parallel for
    for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
      item->approximateNaNDerivatives(&ptsAxis);
    }
  };

  copyPointInformation();
  assignBoundaryValue();
  makeMatrix();

  for (auto item = ptsAxis.begin(); item != ptsAxis.end(); ++item) {
    for (const auto &row : item->getSolMatrixRow()[0]) {
      if (std::isnan(row.value)) {
        std::cout << "(x, y) = (" << item->getXy()[0] << ", " << item->getXy()[1] << ")\n";
        std::cout << "condition = " << item->getCondition() << "\n";
        printError("find NaN0 value");
      }
    }
    for (const auto &row : item->getSolMatrixRow()[1]) {
      if (std::isnan(row.value)) {
        std::cout << "(x, y) = (" << item->getXy()[0] << ", " << item->getXy()[1] << ")\n";
        std::cout << "condition = " << item->getCondition() << "\n";
        printError("find NaN1 value1");
      }
    }
    if (std::isnan(item->getRb()[0])) {
      std::cout << "(x, y) = (" << item->getXy()[0] << ", " << item->getXy()[1] << ")\n";
      std::cout << "condition = " << item->getCondition() << "\n";
      printError("rb0 error");
    }
    if (std::isnan(item->getRb()[1])) {
      std::cout << "(x, y) = (" << item->getXy()[0] << ", " << item->getXy()[1] << ")\n";
      std::cout << "condition = " << item->getCondition() << "\n";
      printError("rb1 error");
    }
  }

  auto matrix = AGM::matrix<pointAxisymmetric>(&ptsAxis);
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();
  matrix.releaseMatrix();

  //    for (const auto &item: ptsAxis) {
  //        std::cout << "condition = " << item.getCondition() << ", " << "value = " << item["sol"] << "\n";
  //    }

  calculateDifferentiation();
  auto wf{AGM::writeFile<pointAxisymmetric>(&ptsAxis)};
  wf.writeResult("/home/jjhong0608/docker/AGM2D/Axisymmetric/AGM_Result");
  std::cout << "Relative L2-error = " << wf.calculateError("sol") << "\n";
}

void AGM::solver::heatSolver() {
  auto f{AGM::heatFunction()};
  point::setNPts(int(pts->size()));
  auto ptsHeat{std::vector<pointHeat>(point::getNPts())};
  auto previousValues{std::vector<value>(point::getNPts())};
  pointHeat::setTime(AGM::heatFunction::initialTime());
  pointHeat::setDelta(AGM::heatFunction::deltaTime());
  auto rhsX = [&](int i) -> double {
    return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) + previousValues.at(i)["phi"] + previousValues.at(i)["sol"] / pointHeat::getDelta();
  };
  auto rhsY = [&](int i) -> double {
    return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) - previousValues.at(i)["phi"] + previousValues.at(i)["sol"] / pointHeat::getDelta();
  };
  auto rhsXp = [&](int i) -> double {
    return -ptsHeat.at(i).getMp() * previousValues.at(i)["dx"];
  };
  auto rhsYp = [&](int i) -> double {
    return -ptsHeat.at(i).getMp() * previousValues.at(i)["dy"];
  };
  auto copyPointInformation = [this, &ptsHeat]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      ptsHeat.at(i).point::operator=(pts->at(i));
      ptsHeat.at(i).findStencil(&(pts->at(i).getElement()), &ptsHeat);
    }
  };
  auto assignInitial = [&f, &ptsHeat, &previousValues]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      f.assignPreviousValue(previousValues.at(i), ptsHeat.at(i));
      f.assignBoundaryValue(ptsHeat.at(i));
    }
  };
  auto assignBoundaryValue = [&f, &ptsHeat]() -> void {
#pragma omp parallel for
    for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
      f.assignBoundaryValue(*item);
    }
  };
  auto makeMatrix = [&rhsX, &rhsY, &rhsXp, &rhsYp, &ptsHeat]() -> void {
#pragma omp parallel for
    for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
      item->calculateRepresentationFormula(1);
      item->makeDerivatives();
      item->updateRightHandSide(rhsX, rhsY);
      item->updateRightHandSidePart(rhsXp, rhsYp);
    }
  };
  auto calculateDifferentiation = [&rhsX, &rhsY, &rhsXp, &rhsYp, &ptsHeat]() -> void {
#pragma omp parallel for
    for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
      item->calculateDerivatives(&ptsHeat, rhsX, rhsY, rhsXp, rhsYp);
    }
#pragma omp parallel for
    for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
      item->approximateNaNDerivatives(&ptsHeat);
    }
  };
  auto updateRHS = [&rhsX, &rhsY, &rhsXp, &rhsYp, &ptsHeat]() -> void {
#pragma omp parallel for
    for (auto item = ptsHeat.begin(); item != ptsHeat.end(); ++item) {
      item->updateRightHandSide(rhsX, rhsY);
      item->updateRightHandSidePart(rhsXp, rhsYp);
    }
  };
  auto updateValue = [&ptsHeat, &previousValues]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      previousValues.at(i) = ptsHeat.at(i).getValue();
    }
  };
  auto updateTime = []() -> void {
    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
    std::cout << "current time = " << pointHeat::getTime() << "\n";
  };
  std::cout << "# of the points = " << point::getNPts() << "\n";
  copyPointInformation();
  assignInitial();
  makeMatrix();
  auto matrix = AGM::matrix<pointHeat>(&ptsHeat);
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();
  calculateDifferentiation();
  updateValue();
  while (pointHeat::getTime() + pointHeat::getDelta() < AGM::heatFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
    updateTime();
    assignBoundaryValue();
    updateRHS();
    matrix.calculateMatrix();
    calculateDifferentiation();
    updateValue();
  }
  updateTime();
  matrix.releaseMatrix();
  for (int i = 0; i < point::getNPts(); ++i) {
    pts->at(i).setValue(ptsHeat.at(i).getValue());
  }
}

void AGM::solver::StokesSolver() {
  auto f{AGM::NavierStokesFunction()};
  int presentIter{};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointStokes>(point::getNPts())};
  auto vvel{std::vector<pointStokes>(point::getNPts())};
  auto pvel{std::vector<pointStokes>(point::getNPts())};
  auto puvel{std::vector<value>(point::getNPts())};
  auto pvvel{std::vector<value>(point::getNPts())};
  auto ppvel{std::vector<double>(point::getNPts())};
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto uRhsXp = [&](int i) -> double {
    return ppvel.at(i);
    //        return ppvel.at(i) + puvel.at(i)["sol"] * puvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return ZEROVALUE;
    //        return puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp = [&](int i) -> double {
    return ZEROVALUE;
    //        return puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp = [&](int i) -> double {
    return ppvel.at(i);
    //        return ppvel.at(i) + pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      pvel.at(i).point::operator=(pts->at(i));
      pvel.at(i).findStencil(&(pts->at(i).getElement()), &pvel);
    }
  };
  auto bd_idx{std::vector<int>()};
  auto assignInitial = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      if (isclose(pts->at(i).getValue()["bdv"], UNITVALUE)) {
        bd_idx.emplace_back(pts->at(i).getIdx());
      }
      uvel.at(i)["bdv"] = ZEROVALUE;
      vvel.at(i)["bdv"] = ZEROVALUE;
    }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), 0);
    }

    //        #pragma omp parallel for
    //        for (int i = 0; i < point::getNPts(); ++i) {
    //            f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    //        }
  };
  auto assignBoundaryValue = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(2);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      uvel.at(i).makeRepresentationFormulaStokes();
      vvel.at(i).calculateRepresentationFormula(2);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
      vvel.at(i).makeRepresentationFormulaStokes();
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculatePressure = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pvel.at(i)["sol"] = -HALFVALUE * pts->at(i).getMp() * (uvel.at(i).calculateRepresentationFormulaStokes(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp, 0) + vvel.at(i).calculateRepresentationFormulaStokes(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp, 1));
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pvel.at(i).approximateNaNPressure(&pvel);
    }
  };
  auto updateRHSVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
    }
  };
  auto updateValue = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pvel.at(i)["sol"] = ppvel.at(i) + 1.5 * (pvel.at(i)["sol"] - ppvel.at(i));
    }
  };
  auto updateValue1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      puvel.at(i)["sol"] = uvel.at(i)["sol"];
      pvvel.at(i)["sol"] = vvel.at(i)["sol"];
      ppvel.at(i) = pvel.at(i)["sol"];
    }
  };
  auto tolerance{1e-5};
  auto calculateError = [&]() -> double {
    ++presentIter;
    auto perror{ZEROVALUE}, maxperror{ZEROVALUE};
    auto maxpvel{ZEROVALUE};
    for (int i = 0; i < point::getNPts(); ++i) {
      perror = std::fabs(pvel.at(i)["sol"] - ppvel.at(i));
      if (maxpvel < pvel.at(i)["sol"]) {
        maxpvel = pvel.at(i)["sol"];
      }
      if (maxperror < perror) {
        maxperror = perror;
      }
    }
    auto err{maxperror / maxpvel};
    std::cout << presentIter << "-th iteration, "
              << "current error = [" << err << " / "
              << tolerance << "]\n";
    return err;
  };
  copyPointInformation();
  assignInitial();
  makeMatrixVelocity();
  auto matrixVelocity{AGM::matrixMulti<pointStokes>(&uvel, &vvel)};
  matrixVelocity.makeMatrix();
  matrixVelocity.factorizeMatrix();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity();
  calculatePressure();
  updateValue();
  while (calculateError() > tolerance) {
    updateValue1();
    updateRHSVelocity();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    calculatePressure();
    updateValue();
  }
  auto wf{writeFileMultiple<pointStokes, pointStokes, pointStokes>(&uvel, &vvel, &pvel)};
  wf.writeResult(
      "/home/jhjo/extHD2/galton_board/cut/test15/AGM_Stokes_Results");
  matrixVelocity.releaseMatrix();
}

void AGM::solver::StokesSolverFull() {
  auto f{AGM::NavierStokesFunction()};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointStokes>(point::getNPts())};
  auto vvel{std::vector<pointStokes>(point::getNPts())};
  auto pvel{std::vector<pointStokes>(point::getNPts())};
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto uRhsXp = [&](int i) -> double {
    return pvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsXp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsYp = [&](int i) -> double {
    return pvel.at(i)["sol"];
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      pvel.at(i).point::operator=(pts->at(i));
      pvel.at(i).findStencil(&(pts->at(i).getElement()), &pvel);
    }
  };
  auto bd_idx{std::vector<int>()};
  auto assignBoundaryValue = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      if (isclose(pts->at(i).getValue()["bdv"], UNITVALUE)) {
        bd_idx.emplace_back(pts->at(i).getIdx());
      }
      uvel.at(i)["bdv"] = ZEROVALUE;
      vvel.at(i)["bdv"] = ZEROVALUE;
    }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), 0);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(2);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).makeRepresentationFormulaStokesFull();
      uvel.at(i).updateRightHandSideFull(uRhsX, uRhsY);
      vvel.at(i).calculateRepresentationFormula(2);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).makeRepresentationFormulaStokesFull();
      vvel.at(i).updateRightHandSideFull(vRhsXp, vRhsY);
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  copyPointInformation();
  assignBoundaryValue();
  makeMatrixVelocity();
  auto matrix{AGM::matrixStokesNormal<pointStokes>(&uvel, &vvel, &pvel, 150000)};
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();
  calculateDifferentiationVelocity();
  auto wf{writeFileMultiple<pointStokes, pointStokes, pointStokes>(&uvel, &vvel, &pvel)};
  wf.writeResult(
      "/home/jhjo/extHD2/galton_board/cut/test6/AGM_Stokes_Results_test");
  matrix.releaseMatrix();
}

void AGM::solver::axisymmetricStokesSolverFull() {
  auto f{NavierStokesFunction()};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointAxisymmetricStokes>(point::getNPts())};
  auto vvel{std::vector<pointAxisymmetricStokes>(point::getNPts())};
  auto pvel{std::vector<pointAxisymmetricStokes>(point::getNPts())};
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * uvel.at(i)["rhs"];
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * vvel.at(i)["rhs"];
  };
  auto uRhsXp = [&](int i) -> double {
    return pvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsXp = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsYp = [&](int i) -> double {
    return pvel.at(i)["sol"];
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      pvel.at(i).point::operator=(pts->at(i));
      pvel.at(i).findStencil(&(pts->at(i).getElement()), &pvel);
    }
  };
  //    auto bd_idx{std::vector<int>{}};
  auto assignBoundaryValue = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      //            if (isclose(pts->at(i).getValue()["bdv"], UNITVALUE)) {
      //                bd_idx.emplace_back(pts->at(i).getIdx());
      //            }
      //            uvel.at(i)["bdv"] = ZEROVALUE;
      //            vvel.at(i)["bdv"] = ZEROVALUE;
      if (iszero(vvel.at(i).getXy()[0])) {
        vvel.at(i).setCondition('N');
      }

      f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
    //        for (const auto &item: bd_idx) {
    //            f.assignBoundaryValue(uvel.at(item), vvel.at(item), 0);
    //        }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula('r', 2);
      //            uvel.at(i).makeDerivatives('r');
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).makeRepresentationFormulaStokesFull('r');
      uvel.at(i).updateRightHandSideFull(uRhsX, uRhsY);
      vvel.at(i).calculateRepresentationFormula('z', 2);
      //            vvel.at(i).makeDerivatives('z');
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).makeRepresentationFormulaStokesFull('z');
      vvel.at(i).updateRightHandSideFull(vRhsXp, vRhsY);
    }

#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      vvel.at(i).checkOnAxis();
      vvel.at(i).EquationOnAxis();
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  copyPointInformation();
  assignBoundaryValue();
  makeMatrixVelocity();
  auto matrix{AGM::matrixStokesNormal<pointAxisymmetricStokes>(&uvel, &vvel, &pvel, 10000)};
  matrix.makeMatrix();
  matrix.factorizeMatrix();
  matrix.calculateMatrix();

  //    calculateDifferentiationVelocity();
  auto wf{writeFileMultiple<pointAxisymmetricStokes, pointAxisymmetricStokes, pointAxisymmetricStokes>(&uvel, &vvel, &pvel)};
  wf.writeResult(
      "/home/jhjo/extHD2/axisymmetric_stokes/test/AGM_Results");
  matrix.releaseMatrix();
}

void AGM::solver::NavierStokesSolver() {
  auto f{AGM::NavierStokesFunction()};
  int fixedPointIndex{-1}, presentIter{}, saveIter{};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointHeat>(point::getNPts())};
  auto vvel{std::vector<pointHeat>(point::getNPts())};
  auto pvel{std::vector<point>(point::getNPts())};
  auto puvel{std::vector<value>(point::getNPts())};
  auto pvvel{std::vector<value>(point::getNPts())};
  auto ppvel{std::vector<value>(point::getNPts())};
  auto tildeuvel{std::vector<value>(point::getNPts())};
  auto tildevvel{std::vector<value>(point::getNPts())};
  auto hatuvel{std::vector<value>(point::getNPts())};
  auto hatvvel{std::vector<value>(point::getNPts())};
  pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
  pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
  auto bd_idx{std::vector<int>()};
  auto old_uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsX = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto pRhsY = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto uRhsX1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto uRhsY1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto vRhsX1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto vRhsY1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto old_uRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp = [&](int i) -> double {
    return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto old_uRhsXp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] + 2. * pts->at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp1 = [&](int i) -> double {
    return 2. * hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] + 2. * pts->at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto uRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp = [&](int i) -> double {
    return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto pRhsXp = [&](int i) -> double {
    return uvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsYp = [&](int i) -> double {
    return vvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] - puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto uRhsYp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp1 = [&](int i) -> double {
    return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] - pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto findSaveIter = [&]() -> void {
    saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
    presentIter = int(std::floor((std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) / NavierStokesFunction::deltaTime() + 0.5));
    std::cout << "Initial iteration number = " << presentIter << ",\n"
              << "file will be saved every " << saveIter << " iterations\n";
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      //            if (isclose(pts->at(i)[0], HALFVALUE) && isclose(pts->at(i)[1], HALFVALUE)) {
      //                fixedPointIndex = i;
      //                std::cout << "fixed point index = " << fixedPointIndex << "\n";
      //            }
    }
    //        if (fixedPointIndex == -1) {
    //            std::cout << "Not Found Fixed Point Index\n";
    //            exit(1);
    //        }
  };
  auto assignInitial = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      if (isclose(pts->at(i).getValue()["bdv"], UNITVALUE)) {
        bd_idx.emplace_back(pts->at(i).getIdx());
      }
    }

#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pts->at(i).setMp(UNITVALUE);
      if (pts->at(i).getCondition() == 'D') {
        pts->at(i).setCondition('N');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'N') {
        pts->at(i).setCondition('D');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'd') {
        pts->at(i).setCondition('n');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'n') {
        pts->at(i).setCondition('d');
        pts->at(i)["bdv"] = ZEROVALUE;
      }
      //            f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
      //            f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), 0);
      puvel.at(item)["sol"] = uvel.at(item)["bdv"];
      pvvel.at(item)["sol"] = vvel.at(item)["bdv"];
    }
  };
  auto loadPreviousValue = [&]() -> void {
    f.loadPreviousValue(
        "/home/jhjo/extHD2/galton_board/cut/test6/AGM_Stokes_Results_test",
        &puvel, &pvvel, &ppvel);
  };
  auto assignBoundaryValue = [&]() -> void {
#pragma omp parallel for
    //        for (int i = 0; i < point::getNPts(); ++i) {
    //            f.assignBoundaryValue(uvel.at(i), vvel.at(i), presentIter);
    //        }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), presentIter);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(2);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateRepresentationFormula(2);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
    }
  };
  auto makeMatrixPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->setOrd(2);
      item->calculateRepresentationFormula(2);
      item->makeDerivatives();
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto calculateDifferentiationVelocityOld = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp, old_vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocityOld1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp1, old_uRhsYp1);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp1, old_vRhsYp1);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX1, uRhsY1, uRhsXp1, uRhsYp1);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX1, vRhsY1, vRhsXp1, vRhsYp1);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivatives(pts, pRhsX, pRhsY, pRhsXp, pRhsYp);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->approximateNaNDerivatives(pts);
    }
  };
  auto updateRHSVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
    }
  };
  auto updateRHSVelocityOld1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp1, old_uRhsYp1);
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp1, old_vRhsYp1);
    }
  };
  auto updateRHSVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX1, uRhsY1);
      uvel.at(i).updateRightHandSidePart(uRhsXp1, uRhsYp1);
      vvel.at(i).updateRightHandSide(vRhsX1, vRhsY1);
      vvel.at(i).updateRightHandSidePart(vRhsXp1, vRhsYp1);
    }
  };
  auto updateRHSPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto updateValues = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      tildeuvel.at(i) = uvel.at(i).getValue();
      tildevvel.at(i) = vvel.at(i).getValue();
      if (uvel.at(i).getCondition() != 'D') {
        uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
        vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
      }
      pts->at(i)["sol"] -= HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"]);
      hatuvel.at(i) = uvel.at(i).getValue();
      hatvvel.at(i) = vvel.at(i).getValue();
      uvel.at(i)["bdv"] = ZEROVALUE;
      vvel.at(i)["bdv"] = ZEROVALUE;
    }
  };
  auto updateValues1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      puvel.at(i) = uvel.at(i).getValue();
      pvvel.at(i) = vvel.at(i).getValue();
      pvel.at(i)["sol"] = HALFVALUE * (pts->at(i)["sol"] + ppvel.at(i)["sol"]);
      pvel.at(i)["phi"] = pts->at(i)["sol"];
      pvel.at(i)["dx"] = pts->at(i)["dx"];
      pvel.at(i)["dy"] = pts->at(i)["dy"];
      ppvel.at(i) = pts->at(i).getValue();
    }
  };
  auto subtractPreviousVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] -= puvel.at(i)["sol"];
      vvel.at(i)["sol"] -= pvvel.at(i)["sol"];
      uvel.at(i)["dx"] -= puvel.at(i)["dx"];
      vvel.at(i)["dx"] -= pvvel.at(i)["dx"];
      uvel.at(i)["dy"] -= puvel.at(i)["dy"];
      vvel.at(i)["dy"] -= pvvel.at(i)["dy"];
    }
  };
  auto addPreviousVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] += tildeuvel.at(i)["sol"];
      vvel.at(i)["sol"] += tildevvel.at(i)["sol"];
      uvel.at(i)["dx"] += tildeuvel.at(i)["dx"];
      vvel.at(i)["dx"] += tildevvel.at(i)["dx"];
      uvel.at(i)["dy"] += tildeuvel.at(i)["dy"];
      vvel.at(i)["dy"] += tildevvel.at(i)["dy"];
    }
  };
  auto stopCondition{ZEROVALUE};
  auto currentDiv{ZEROVALUE};
  auto checkStopCondition = [&]() -> bool {
    if (presentIter % saveIter != 0) {
      return false;
    }
    auto max_vel{ZEROVALUE}, max_err{ZEROVALUE}, mean_div{ZEROVALUE};
    auto vel{ZEROVALUE}, pvel{ZEROVALUE}, div{ZEROVALUE};
    auto tolerance{1e-6};
    auto size{0};
    for (int i = 0; i < point::getNPts(); ++i) {
      vel = std::sqrt(uvel.at(i)["sol"] * uvel.at(i)["sol"] + vvel.at(i)["sol"] * vvel.at(i)["sol"]);
      pvel = std::sqrt(puvel.at(i)["sol"] * puvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"]);
      div = std::abs(uvel.at(i)["dx"] + vvel.at(i)["dy"]);
      if (max_vel < vel) {
        max_vel = vel;
      }
      if (max_err < std::abs(vel - pvel)) {
        max_err = std::abs(vel - pvel);
      }
      if (uvel.at(i).getCondition() != 'D' && uvel.at(i).getCondition() != 'N') {
        mean_div += div;
        size++;
      }
    }
    stopCondition = max_err / pointHeat::getDelta() / max_vel;
    currentDiv = mean_div / size;
    //        stopCondition = max_err / max_vel;
    std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "], Div = "
              << currentDiv << "\n";
    if (stopCondition < tolerance) {
      return true;
    }
    return false;
  };
  auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, &pvel)};
  std::ifstream infile("/home/jhjo/extHD2/galton_board/cut/test6/Stopping_error.txt");
  if (infile.good()) {
    remove("/home/jhjo/extHD2/galton_board/cut/test6/Stopping_error.txt");
  }
  std::ofstream fp("/home/jhjo/extHD2/galton_board/cut/test6/Stopping_error.txt", std::ios_base::app);
  fp << "Reynolds number = " << 1. / pts->at(0).getMp() << "\n";
  fp.close();
  auto updateTime = [&]() -> void {
    ++presentIter;
    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
    std::cout << presentIter << "-th iteration, "
              << "current time = [" << pointHeat::getTime() << " / "
              << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]"
              << ", Current Div = [" << currentDiv << "]\n";

    if (presentIter % saveIter == 0) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(8) << pointHeat::getTime();
      std::string time_string = ss.str();
      wf.writeResult(
          "/home/jhjo/extHD2/galton_board/cut/test6/AGM_Results_"
          + time_string);
      std::ofstream fp("/home/jhjo/extHD2/galton_board/cut/test6/Stopping_error.txt", std::ios_base::app);
      fp << presentIter << "-th iteration, "
         << "current time = [" << pointHeat::getTime() << " / "
         << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]"
         << ", Current Div = [" << currentDiv << "]\n";
      fp.close();
    }
  };
  findSaveIter();
  copyPointInformation();
  assignInitial();
  makeMatrixVelocity();
  auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel)};
  auto matrixPressure{AGM::matrix<point>(pts)};
  //    auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPointIndex)};
  matrixVelocity.makeMatrix();
  matrixVelocity.factorizeMatrix();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocityOld();
  makeMatrixPressure();
  matrixPressure.makeMatrix();
  matrixPressure.factorizeMatrix();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  updateValues();
  updateRHSVelocity1();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity1();
  addPreviousVelocity();
  //    updateValues1();

  for (int i = 0; i < point::getNPts(); ++i) {
    hatuvel.at(i)["dx"] = tildeuvel.at(i)["dx"];
    hatuvel.at(i)["dy"] = tildeuvel.at(i)["dy"];
    hatvvel.at(i)["dx"] = tildevvel.at(i)["dx"];
    hatvvel.at(i)["dy"] = tildevvel.at(i)["dy"];
  }

  //--------
  ++presentIter;
  assignBoundaryValue();
  updateRHSVelocity();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity();
  subtractPreviousVelocity();
  updateRHSPressure();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  for (int i = 0; i < point::getNPts(); ++i) {
    uvel.at(i)["dx"] = hatuvel.at(i)["dx"];
    uvel.at(i)["dy"] = hatuvel.at(i)["dy"];
    vvel.at(i)["dx"] = hatvvel.at(i)["dx"];
    vvel.at(i)["dy"] = hatvvel.at(i)["dy"];
  }
  updateValues();
  for (int i = 0; i < point::getNPts(); ++i) {
    ppvel.at(i)["sol"] = pts->at(i)["sol"];
  }
  updateRHSVelocity1();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity1();
  addPreviousVelocity();
  updateValues1();
  --presentIter;
  //--------

  //    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
  //    for (int i = 0; i < point::getNPts(); ++i) {
  //        f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
  //    }
  //    pointHeat::setTime(pointHeat::getTime() - pointHeat::getDelta());
  loadPreviousValue();
  while (pointHeat::getTime() + pointHeat::getDelta() < AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
    updateTime();
    assignBoundaryValue();
    updateRHSVelocity();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    subtractPreviousVelocity();
    updateRHSPressure();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    updateRHSVelocity1();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity1();
    addPreviousVelocity();
    if (checkStopCondition()) {
      break;
    }
    updateValues1();
  }
  wf.writeResult("/home/jhjo/extHD2/galton_board/cut/test6/AGM_Results");
  updateTime();
  matrixVelocity.releaseMatrix();
  matrixPressure.releaseMatrix();
}

void AGM::solver::ModifyNavierStokesSolver() {
  auto f{AGM::NavierStokesFunction()};
  int fixedPointIndex{-1}, presentIter{}, saveIter{};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointHeat>(point::getNPts())};
  auto vvel{std::vector<pointHeat>(point::getNPts())};
  auto pvel{std::vector<point>(point::getNPts())};
  auto puvel{std::vector<value>(point::getNPts())};
  auto pvvel{std::vector<value>(point::getNPts())};
  auto ppvel{std::vector<value>(point::getNPts())};
  auto tildeuvel{std::vector<value>(point::getNPts())};
  auto tildevvel{std::vector<value>(point::getNPts())};
  auto hatuvel{std::vector<pointHeat>(point::getNPts())};
  auto hatvvel{std::vector<pointHeat>(point::getNPts())};
  auto phatuvel{std::vector<value>(point::getNPts())};
  auto phatvvel{std::vector<value>(point::getNPts())};
  pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
  pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
  auto bd_idx{std::vector<int>()};
  auto old_uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsX = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto pRhsY = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto uRhsX1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto uRhsY1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto vRhsX1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto vRhsY1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto old_uRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp = [&](int i) -> double {
    return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto old_uRhsXp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] + 2. * pts->at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp1 = [&](int i) -> double {
    return 2. * hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp1 = [&](int i) -> double {
    return 2. * hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] + 2. * pts->at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto uRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp = [&](int i) -> double {
    return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto pRhsXp = [&](int i) -> double {
    return uvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsYp = [&](int i) -> double {
    return vvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsX0 = [&](int i) -> double {
    return -uvel.at(i)["dx"] / pointHeat::getDelta();
  };
  auto pRhsY0 = [&](int i) -> double {
    return -vvel.at(i)["dy"] / pointHeat::getDelta();
  };
  auto uRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] - puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto uRhsYp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp1 = [&](int i) -> double {
    return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] - pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto findSaveIter = [&]() -> void {
    saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
    presentIter = int(std::floor(
        (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) / NavierStokesFunction::deltaTime() + 0.5));
    std::cout << "Initial iteration number = " << presentIter << ",\n"
              << "file will be saved every " << saveIter
              << " iterations\n";
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      hatuvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      hatvvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      //            if (isclose(pts->at(i)[0], HALFVALUE) && isclose(pts->at(i)[1], HALFVALUE)) {
      //                fixedPointIndex = i;
      //                std::cout << "fixed point index = " << fixedPointIndex << "\n";
      //            }
    }
    //        if (fixedPointIndex == -1) {
    //            std::cout << "Not Found Fixed Point Index\n";
    //            exit(1);
    //        }
  };
  auto assignInitial = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      if (isclose(pts->at(i).getValue()["bdv"], UNITVALUE)) {
        bd_idx.emplace_back(pts->at(i).getIdx());
      }
    }

#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pts->at(i).setMp(UNITVALUE);
      if (pts->at(i).getCondition() == 'D') {
        pts->at(i).setCondition('N');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'N') {
        pts->at(i).setCondition('D');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'd') {
        pts->at(i).setCondition('n');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'n') {
        pts->at(i).setCondition('d');
        pts->at(i)["bdv"] = ZEROVALUE;
      }
      //            f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
      //            f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), 0);
      puvel.at(item)["sol"] = uvel.at(item)["bdv"];
      pvvel.at(item)["sol"] = vvel.at(item)["bdv"];
    }
  };
  auto loadPreviousValue = [&]() -> void {
    f.loadPreviousValue(
        "/home/jhjo/extHD2/galton_board/test2/AGM_Results_8.37450000",
        &puvel, &pvvel, &ppvel);
  };
  auto assignBoundaryValue = [&]() -> void {
#pragma omp parallel for
    //        for (int i = 0; i < point::getNPts(); ++i) {
    //            f.assignBoundaryValue(uvel.at(i), vvel.at(i), presentIter);
    //        }
    for (const auto &item : bd_idx) {
      f.assignBoundaryValue(uvel.at(item), vvel.at(item), presentIter);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(2);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateRepresentationFormula(2);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
    }
  };
  auto makeMatrixPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->setOrd(2);
      item->calculateRepresentationFormula(2);
      item->makeDerivatives();
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto calculateDifferentiationVelocityOld = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp, old_vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocityOld1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp1, old_uRhsYp1);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp1, old_vRhsYp1);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX1, uRhsY1, uRhsXp1, uRhsYp1);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX1, vRhsY1, vRhsXp1, vRhsYp1);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivatives(pts, pRhsX, pRhsY, pRhsXp, pRhsYp);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivativesTwice(pRhsX0, pRhsY0);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->approximateNaNDerivatives(pts);
    }
  };
  auto updateRHSVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
    }
  };
  auto updateRHSVelocityOld1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp1, old_uRhsYp1);
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp1, old_vRhsYp1);
    }
  };
  auto updateRHSVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX1, uRhsY1);
      uvel.at(i).updateRightHandSidePart(uRhsXp1, uRhsYp1);
      vvel.at(i).updateRightHandSide(vRhsX1, vRhsY1);
      vvel.at(i).updateRightHandSidePart(vRhsXp1, vRhsYp1);
    }
  };
  auto updateRHSPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto updateValues = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      tildeuvel.at(i) = uvel.at(i).getValue();
      tildevvel.at(i) = vvel.at(i).getValue();
      if (uvel.at(i).getCondition() != 'D') {
        uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
        vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
      }
      pts->at(i)["sol"] -= HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"]);
      uvel.at(i)["dx"] -= pts->at(i)["dxx"] * pointHeat::getDelta();
      vvel.at(i)["dy"] -= pts->at(i)["dyy"] * pointHeat::getDelta();
      phatuvel.at(i) = hatuvel.at(i).getValue();
      phatvvel.at(i) = hatvvel.at(i).getValue();
      hatuvel.at(i).setValue(uvel.at(i).getValue());
      hatvvel.at(i).setValue(vvel.at(i).getValue());
      uvel.at(i)["bdv"] = ZEROVALUE;
      vvel.at(i)["bdv"] = ZEROVALUE;
    }
  };
  auto updateValues1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      puvel.at(i) = uvel.at(i).getValue();
      pvvel.at(i) = vvel.at(i).getValue();
      pvel.at(i)["sol"] = HALFVALUE * (pts->at(i)["sol"] + ppvel.at(i)["sol"]);
      pvel.at(i)["phi"] = pts->at(i)["sol"];
      pvel.at(i)["dx"] = pts->at(i)["dx"];
      pvel.at(i)["dy"] = pts->at(i)["dy"];
      ppvel.at(i) = pts->at(i).getValue();
      hatuvel.at(i)["dy"] = uvel.at(i)["dy"];
      hatvvel.at(i)["dx"] = vvel.at(i)["dx"];
    }
  };
  auto subtractPreviousVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] -= puvel.at(i)["sol"];
      vvel.at(i)["sol"] -= pvvel.at(i)["sol"];
      uvel.at(i)["dx"] -= puvel.at(i)["dx"];
      vvel.at(i)["dx"] -= pvvel.at(i)["dx"];
      uvel.at(i)["dy"] -= puvel.at(i)["dy"];
      vvel.at(i)["dy"] -= pvvel.at(i)["dy"];
    }
  };
  auto addPreviousVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] += tildeuvel.at(i)["sol"];
      vvel.at(i)["sol"] += tildevvel.at(i)["sol"];
      uvel.at(i)["dx"] += tildeuvel.at(i)["dx"];
      vvel.at(i)["dx"] += tildevvel.at(i)["dx"];
      uvel.at(i)["dy"] += tildeuvel.at(i)["dy"];
      vvel.at(i)["dy"] += tildevvel.at(i)["dy"];
    }
  };
  auto stopCondition{ZEROVALUE};
  auto currentDiv{ZEROVALUE};
  auto checkStopCondition = [&]() -> bool {
    if (presentIter % saveIter != 0) {
      return false;
    }
    auto max_vel{ZEROVALUE}, max_err{ZEROVALUE}, mean_div{ZEROVALUE};
    auto vel{ZEROVALUE}, pvel{ZEROVALUE}, div{ZEROVALUE};
    auto tolerance{1e-6};
    auto size{0};
    for (int i = 0; i < point::getNPts(); ++i) {
      vel = std::sqrt(hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] + hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"]);
      pvel = std::sqrt(
          phatuvel.at(i)["sol"] * phatuvel.at(i)["sol"] + phatvvel.at(i)["sol"] * phatvvel.at(i)["sol"]);
      div = std::abs(hatuvel.at(i)["dx"] + hatvvel.at(i)["dy"]);
      if (max_vel < vel) {
        max_vel = vel;
      }
      if (max_err < std::abs(vel - pvel)) {
        max_err = std::abs(vel - pvel);
      }
      if (uvel.at(i).getCondition() != 'D' && uvel.at(i).getCondition() != 'N') {
        mean_div += div;
        size++;
      }
    }
    stopCondition = max_err / pointHeat::getDelta() / max_vel;
    currentDiv = mean_div / size;
    //        stopCondition = max_err / max_vel;
    std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "], Div = "
              << currentDiv << "\n";
    if (stopCondition < tolerance) {
      return true;
    }
    return false;
  };
  auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&hatuvel, &hatvvel, &pvel)};
  //    auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, &pvel)};
  std::ifstream infile("/home/jhjo/extHD2/galton_board/test3/Stopping_error.txt");
  if (infile.good()) {
    remove("/home/jhjo/extHD2/galton_board/test3/Stopping_error.txt");
  }
  std::ofstream fp("/home/jhjo/extHD2/galton_board/test3/Stopping_error.txt", std::ios_base::app);
  fp << "Reynolds number = " << 1. / pts->at(0).getMp() << "\n";
  fp.close();
  auto updateTime = [&]() -> void {
    ++presentIter;
    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
    std::cout << presentIter << "-th iteration, "
              << "current time = [" << pointHeat::getTime() << " / "
              << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]"
              << ", Current Div = [" << currentDiv << "]\n";

    if (presentIter % saveIter == 0) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(8) << pointHeat::getTime();
      std::string time_string = ss.str();
      wf.writeResult(
          "/home/jhjo/extHD2/galton_board/test3/AGM_Results_"
          + time_string);
      std::ofstream fp("/home/jhjo/extHD2/galton_board/test3/Stopping_error.txt", std::ios_base::app);
      fp << presentIter << "-th iteration, "
         << "current time = [" << pointHeat::getTime() << " / "
         << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]"
         << ", Current Div = [" << currentDiv << "]\n";
      fp.close();
    }
  };
  findSaveIter();
  copyPointInformation();
  assignInitial();
  makeMatrixVelocity();
  auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel)};
  auto matrixPressure{AGM::matrix<point>(pts)};
  //    auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPointIndex)};
  matrixVelocity.makeMatrix();
  matrixVelocity.factorizeMatrix();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocityOld();
  makeMatrixPressure();
  matrixPressure.makeMatrix();
  matrixPressure.factorizeMatrix();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  updateValues();
  updateRHSVelocity1();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity1();
  addPreviousVelocity();
  updateValues1();

  for (int i = 0; i < point::getNPts(); ++i) {
    hatuvel.at(i)["dx"] = tildeuvel.at(i)["dx"];
    hatuvel.at(i)["dy"] = tildeuvel.at(i)["dy"];
    hatvvel.at(i)["dx"] = tildevvel.at(i)["dx"];
    hatvvel.at(i)["dy"] = tildevvel.at(i)["dy"];
  }

  //--------
  ++presentIter;
  assignBoundaryValue();
  updateRHSVelocity();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity();
  subtractPreviousVelocity();
  updateRHSPressure();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  for (int i = 0; i < point::getNPts(); ++i) {
    uvel.at(i)["dx"] = hatuvel.at(i)["dx"];
    uvel.at(i)["dy"] = hatuvel.at(i)["dy"];
    vvel.at(i)["dx"] = hatvvel.at(i)["dx"];
    vvel.at(i)["dy"] = hatvvel.at(i)["dy"];
  }
  updateValues();
  for (int i = 0; i < point::getNPts(); ++i) {
    ppvel.at(i)["sol"] = pts->at(i)["sol"];
  }
  updateRHSVelocity1();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity1();
  addPreviousVelocity();
  updateValues1();
  --presentIter;
  //--------

  //    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
  //    for (int i = 0; i < point::getNPts(); ++i) {
  //        f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
  //    }
  //    pointHeat::setTime(pointHeat::getTime() - pointHeat::getDelta());
  loadPreviousValue();
  while (pointHeat::getTime() + pointHeat::getDelta() < AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
    updateTime();
    assignBoundaryValue();
    updateRHSVelocity();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    subtractPreviousVelocity();
    updateRHSPressure();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    updateRHSVelocity1();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity1();
    addPreviousVelocity();
    if (checkStopCondition()) {
      break;
    }
    updateValues1();
  }
  wf.writeResult("/home/jhjo/extHD2/galton_board/test3/AGM_Results");
  updateTime();
  matrixVelocity.releaseMatrix();
  matrixPressure.releaseMatrix();
}

void AGM::solver::TwoStepNS() {
  auto f{AGM::NavierStokesFunction()};
  int fixedPointIndex{}, presentIter{}, saveIter{};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointHeat>(point::getNPts())};
  auto vvel{std::vector<pointHeat>(point::getNPts())};
  auto hatuvel{std::vector<pointHeat>(point::getNPts())};
  auto hatvvel{std::vector<pointHeat>(point::getNPts())};
  auto puvel{std::vector<value>(point::getNPts())};
  auto pvvel{std::vector<value>(point::getNPts())};
  auto ppvel{std::vector<value>(point::getNPts())};
  auto ppuvel{std::vector<value>(point::getNPts())};
  auto ppvvel{std::vector<value>(point::getNPts())};
  pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
  pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
  auto old_uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2. * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsX = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto pRhsY = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto old_uRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp = [&](int i) -> double {
    return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp = [&](int i) -> double {
    return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto uRhsXp = [&](int i) -> double {
    if (presentIter < 2) {
      return 2. * puvel.at(i)["sol"] * puvel.at(i)["sol"];
    } else {
      return 3. * puvel.at(i)["sol"] * puvel.at(i)["sol"] - ppuvel.at(i)["sol"] * ppuvel.at(i)["sol"];
    }
  };
  auto uRhsYp = [&](int i) -> double {
    if (presentIter < 2) {
      return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    } else {
      return 3. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - ppuvel.at(i)["sol"] * ppvvel.at(i)["sol"];
    }
  };
  auto vRhsXp = [&](int i) -> double {
    if (presentIter < 2) {
      return 2. * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    } else {
      return 3. * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - ppuvel.at(i)["sol"] * ppvvel.at(i)["sol"];
    }
  };
  auto vRhsYp = [&](int i) -> double {
    if (presentIter < 2) {
      return 2. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
    } else {
      return 3. * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - ppvvel.at(i)["sol"] * ppvvel.at(i)["sol"];
    }
  };
  auto pRhsXp = [&](int i) -> double {
    return uvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsYp = [&](int i) -> double {
    return vvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsXDiff = [&](int i) -> double {
    return -uvel.at(i)["dx"] / pointHeat::getDelta();
  };
  auto pRhsYDiff = [&](int i) -> double {
    return -vvel.at(i)["dy"] / pointHeat::getDelta();
  };
  auto uRhsXPressure = [&](int i) -> double {
    return -2. * pts->at(i)["sol"];
  };
  auto uRhsYPressure = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsXPressure = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto vRhsYPressure = [&](int i) -> double {
    return -2. * pts->at(i)["sol"];
  };
  auto findSaveIter = [&presentIter, &saveIter]() -> void {
    saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
    presentIter = int(std::floor(
        (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) / NavierStokesFunction::deltaTime() + 0.5));
    std::cout << "Initial iteration number = " << presentIter << ",\n"
              << "file will be saved every " << saveIter
              << " iterations\n";
  };
  auto copyPointInformation = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      hatuvel.at(i).point::operator=(pts->at(i));
      hatuvel.at(i).findStencil(&(pts->at(i).getElement()), &hatuvel);
      hatvvel.at(i).point::operator=(pts->at(i));
      hatvvel.at(i).findStencil(&(pts->at(i).getElement()), &hatvvel);
      if (isclose(pts->at(i)[0], HALFVALUE) && isclose(pts->at(i)[1], HALFVALUE)) {
        fixedPointIndex = i;
        std::cout << "fixed point index = " << fixedPointIndex << "\n";
      }
    }
  };
  auto assignInitial = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pts->at(i).setMp(UNITVALUE);
      if (pts->at(i).getCondition() == 'D') {
        pts->at(i).setCondition('N');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'N') {
        pts->at(i).setCondition('D');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'd') {
        pts->at(i).setCondition('n');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'n') {
        pts->at(i).setCondition('d');
        pts->at(i)["bdv"] = ZEROVALUE;
      }
      f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
      f.assignBoundaryValue(uvel.at(i), vvel.at(i), presentIter);
    }
    //        f.loadPreviousValue(
    //                "/home/jjhong0608/docker/Navier-Stokes_Result/2D/1.Lid-driven_cavity_flow/Re_7500-2/AGM_Result_380.000000",
    //                "/home/jjhong0608/docker/Navier-Stokes_Result/2D/2.Backward-facing_step_flow/Re_1400/AGM_Result_640.000000",
    //                &puvel, &pvvel, &ppvel);
  };
  auto assignBoundaryValue = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      f.assignBoundaryValue(uvel.at(i), vvel.at(i), presentIter);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(1);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateRepresentationFormula(1);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
      hatuvel.at(i).calculateRepresentationFormula(1);
      hatvvel.at(i).calculateRepresentationFormula(1);
      hatuvel.at(i).calculateRepresentationFormulaPhiPressure('u');
      hatvvel.at(i).calculateRepresentationFormulaPhiPressure('v');
    }
  };
  auto makeMatrixPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateRepresentationFormula(1);
      item->makeDerivatives();
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto calculateDifferentiationVelocityOld = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp, old_vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivatives(pts, pRhsX, pRhsY, pRhsXp, pRhsYp);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->approximateNaNDerivatives(pts);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivativesTwice(pRhsXDiff, pRhsYDiff);
    }
  };
  auto updateRHSVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
    }
  };
  auto updateRHSPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto updateValues = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      if (pts->at(i).getCondition() != 'N') {
        uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
        vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
      } else {
        f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
        uvel.at(i)["sol"] = uvel.at(i)["bdv"];
        vvel.at(i)["sol"] = vvel.at(i)["bdv"];
      }
      hatuvel.at(i).updateRightHandSidePhiPressure(uRhsXPressure, uRhsYPressure, &uvel);
      hatvvel.at(i).updateRightHandSidePhiPressure(vRhsXPressure, vRhsYPressure, &vvel);
      //            pts->at(i)["sol"] -= HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"]);
      uvel.at(i)["dx"] -= pts->at(i)["dxx"] * pointHeat::getDelta();
      vvel.at(i)["dy"] -= pts->at(i)["dyy"] * pointHeat::getDelta();
    }
  };
  auto updateValues1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      ppuvel.at(i) = puvel.at(i);
      ppvvel.at(i) = pvvel.at(i);
      puvel.at(i) = uvel.at(i).getValue();
      pvvel.at(i) = vvel.at(i).getValue();
      if (std::isnan(uvel.at(i)["sol"]) || std::isnan(vvel.at(i)["sol"])) {
        printError("NaN value found.");
      }
    }
  };
  auto subtractPreviousVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] -= puvel.at(i)["sol"];
      vvel.at(i)["sol"] -= pvvel.at(i)["sol"];
      uvel.at(i)["dx"] -= puvel.at(i)["dx"];
      vvel.at(i)["dx"] -= pvvel.at(i)["dx"];
      uvel.at(i)["dy"] -= puvel.at(i)["dy"];
      vvel.at(i)["dy"] -= pvvel.at(i)["dy"];
    }
  };
  auto stopCondition{ZEROVALUE};
  auto checkStopCondition = [&]() -> bool {
    if (presentIter % saveIter != 0) {
      return false;
    }
    auto max_vel{ZEROVALUE}, max_err{ZEROVALUE};
    auto vel{ZEROVALUE}, pvel{ZEROVALUE};
    auto tolerance{1e-6};
    for (int i = 0; i < point::getNPts(); ++i) {
      vel = std::sqrt(
          uvel.at(i)["sol"] * uvel.at(i)["sol"] + vvel.at(i)["sol"] * vvel.at(i)["sol"]);
      pvel = std::sqrt(
          puvel.at(i)["sol"] * puvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"]);
      if (max_vel < vel) {
        max_vel = vel;
      }
      if (max_err < std::abs(vel - pvel)) {
        max_err = std::abs(vel - pvel);
      }
    }
    stopCondition = max_err / max_vel;
    std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "]\n";
    if (stopCondition < tolerance) {
      return true;
    }
    return false;
  };
  auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, pts)};
  auto updateTime = [&]() -> void {
    ++presentIter;
    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
    std::cout << presentIter << "-th iteration, "
              << "current time = [" << pointHeat::getTime() << " / "
              << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
    if (presentIter % saveIter == 0) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(8) << pointHeat::getTime();
      std::string time_string = ss.str();
      wf.writeResult(
          "/home/jhjo/extHD1/tesla_valve/test/AGM_Results_"
          + time_string);
    }
  };
  auto check_neumann = [&]() -> void {
    auto max_grad{0.};
    auto max_grad_idx{0};
    for (int i = 0; i < point::getNPts(); ++i) {
      if (uvel.at(i).getCondition() == 'N') {
        printf(
            "Neumann boundary at (%f, %f) = %23.16e\n"
            "normal = (%23.16e, %23.16e), grad = (%23.16e, %23.16e)\n",
            uvel.at(i).getXy()[0],
            uvel.at(i).getXy()[1],
            uvel.at(i)["dx"] * uvel.at(i).getNormal()[0] + uvel.at(i)["dy"] * uvel.at(i).getNormal()[1],
            uvel.at(i).getNormal()[0],
            uvel.at(i).getNormal()[1],
            uvel.at(i)["dx"],
            uvel.at(i)["dy"]);
        if (max_grad < std::abs(uvel.at(i)["dx"] * uvel.at(i).getNormal()[0] + uvel.at(i)["dy"] * uvel.at(i).getNormal()[1])) {
          max_grad =
              uvel.at(i)["dx"] * uvel.at(i).getNormal()[0] + uvel.at(i)["dy"] * uvel.at(i).getNormal()[1];
          max_grad_idx = i;
        }
      }
    }
    printf(
        "maximum grad = %23.16e, at (%f, %f)\n",
        max_grad,
        pts->at(max_grad_idx).getXy()[0],
        pts->at(max_grad_idx).getXy()[1]);
    exit(1);
  };
  findSaveIter();
  copyPointInformation();
  assignInitial();
  makeMatrixVelocity();
  auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel)};
  //    auto matrixPressure{AGM::matrix<point>(pts)};
  auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPointIndex)};
  auto matrixPhiPressureU{AGM::matrixPhi<pointHeat>(&hatuvel)};
  auto matrixPhiPressureV{AGM::matrixPhi<pointHeat>(&hatvvel)};
  matrixVelocity.makeMatrix();
  matrixVelocity.factorizeMatrix();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocityOld();
  makeMatrixPressure();
  matrixPressure.makeMatrix();
  matrixPressure.factorizeMatrix();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  updateValues();
  matrixPhiPressureU.makeMatrix();
  matrixPhiPressureU.factorizeMatrix();
  matrixPhiPressureU.calculateMatrix();
  matrixPhiPressureV.makeMatrix();
  matrixPhiPressureV.factorizeMatrix();
  matrixPhiPressureV.calculateMatrix();
  updateValues1();
  //    check_neumann();
  while (pointHeat::getTime() + pointHeat::getDelta() < AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
    updateTime();
    assignBoundaryValue();
    updateRHSVelocity();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    subtractPreviousVelocity();
    updateRHSPressure();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    if (checkStopCondition()) {
      break;
    }
    updateValues1();
  }
  updateTime();
  matrixVelocity.releaseMatrix();
  matrixPressure.releaseMatrix();
  wf.writeResult("/home/jhjo/extHD1/tesla_valve/test/AGM_Results_" + std::to_string(pointHeat::getTime()));
}

void AGM::solver::FluidStructureInteraction() {
  auto f{AGM::NavierStokesFunction()};
  int fixedPointIndex{}, presentIter{}, saveIter{};
  point::setNPts(int(pts->size()));
  auto uvel{std::vector<pointHeat>(point::getNPts())};
  auto vvel{std::vector<pointHeat>(point::getNPts())};
  auto puvel{std::vector<value>(point::getNPts())};
  auto pvvel{std::vector<value>(point::getNPts())};
  auto ppvel{std::vector<value>(point::getNPts())};
  auto tildeuvel{std::vector<value>(point::getNPts())};
  auto tildevvel{std::vector<value>(point::getNPts())};
  auto hatuvel{std::vector<value>(point::getNPts())};
  auto hatvvel{std::vector<value>(point::getNPts())};
  pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
  pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
  auto old_uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto old_vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsX = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2 * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsY = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + 2 * puvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsX = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto vRhsY = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + 2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsX = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto pRhsY = [&](int i) -> double {
    return ZEROVALUE;
  };
  auto uRhsX1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto uRhsY1 = [&](int i) -> double {
    return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]);
  };
  auto vRhsX1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto vRhsY1 = [&](int i) -> double {
    return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]);
  };
  auto old_uRhsXp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * puvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dx"];
  };
  auto old_uRhsYp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - uvel.at(i).getMp() * puvel.at(i)["dy"];
  };
  auto old_vRhsXp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dx"];
  };
  auto old_vRhsYp = [&](int i) -> double {
    return 2 * pvvel.at(i)["sol"] * pvvel.at(i)["sol"] - vvel.at(i).getMp() * pvvel.at(i)["dy"];
  };
  auto uRhsXp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * puvel.at(i)["sol"];
  };
  auto uRhsYp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp = [&](int i) -> double {
    return 2 * puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp = [&](int i) -> double {
    return 2 * pvvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto pRhsXp = [&](int i) -> double {
    return uvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto pRhsYp = [&](int i) -> double {
    return vvel.at(i)["sol"] / pointHeat::getDelta();
  };
  auto uRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] - puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto uRhsYp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsXp1 = [&](int i) -> double {
    return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
  };
  auto vRhsYp1 = [&](int i) -> double {
    return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] - pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] + ppvel.at(i)["sol"];
  };
  auto findSaveIter = [&presentIter, &saveIter]() -> void {
    saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
    presentIter = int(std::floor(
        (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) / NavierStokesFunction::deltaTime() + 0.5));
    std::cout << "Initial iteration number = " << presentIter << ",\n"
              << "file will be saved every " << saveIter
              << " iterations\n";
  };
  auto copyPointInformation = [this, &uvel, &vvel, &fixedPointIndex]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).point::operator=(pts->at(i));
      uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
      vvel.at(i).point::operator=(pts->at(i));
      vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
      if (isclose(pts->at(i)[0], HALFVALUE) && isclose(pts->at(i)[1], HALFVALUE)) {
        fixedPointIndex = i;
        std::cout << "fixed point index = " << fixedPointIndex << "\n";
      }
    }
  };
  auto assignInitial = [this, &f, &uvel, &vvel, &puvel, &pvvel, &ppvel, &tildeuvel, &tildevvel]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      pts->at(i).setMp(UNITVALUE);
      if (pts->at(i).getCondition() == 'D') {
        pts->at(i).setCondition('N');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'N') {
        pts->at(i).setCondition('D');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'd') {
        pts->at(i).setCondition('n');
        pts->at(i)["bdv"] = ZEROVALUE;
      } else if (pts->at(i).getCondition() == 'n') {
        pts->at(i).setCondition('d');
        pts->at(i)["bdv"] = ZEROVALUE;
      }
      f.assignPreviousValue(puvel.at(i), pvvel.at(i), ppvel.at(i), uvel.at(i), vvel.at(i), pts->at(i));
      f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
    //        f.loadPreviousValue(
    //                "/home/jjhong0608/docker/Navier-Stokes_Result/2D/1.Lid-driven_cavity_flow/Re_7500-2/AGM_Result_380.000000",
    //                "/home/jjhong0608/docker/Navier-Stokes_Result/2D/2.Backward-facing_step_flow/Re_1400/AGM_Result_640.000000",
    //                &puvel, &pvvel, &ppvel);
  };
  auto assignBoundaryValue = [&f, &uvel, &vvel]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      f.assignBoundaryValue(uvel.at(i), vvel.at(i), 0);
    }
  };
  auto makeMatrixVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateRepresentationFormula(1);
      uvel.at(i).makeDerivatives();
      uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
      uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateRepresentationFormula(1);
      vvel.at(i).makeDerivatives();
      vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
      vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
    }
  };
  auto makeMatrixPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateRepresentationFormula(1);
      item->makeDerivatives();
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto calculateDifferentiationVelocityOld = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, old_uRhsX, old_uRhsY, old_uRhsXp, old_uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, old_vRhsX, old_vRhsY, old_vRhsXp, old_vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX, uRhsY, uRhsXp, uRhsYp);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX, vRhsY, vRhsXp, vRhsYp);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).calculateDerivatives(&uvel, uRhsX1, uRhsY1, uRhsXp1, uRhsYp1);
      vvel.at(i).calculateDerivatives(&vvel, vRhsX1, vRhsY1, vRhsXp1, vRhsYp1);
    }
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).approximateNaNDerivatives(&uvel);
      vvel.at(i).approximateNaNDerivatives(&vvel);
    }
  };
  auto calculateDifferentiationPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->calculateDerivatives(pts, pRhsX, pRhsY, pRhsXp, pRhsYp);
    }
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->approximateNaNDerivatives(pts);
    }
  };
  auto updateRHSVelocity = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
      uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
      vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
      vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
    }
  };
  auto updateRHSVelocity1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i).updateRightHandSide(uRhsX1, uRhsY1);
      uvel.at(i).updateRightHandSidePart(uRhsXp1, uRhsYp1);
      vvel.at(i).updateRightHandSide(vRhsX1, vRhsY1);
      vvel.at(i).updateRightHandSidePart(vRhsXp1, vRhsYp1);
    }
  };
  auto updateRHSPressure = [&]() -> void {
#pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
      item->updateRightHandSide(pRhsX, pRhsY);
      item->updateRightHandSidePart(pRhsXp, pRhsYp);
    }
  };
  auto updateValues = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      tildeuvel.at(i) = uvel.at(i).getValue();
      tildevvel.at(i) = vvel.at(i).getValue();
      if (pts->at(i).getCondition() != 'N') {
        uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
        vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
      }
      pts->at(i)["sol"] -= HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"]);
      hatuvel.at(i) = uvel.at(i).getValue();
      hatvvel.at(i) = vvel.at(i).getValue();
      uvel.at(i)["bdv"] = ZEROVALUE;
      vvel.at(i)["bdv"] = ZEROVALUE;
    }
  };
  auto updateValues1 = [&]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      puvel.at(i) = uvel.at(i).getValue();
      pvvel.at(i) = vvel.at(i).getValue();
      ppvel.at(i) = pts->at(i).getValue();
    }
  };
  auto subtractPreviousVelocity = [&uvel, &vvel, &puvel, &pvvel]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] -= puvel.at(i)["sol"];
      vvel.at(i)["sol"] -= pvvel.at(i)["sol"];
      uvel.at(i)["dx"] -= puvel.at(i)["dx"];
      vvel.at(i)["dx"] -= pvvel.at(i)["dx"];
      uvel.at(i)["dy"] -= puvel.at(i)["dy"];
      vvel.at(i)["dy"] -= pvvel.at(i)["dy"];
    }
  };
  auto addPreviousVelocity = [&uvel, &vvel, &tildeuvel, &tildevvel]() -> void {
#pragma omp parallel for
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["sol"] += tildeuvel.at(i)["sol"];
      vvel.at(i)["sol"] += tildevvel.at(i)["sol"];
      uvel.at(i)["dx"] += tildeuvel.at(i)["dx"];
      vvel.at(i)["dx"] += tildevvel.at(i)["dx"];
      uvel.at(i)["dy"] += tildeuvel.at(i)["dy"];
      vvel.at(i)["dy"] += tildevvel.at(i)["dy"];
    }
  };
  auto stopCondition{ZEROVALUE};
  auto checkStopCondition = [&presentIter, &saveIter, &stopCondition, &uvel, &vvel, &puvel, &pvvel]() -> bool {
    if (presentIter % saveIter != 0) {
      return false;
    }
    auto max_vel{ZEROVALUE}, max_err{ZEROVALUE};
    auto vel{ZEROVALUE}, pvel{ZEROVALUE};
    auto tolerance{1e-6};
    for (int i = 0; i < point::getNPts(); ++i) {
      vel = std::sqrt(uvel.at(i)["sol"] * uvel.at(i)["sol"] + vvel.at(i)["sol"] * vvel.at(i)["sol"]);
      pvel = std::sqrt(puvel.at(i)["sol"] * puvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"]);
      if (max_vel < vel) {
        max_vel = vel;
      }
      if (max_err < std::abs(vel - pvel)) {
        max_err = std::abs(vel - pvel);
      }
    }
    stopCondition = max_err / max_vel;
    std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "]\n";
    if (stopCondition < tolerance) {
      return true;
    }
    return false;
  };
  auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, pts)};
  //    auto updateTime = [&]() -> void {
  //        ++presentIter;
  //        pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
  //        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
  //                  << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
  //        if (presentIter % saveIter == 0) {
  //            for (int i = 0; i < point::getNPts(); ++i) {
  //                vvel.at(i)["dx"] = tildevvel.at(i)["dx"];
  //            }
  //
  //            wf.writeResult(
  //                    "/home/jhjo/extHD1/Navier-Stokes_Result/2D/0.Taylor-Green_vortex/Case5-2/AGM_Result_"
  //                    + std::to_string(pointHeat::getTime()));
  //        }
  //    };

  structure::setNStructures(40);
  structure::setHf(0.02 * 2.5);
  structure::setHs(0.02 * 5.1);
  auto ustruct{std::vector<structure>(structure::getNStructures())};
  auto vstruct{std::vector<structure>(structure::getNStructures())};
  auto upstruct{std::vector<structure>(structure::getNStructures())};
  auto vpstruct{std::vector<structure>(structure::getNStructures())};
  auto wf_struct{writeFileMultiple<structure, structure, point>(&ustruct, &vstruct, pts)};
  auto structArray{std::array<std::vector<structure> *, 2>{&ustruct, &vstruct}};
  structure::setStructures(&structArray);
  auto makeStruct = [&](int nStruct = 360) -> void {
    auto theta{std::vector<double>(nStruct)};
    auto dTheta{2. * M_PI / nStruct};
    auto r{HALFVALUE};
    theta.at(0) = -2. * M_PI;
    for (int i = 1; i < nStruct; ++i) {
      theta.at(i) = theta.at(i - 1) + dTheta;
    }
    for (int i = 0; i < nStruct; ++i) {
      ustruct.at(i).setIdx(i);
      vstruct.at(i).setIdx(i);
      upstruct.at(i).setIdx(i);
      vpstruct.at(i).setIdx(i);
      ustruct.at(i).setXy(coordinate(r * std::cos(theta.at(i)), r * std::sin(theta.at(i))));
      vstruct.at(i).setXy(coordinate(r * std::cos(theta.at(i)), r * std::sin(theta.at(i))));
      upstruct.at(i).setXy(ustruct.at(i).getXy());
      vpstruct.at(i).setXy(vstruct.at(i).getXy());
    }
  };
  auto initializeStructure = [&]() -> void {
    for (int i = 0; i < structure::getNStructures(); ++i) {
      ustruct.at(i)["sol"] = ZEROVALUE;
      ustruct.at(i)["rhs"] = ZEROVALUE;
      vstruct.at(i)["sol"] = ZEROVALUE;
      vstruct.at(i)["rhs"] = ZEROVALUE;
    }
  };
  auto initializeFluidForce = [&]() -> void {
    for (int i = 0; i < point::getNPts(); ++i) {
      uvel.at(i)["rhs"] = ZEROVALUE;
      vvel.at(i)["rhs"] = ZEROVALUE;
    }
  };
  auto copyFluidToStructure = [&]() -> void {
    for (int i = 0; i < structure::getNStructures(); ++i) {
      for (int j = 0; j < point::getNPts(); ++j) {
        ustruct.at(i).copyFluidVelocity(&uvel.at(j));
        vstruct.at(i).copyFluidVelocity(&vvel.at(j));
      }
    }
  };
  auto updateStructureToFluid = [&]() -> void {
    for (int i = 0; i < structure::getNStructures(); ++i) {
      for (int j = 0; j < point::getNPts(); ++j) {
        ustruct.at(i).structureForceUpdateToFluid(&uvel.at(j));
        vstruct.at(i).structureForceUpdateToFluid(&vvel.at(j));
      }
    }
  };
  auto updateStructureForce = [&]() -> void {
    for (int i = 0; i < structure::getNStructures(); ++i) {
      ustruct.at(i).updateForce(&upstruct, 0);
      vstruct.at(i).updateForce(&vpstruct, 1);
    }
  };
  auto structurePositionUpdate = [&]() -> void {
    for (int i = 0; i < structure::getNStructures(); ++i) {
      ustruct.at(i)[0] += ustruct.at(i)["sol"] * pointHeat::getDelta();
      ustruct.at(i)[1] += vstruct.at(i)["sol"] * pointHeat::getDelta();
      vstruct.at(i).setXy(ustruct.at(i).getXy());
    }
  };
  auto FSI = [&]() -> void {
    initializeStructure();
    copyFluidToStructure();
    updateStructureForce();
    structurePositionUpdate();
    initializeFluidForce();
    updateStructureToFluid();
  };

  auto updateTime = [&]() -> void {
    ++presentIter;
    pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
    std::cout << presentIter << "-th iteration, "
              << "current time = [" << pointHeat::getTime() << " / "
              << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
    if (presentIter % saveIter == 0) {
      for (int i = 0; i < point::getNPts(); ++i) {
        vvel.at(i)["dx"] = tildevvel.at(i)["dx"];
      }
      wf.writeResult(
          "/home/jhjo/extHD1/FSI/test0/AGM_Result_"
          + std::to_string(pointHeat::getTime()));
      wf_struct.writeStruct(
          "/home/jhjo/extHD1/FSI/test0/AGM_Struct_"
          + std::to_string(pointHeat::getTime()));
    }
  };

  findSaveIter();
  copyPointInformation();
  assignInitial();
  makeMatrixVelocity();
  auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel)};
  //    auto matrixPressure{AGM::matrix<point>(pts)};
  auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPointIndex)};
  matrixVelocity.makeMatrix();
  matrixVelocity.factorizeMatrix();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocityOld();
  makeMatrixPressure();
  matrixPressure.makeMatrix();
  matrixPressure.factorizeMatrix();
  matrixPressure.calculateMatrix();
  calculateDifferentiationPressure();
  updateValues();
  updateRHSVelocity1();
  matrixVelocity.calculateMatrix();
  calculateDifferentiationVelocity1();
  addPreviousVelocity();
  updateValues1();
  makeStruct(structure::getNStructures());
  while (pointHeat::getTime() + pointHeat::getDelta() < AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
    updateTime();
    assignBoundaryValue();
    updateRHSVelocity();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    subtractPreviousVelocity();
    updateRHSPressure();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    updateRHSVelocity1();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity1();
    addPreviousVelocity();
    if (checkStopCondition()) {
      break;
    }
    updateValues1();
    FSI();
  }
  updateTime();
  matrixVelocity.releaseMatrix();
  matrixPressure.releaseMatrix();
}
