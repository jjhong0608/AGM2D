//
// Created by 조준홍 on 2022/02/19.
//

#include "solver.h"

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
        f.assignBoundaryValue(*item);
    }
    point::setNPts(int(pts->size()));

    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        item->calculateRepresentationFormula();
        item->makeDerivatives();
        item->updateRightHandSide(rhsX, rhsY);
        item->updateRightHandSidePart(rhsXp, rhsXp);
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
    wf.writeResult("/home/jjhong0608/docker/AGM2D/Axisymmetric/AGM_Result");
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
            "/home/jjhong0608/docker/AGM2D/New-nonlinear_algorithm/BFS_flow/Re800-4/AGM_Result_250.000000")};
    if (NS.fail()) {
        printError("file is not opened");
    }
    std::cout << "file is opened\n";
    int idx{}, bc{}, bdnum{};
    double x{}, y{}, p{};
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
        NS >> idx >> x >> y >> uvel.at(idx)["sol"] >> vvel.at(idx)["sol"] >> p >> uvel.at(idx)["dx"]
           >> uvel.at(idx)["dy"] >> vvel.at(idx)["dx"] >> vvel.at(idx)["dy"] >> uvel.at(idx)["phi"]
           >> vvel.at(idx)["phi"] >> bc;
        if (isclose(x, 1.5e1)) {
            rightbound.emplace_back(uvel.at(idx)["sol"]);
            righty.emplace_back(y);
        }
        if (pt.getCondition() == 'D') {
            pt["bdv"] = f.u(pt);
        } else if (pt.getCondition() == 'N') {
            pt.setCondition('D');
            pt["bdv"] = assignBoundary(++bdnum);
//            pt["bdv"] = f.u(pt);
        }
        pt.setMp(UNITVALUE);
    });

    NS.close();

    #pragma omp parallel for
    for (auto item = pts->begin(); item != pts->end(); ++item) {
        item->calculateRepresentationFormula();
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
            "/home/jjhong0608/docker/AGM2D/New-nonlinear_algorithm/BFS_flow/Re800-4/AGM_Result_stream");
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
            item->calculateRepresentationFormula();
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
        for (const auto &row: item->getSolMatrixRow()[0]) {
            if (std::isnan(row.value)) {
                std::cout << "(x, y) = (" << item->getXy()[0] << ", " << item->getXy()[1] << ")\n";
                std::cout << "condition = " << item->getCondition() << "\n";
                printError("find NaN0 value");
            }
        }
        for (const auto &row: item->getSolMatrixRow()[1]) {
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
        return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) + previousValues.at(i)["phi"] +
               previousValues.at(i)["sol"] / pointHeat::getDelta();
    };
    auto rhsY = [&](int i) -> double {
        return HALFVALUE * (ptsHeat.at(i)["rhs"] + previousValues.at(i)["rhs"]) - previousValues.at(i)["phi"] +
               previousValues.at(i)["sol"] / pointHeat::getDelta();
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
            item->calculateRepresentationFormula();
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
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::heatFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
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

void AGM::solver::NavierStokesSolver() {
    auto f{AGM::NavierStokesFunction()};
    int fixedPointIndex{}, presentIter{}, saveIter{};
    point::setNPts(int(pts->size()));
    auto uvel{std::vector<pointHeat>(point::getNPts())};
    auto vvel{std::vector<pointHeat>(point::getNPts())};
    auto puvel{std::vector<value>(point::getNPts())};
    auto pvvel{std::vector<value>(point::getNPts())};
    auto ppvel{std::vector<value>(point::getNPts())};
    auto uvalue{std::vector<value>(point::getNPts())};
    auto vvalue{std::vector<value>(point::getNPts())};
    pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
    pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
    auto uRhsX = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) + puvel.at(i)["phi"] +
               2 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) - puvel.at(i)["phi"] +
               2 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) + pvvel.at(i)["phi"] +
               2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) - pvvel.at(i)["phi"] +
               2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto pRhsX = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto pRhsY = [&](int i) -> double {
        return ZEROVALUE;
    };
    auto uRhsX1 = [&](int i) -> double {
        return uRhsX(i);
    };
    auto uRhsY1 = [&](int i) -> double {
        return uRhsY(i);
    };
    auto vRhsX1 = [&](int i) -> double {
        return vRhsX(i);
    };
    auto vRhsY1 = [&](int i) -> double {
        return vRhsY(i);
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
        return uvalue.at(i)["sol"] * uvalue.at(i)["sol"] + puvel.at(i)["sol"] * puvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"];
    };
    auto uRhsYp1 = [&](int i) -> double {
        return uvalue.at(i)["sol"] * vvalue.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsXp1 = [&](int i) -> double {
        return uvalue.at(i)["sol"] * vvalue.at(i)["sol"] + puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsYp1 = [&](int i) -> double {
        return vvalue.at(i)["sol"] * vvalue.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"] + pts->at(i)["sol"] +
               ppvel.at(i)["sol"];
    };
    auto pRhsX1 = [&](int i) -> double {
        return -uvel.at(i)["dx"] / pointHeat::getDelta();
    };
    auto pRhsY1 = [&](int i) -> double {
        return -vvel.at(i)["dy"] / pointHeat::getDelta();
    };
    auto findSaveIter = [&presentIter, &saveIter]() -> void {
        saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
        presentIter = int(std::floor(
                (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) /
                NavierStokesFunction::deltaTime() + 0.5));
        std::cout << "Initial iteration number = " << presentIter << "\n";
    };
    auto copyPointInformation = [this, &uvel, &vvel, &fixedPointIndex]() -> void {
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).point::operator=(pts->at(i));
            uvel.at(i).findStencil(&(pts->at(i).getElement()), &uvel);
            vvel.at(i).point::operator=(pts->at(i));
            vvel.at(i).findStencil(&(pts->at(i).getElement()), &vvel);
            if (isclose(pts->at(i)[0], 0.5) && isclose(pts->at(i)[1], 0.5)) {
                fixedPointIndex = i;
                std::cout << "fixed point index = " << fixedPointIndex << "\n";
            }
        }
    };
    auto assignInitial = [this, &f, &uvel, &vvel, &puvel, &pvvel, &ppvel]() -> void {
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
            f.assignBoundaryValue(uvel.at(i), vvel.at(i));
        }
//        f.loadPreviousValue(
//                "/home/jjhong0608/docker/AGM2D/New-nonlinear_algorithm/BFS_flow/Re800-3/AGM_Result_3.000000", &puvel,
//                &pvvel, &ppvel);
    };
    auto assignBoundaryValue = [&f, &uvel, &vvel]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            f.assignBoundaryValue(uvel.at(i), vvel.at(i));
        }
    };
    auto makeMatrixVelocity = [&]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            uvel.at(i).calculateRepresentationFormula();
            uvel.at(i).makeDerivatives();
            uvel.at(i).updateRightHandSide(uRhsX, uRhsY);
            uvel.at(i).updateRightHandSidePart(uRhsXp, uRhsYp);
            vvel.at(i).calculateRepresentationFormula();
            vvel.at(i).makeDerivatives();
            vvel.at(i).updateRightHandSide(vRhsX, vRhsY);
            vvel.at(i).updateRightHandSidePart(vRhsXp, vRhsYp);
        }
    };
    auto makeMatrixPressure = [&]() -> void {
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->calculateRepresentationFormula();
            item->makeDerivatives();
            item->updateRightHandSide(pRhsX, pRhsY);
            item->updateRightHandSidePart(pRhsXp, pRhsYp);
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
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->calculateDerivativesTwice(pRhsX1, pRhsY1);
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
    auto updateRHSPressrue = [&]() -> void {
        #pragma omp parallel for
        for (auto item = pts->begin(); item != pts->end(); ++item) {
            item->updateRightHandSide(pRhsX, pRhsY);
            item->updateRightHandSidePart(pRhsXp, pRhsYp);
        }
    };
    auto updateValues = [this, &uvel, &vvel, &uvalue, &vvalue]() -> void {
        #pragma omp parallel for
        for (int i = 0; i < point::getNPts(); ++i) {
            if (pts->at(i).getCondition() != 'N') {
                uvel.at(i)["sol"] -= pts->at(i)["dx"] * pointHeat::getDelta();
                vvel.at(i)["sol"] -= pts->at(i)["dy"] * pointHeat::getDelta();
            }
            pts->at(i)["sol"] -= HALFVALUE * uvel.at(i).getMp() * (uvel.at(i)["dx"] + vvel.at(i)["dy"]);
            uvel.at(i)["dx"] -= pts->at(i)["dxx"] * pointHeat::getDelta();
            vvel.at(i)["dy"] -= pts->at(i)["dyy"] * pointHeat::getDelta();
            uvalue.at(i) = uvel.at(i).getValue();
            vvalue.at(i) = vvel.at(i).getValue();
        }
    };
    auto updateValues1 = [this, &uvel, &vvel, &puvel, &pvvel, &ppvel]() -> void {
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
    auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, pts)};
    auto updateTime = [&presentIter, &saveIter, &wf]() -> void {
        ++presentIter;
        pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
                  << AGM::NavierStokesFunction::terminalTime() << "]\n";
//        if (presentIter % saveIter == 0) {
//            wf.writeResult(
//                    "/home/jjhong0608/docker/AGM2D/New-nonlinear_algorithm/BFS_flow/Re800-4/AGM_Result_" +
//                    std::to_string(pointHeat::getTime()));
//        }
    };
    findSaveIter();
    copyPointInformation();
    assignInitial();
    makeMatrixVelocity();
    auto matrixVelocity{AGM::matrixMulti<pointHeat>(&uvel, &vvel)};
    auto matrixPressure{AGM::matrix<point>(pts)};
//    auto matrixPressure{AGM::matrixNormal<point>(pts, fixedPoin`tIndex)};
    matrixVelocity.makeMatrix();
    matrixVelocity.factorizeMatrix();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity();
    subtractPreviousVelocity();
    makeMatrixPressure();
    matrixPressure.makeMatrix();
    matrixPressure.factorizeMatrix();
    matrixPressure.calculateMatrix();
    calculateDifferentiationPressure();
    updateValues();
    updateRHSVelocity1();
    matrixVelocity.calculateMatrix();
    calculateDifferentiationVelocity1();
    subtractPreviousVelocity();
    updateValues1();
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
        updateTime();
        assignBoundaryValue();
        updateRHSVelocity();
        matrixVelocity.calculateMatrix();
        calculateDifferentiationVelocity();
        subtractPreviousVelocity();
        updateRHSPressrue();
        matrixPressure.calculateMatrix();
        calculateDifferentiationPressure();
        updateValues();
        updateRHSVelocity1();
        matrixVelocity.calculateMatrix();
        calculateDifferentiationVelocity1();
        subtractPreviousVelocity();
        updateValues1();
    }
    updateTime();
    matrixVelocity.releaseMatrix();
    matrixPressure.releaseMatrix();

//    wf.writeResult("/home/jjhong0608/docker/AGM2D/New-nonlinear_algorithm/BFS_flow/Re800-4/AGM_Result");
}
