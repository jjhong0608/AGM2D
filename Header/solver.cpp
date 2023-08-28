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
            "/home/jhjo/extHD1/Navier-Stokes_Result/2D/1.Lid-driven_cavity_flow_test/Re_5000/AGM_Result_266.000000"
    )};
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
//        if (isclose(x, 6e1)) {
//            rightbound.emplace_back(uvel.at(idx)["sol"]);
//            righty.emplace_back(y);
//        }
        if (pt.getCondition() == 'D') {
            pt["bdv"] = f.u(pt);
        } else if (pt.getCondition() == 'N') {
            pt.setCondition('D');
            pt["bdv"] = assignBoundary(++bdnum);
//            pt["bdv"] = f.u(pt);
        }
//        pt["bdv"] = ZEROVALUE;
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
            "/home/jhjo/extHD1/Navier-Stokes_Result/2D/1.Lid-driven_cavity_flow_test/Re_5000/AGM_Result_stream"
    );
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
    auto tildeuvel{std::vector<value>(point::getNPts())};
    auto tildevvel{std::vector<value>(point::getNPts())};
    auto hatuvel{std::vector<value>(point::getNPts())};
    auto hatvvel{std::vector<value>(point::getNPts())};
    auto phatuvel{std::vector<value>(point::getNPts())};
    auto phatvvel{std::vector<value>(point::getNPts())};
    pointHeat::setTime(AGM::NavierStokesFunction::initialTime());
    pointHeat::setDelta(AGM::NavierStokesFunction::deltaTime());
    auto old_uRhsX = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsX = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2. * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2. * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
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
        return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] - phatuvel.at(i)["sol"] * phatuvel.at(i)["sol"] +
               pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto uRhsYp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - phatuvel.at(i)["sol"] * phatvvel.at(i)["sol"];
    };
    auto vRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - phatuvel.at(i)["sol"] * phatvvel.at(i)["sol"];
    };
    auto vRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] - phatvvel.at(i)["sol"] * phatvvel.at(i)["sol"] +
               pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto findSaveIter = [&]() -> void {
        saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
        presentIter = int(std::floor(
                (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) /
                NavierStokesFunction::deltaTime() + 0.5));
        std::cout << "Initial iteration number = " << presentIter << ",\n" << "file will be saved every " << saveIter
                  << " iterations\n";
    };
    auto copyPointInformation = [&]() -> void {
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
            uvel.at(i).calculateRepresentationFormula();
            uvel.at(i).makeDerivatives();
            uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
            uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
            vvel.at(i).calculateRepresentationFormula();
            vvel.at(i).makeDerivatives();
            vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
            vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
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
            phatuvel.at(i) = hatuvel.at(i);
            phatvvel.at(i) = hatvvel.at(i);
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
    auto addPreviousVelocity = [&uvel, &vvel, &tildeuvel, &tildevvel, &hatuvel, &hatvvel]() -> void {
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
        stopCondition = max_err / max_vel;
        std::cout << "Stop Condition: [" << stopCondition << " / " << tolerance << "], Div = "
                  << mean_div / size << "\n";
        if (stopCondition < tolerance) {
            return true;
        }
        return false;
    };
    auto wf{writeFileMultiple<pointHeat, pointHeat, point>(&uvel, &vvel, pts)};
    auto updateTime = [&]() -> void {
        ++presentIter;
        pointHeat::setTime(pointHeat::getTime() + pointHeat::getDelta());
        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
                  << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
        if (presentIter % saveIter == 0) {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(8) << pointHeat::getTime();
            std::string time_string = ss.str();
            wf.writeResult(
                    "/home/jhjo/extHD1/tesla_valve/test/AGM_Results_"
                    + time_string
            );
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
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
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
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsX = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2. * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2. * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2. * pvvel.at(i)["sol"] / pointHeat::getDelta();
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
        return ZEROVALUE;
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
                (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) /
                NavierStokesFunction::deltaTime() + 0.5));
        std::cout << "Initial iteration number = " << presentIter << ",\n" << "file will be saved every " << saveIter
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
            uvel.at(i).calculateRepresentationFormula();
            uvel.at(i).makeDerivatives();
            uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
            uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
            vvel.at(i).calculateRepresentationFormula();
            vvel.at(i).makeDerivatives();
            vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
            vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
            hatuvel.at(i).calculateRepresentationFormula();
            hatvvel.at(i).calculateRepresentationFormula();
            hatuvel.at(i).calculateRepresentationFormulaPhiPressure('u');
            hatvvel.at(i).calculateRepresentationFormulaPhiPressure('v');
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
                    uvel.at(i)["sol"] * uvel.at(i)["sol"] + vvel.at(i)["sol"] * vvel.at(i)["sol"]
            );
            pvel = std::sqrt(
                    puvel.at(i)["sol"] * puvel.at(i)["sol"] + pvvel.at(i)["sol"] * pvvel.at(i)["sol"]
            );
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
        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
                  << AGM::NavierStokesFunction::terminalTime() << "], Stopping Error = [" << stopCondition << "]\n";
        if (presentIter % saveIter == 0) {
            std::stringstream ss;
            ss << std::fixed << std::setprecision(8) << pointHeat::getTime();
            std::string time_string = ss.str();
            wf.writeResult(
                    "/home/jhjo/extHD1/tesla_valve/test/AGM_Results_"
                    + time_string
            );
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
                        uvel.at(i)["dy"]
                );
                if (max_grad <
                    std::abs(uvel.at(i)["dx"] * uvel.at(i).getNormal()[0] +
                             uvel.at(i)["dy"] * uvel.at(i).getNormal()[1])) {
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
                pts->at(max_grad_idx).getXy()[1]
        );
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
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
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
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto old_vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsX = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto uRhsY = [&](int i) -> double {
        return HALFVALUE * (uvel.at(i)["rhs"] + puvel.at(i)["rhs"]) +
               2 * puvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsX = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
    };
    auto vRhsY = [&](int i) -> double {
        return HALFVALUE * (vvel.at(i)["rhs"] + pvvel.at(i)["rhs"]) +
               2 * pvvel.at(i)["sol"] / pointHeat::getDelta();
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
        return hatuvel.at(i)["sol"] * hatuvel.at(i)["sol"] - puvel.at(i)["sol"] * puvel.at(i)["sol"] +
               pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto uRhsYp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsXp1 = [&](int i) -> double {
        return hatuvel.at(i)["sol"] * hatvvel.at(i)["sol"] - puvel.at(i)["sol"] * pvvel.at(i)["sol"];
    };
    auto vRhsYp1 = [&](int i) -> double {
        return hatvvel.at(i)["sol"] * hatvvel.at(i)["sol"] - pvvel.at(i)["sol"] * pvvel.at(i)["sol"] +
               pts->at(i)["sol"] + ppvel.at(i)["sol"];
    };
    auto findSaveIter = [&presentIter, &saveIter]() -> void {
        saveIter = int(std::floor(NavierStokesFunction::writeTime() / NavierStokesFunction::deltaTime() + 0.5));
        presentIter = int(std::floor(
                (std::fmod(NavierStokesFunction::initialTime(), NavierStokesFunction::writeTime())) /
                NavierStokesFunction::deltaTime() + 0.5));
        std::cout << "Initial iteration number = " << presentIter << ",\n" << "file will be saved every " << saveIter
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
            uvel.at(i).calculateRepresentationFormula();
            uvel.at(i).makeDerivatives();
            uvel.at(i).updateRightHandSide(old_uRhsX, old_uRhsY);
            uvel.at(i).updateRightHandSidePart(old_uRhsXp, old_uRhsYp);
            vvel.at(i).calculateRepresentationFormula();
            vvel.at(i).makeDerivatives();
            vvel.at(i).updateRightHandSide(old_vRhsX, old_vRhsY);
            vvel.at(i).updateRightHandSidePart(old_vRhsXp, old_vRhsYp);
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
        std::cout << presentIter << "-th iteration, " << "current time = [" << pointHeat::getTime() << " / "
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
    while (pointHeat::getTime() + pointHeat::getDelta() <
           AGM::NavierStokesFunction::terminalTime() - HALFVALUE * pointHeat::getDelta()) {
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
