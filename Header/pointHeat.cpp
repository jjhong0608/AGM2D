//
// Created by 조준홍 on 2022/02/13.
//

#include "pointHeat.h"

double AGM::pointHeat::time;
double AGM::pointHeat::delta;

auto AGM::pointHeat::getTime() -> double {
  return time;
}

void AGM::pointHeat::setTime(double d) {
  pointHeat::time = d;
}

auto AGM::pointHeat::getDelta() -> double {
  return delta;
}

void AGM::pointHeat::setDelta(double d) {
  pointHeat::delta = d;
}

void AGM::pointHeat::findStencil(const axialElement *axialElement1, std::vector<AGM::pointHeat> *vector) {
  int n{};
  for (const auto &item : *axialElement1) {
    if (item) {
      element.at(n) = &(vector->at(item->getIdx()));
    }
    ++n;
  }
}

void AGM::pointHeat::calculateRepresentationFormulaCross() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gfuncX{GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta)};
  auto gfuncY{GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta)};
  std::array<matrixRow, 2> row{};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::point::calculateRepresentationFormulaCross", "getEachMp");
          return ZEROVALUE;
        };
        double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                                   : Error();
        return rtv;
      };
      double mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      double mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      double mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      double mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[i][getIdx() + getNPts()] += row[i][pt->getIdx() + getNPts()];
      row[i].remove(pt->getIdx() + getNPts());
    }
  };
  row[0][getIdx()] = -UNITVALUE;
  row[0][element[W]->getIdx()] = mp * gfuncX.green_function_t(xm);
  row[0][element[E]->getIdx()] = -mp * gfuncX.green_function_t(xp);

  row[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral('l', ord);
  row[0][getIdx() + getNPts()] = gfuncX.green_integral('c', ord);
  row[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral('r', ord);

  rhsMatrixRow[0][element[W]->getIdx()] = gfuncX.green_integral('l', ord);
  rhsMatrixRow[0][getIdx()] = gfuncX.green_integral('c', ord);
  rhsMatrixRow[0][element[E]->getIdx()] = gfuncX.green_integral('r', ord);

  partMatrixRow[0][element[W]->getIdx()] = gfuncX.green_integral_t('l', ord);
  partMatrixRow[0][getIdx()] = gfuncX.green_integral_t('c', ord);
  partMatrixRow[0][element[E]->getIdx()] = gfuncX.green_integral_t('r', ord);

  row[1][getIdx()] = -UNITVALUE;
  row[1][element[S]->getIdx()] = mp * gfuncY.green_function_t(ym);
  row[1][element[N]->getIdx()] = -mp * gfuncY.green_function_t(yp);

  row[1][element[S]->getIdx() + getNPts()] = -gfuncY.green_integral('l', ord);
  row[1][getIdx() + getNPts()] = -gfuncY.green_integral('c', ord);
  row[1][element[N]->getIdx() + getNPts()] = -gfuncY.green_integral('r', ord);

  rhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gfuncY.green_integral('l', ord);
  rhsMatrixRow[1][getIdx() + getNPts()] = gfuncY.green_integral('c', ord);
  rhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gfuncY.green_integral('r', ord);

  partMatrixRow[1][element[S]->getIdx() + getNPts()] = gfuncY.green_integral_t('l', ord);
  partMatrixRow[1][getIdx() + getNPts()] = gfuncY.green_integral_t('c', ord);
  partMatrixRow[1][element[N]->getIdx() + getNPts()] = gfuncY.green_integral_t('r', ord);

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);
  eraseInterface(element[N], 1);
  eraseInterface(element[S], 1);

  solMatrixRow[0] = row[0] + row[1];
  solMatrixRow[1] = row[0] - row[1];
}

auto AGM::pointHeat::calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointHeat::calculateRepresentationFormulaNeumannOnAxial", "nullptr");
    return ZEROVALUE;
  };
  point *ptc = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(1) : getAxialLine(axis)->back()->getIdx() == getIdx() ? *std::prev(getAxialLine(axis)->end() - 1)
                                                                                                                                                : nullptr;
  point *ptl = getAxialLine(axis)->front()->getIdx() == getIdx() ? this : getAxialLine(axis)->back()->getIdx() == getIdx() ? *std::prev(getAxialLine(axis)->end() - 2)
                                                                                                                           : nullptr;
  point *ptr = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(2) : getAxialLine(axis)->back()->getIdx() == getIdx() ? this
                                                                                                                                                : nullptr;
  std::string string = getAxialLine(axis)->front()->getIdx() == getIdx() ? "ND" : getAxialLine(axis)->back()->getIdx() == getIdx() ? "DN"
                                                                                                                                   : "";
  double tm = ptl ? ptl->getXy()[axisInt] : Error();
  double tb = ptc ? ptc->getXy()[axisInt] : Error();
  double tp = ptr ? ptr->getXy()[axisInt] : Error();
  double signPhi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
  auto gFunc{GreenfunctionReactionDiffusion(tm, tb, tp, mp, mp, UNITVALUE / delta)};

  if (string == "ND") {
    if (iszero(gFunc.green_function_ND(tm))) {
      return calculateRepresentationFormulaNeumannOffAxial(axis, axisInt);
    }
  } else if (string == "DN") {
    if (iszero(gFunc.green_function_DN(tp))) {
      return calculateRepresentationFormulaNeumannOffAxial(axis, axisInt);
    }
  }

  matrixRow row{};
  if (string == "ND") {
    row[ptl->getIdx()] = -mp * gFunc.green_function_ND(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = -mp * gFunc.green_function_t_ND(tp);

    row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + (axisInt + 4) * getNPts()] =
        gFunc.green_integral_t_ND('l', ord) + gFunc.green_function_ND(tm);
    row[ptc->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_ND('c', ord);
    row[ptr->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_ND('r', ord);
  } else if (string == "DN") {
    row[ptl->getIdx()] = mp * gFunc.green_function_t_DN(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = mp * gFunc.green_function_DN(tp);

    row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_DN('l', ord);
    row[ptc->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_DN('c', ord);
    row[ptr->getIdx() + (axisInt + 4) * getNPts()] =
        gFunc.green_integral_t_DN('r', ord) - gFunc.green_function_DN(tp);
  }
  auto c = -row[getIdx()];
  for (auto &item : row) {
    item.value /= c;
  }
  row.remove(getIdx());
  return row;
}

auto AGM::pointHeat::calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> AGM::matrixRow {
  double tm{}, tb{}, tp{};
  if (axis == 'x') {
    tm = element[W] ? element[W]->getXy()[0] : element[WN]->getXy()[0];
    tb = getXy()[0];
    tp = element[E] ? element[E]->getXy()[0] : element[EN]->getXy()[0];
  } else if (axis == 'y') {
    tm = element[S] ? element[S]->getXy()[1] : element[SE]->getXy()[1];
    tb = getXy()[1];
    tp = element[N] ? element[N]->getXy()[1] : element[NE]->getXy()[1];
  }
  char realAxis{};
  for (const auto &item : {'y', 'x'}) {
    if (getAxialLine(item)) realAxis = item;
  }
  double sign = axis == 'x' ? UNITVALUE : -UNITVALUE;
  auto gFunc{GreenfunctionReactionDiffusion(tm, tb, tp, mp, mp, UNITVALUE / delta)};
  auto approximateSol = [this, &realAxis](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{GreenfunctionReactionDiffusion(m, b, p, d, d, UNITVALUE / delta)};
    auto mRow{matrixRow()};
    double sign = realAxis == 'x' ? UNITVALUE : -UNITVALUE;

    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R', ord);

    mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int i, int plus) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl,
                          double mp0, Greenfunction *func, double d,
                          int i, int i0, char c) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += sign * func->green_integral_tau(c, ord);
      row[pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c, ord);
      row[pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
      row += linearApproximation(ptr, ptl, sign * func->green_integral_tau(c, ord), i, 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), i, i0 + 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), i, i0 + 4);
    }
  };
  row[getIdx() + getNPts()] += sign * gFunc.green_integral_tau('c', ord);
  row[getIdx() + (axisInt + 2) * getNPts()] += gFunc.green_integral_tau('c', ord);
  row[getIdx() + (axisInt + 4) * getNPts()] += gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp;

  if (axis == 'x') {
    assignMatrix(element[E], element[EN], element[ES], -mp, &gFunc, tp, 1, 0, 'r');
    assignMatrix(element[W], element[WN], element[WS], mp, &gFunc, tm, 1, 0, 'l');
  } else if (axis == 'y') {
    assignMatrix(element[N], element[NE], element[NW], -mp, &gFunc, tp, 0, 1, 'r');
    assignMatrix(element[S], element[SE], element[SW], mp, &gFunc, tm, 0, 1, 'l');
  }
  return row;
}

void AGM::pointHeat::calculateRepresentationFormulaInterface() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointHeat::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  double mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  double mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFuncX{GreenfunctionReactionDiffusion(xm, xb, xp, mpw, mpe, UNITVALUE / delta)};
  auto gFuncY{GreenfunctionReactionDiffusion(ym, yb, yp, mps, mpn, UNITVALUE / delta)};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    double mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    double mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    double mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    double mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  bool isInterface = checkInterface(this);
  auto checkMatrixRow = [&checkInterface](matrixRow *row, point *ptr, point *ptl) -> void {
    if (checkInterface(ptr)) {
      (*row)[ptl->getIdx() + getNPts()] += (*row)[ptr->getIdx() + getNPts()];
      row->remove(ptr->getIdx() + getNPts());
    } else if (checkInterface(ptl)) {
      (*row)[ptr->getIdx() + getNPts()] += (*row)[ptl->getIdx() + getNPts()];
      row->remove(ptl->getIdx() + getNPts());
    }
  };
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i,
                                                double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{GreenfunctionReactionDiffusion(m, b, p, d, d, UNITVALUE / delta)};
    auto mRow{matrixRow()};
    double sign = i ? -UNITVALUE : UNITVALUE;
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int i, int plus) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  std::array<matrixRow, 2> row{};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl,
                          double mp0,
                          GreenfunctionReactionDiffusion *func,
                          double d, int i, int i0, char c,
                          char C) -> void {
    double sign = i0 ? -UNITVALUE : UNITVALUE;
    if (pt) {
      row[i0][pt->getIdx()] += mp0 * func->green_function_t(d);
      row[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral(C, ord)
                                                       : sign * func->green_integral(c, ord);
      row[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral(c, ord);
      row[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_t(c, ord);
    } else {
      row[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_t(d), i, std::abs(mp0));
      row[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral(C, ord), i, 1)
                             : linearApproximation(ptr, ptl, sign * func->green_integral(c, ord), i, 1);
      row[i0] += linearApproximation(ptr, ptl, func->green_integral(c, ord), i, i0 + 2);
      row[i0] += linearApproximation(ptr, ptl, func->green_integral_t(c, ord), i, i0 + 4);
    }
  };
  row[0][getIdx()] = -UNITVALUE;
  if (!isInterface) row[0][getIdx() + getNPts()] = gFuncX.green_integral('c', ord);
  row[0][getIdx() + 2 * getNPts()] = gFuncX.green_integral('c', ord);
  row[0][getIdx() + 4 * getNPts()] = gFuncX.green_integral_t('c', ord);
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

  row[1][getIdx()] = -UNITVALUE;
  if (!isInterface) row[1][getIdx() + getNPts()] = -gFuncY.green_integral('c', ord);
  row[1][getIdx() + 3 * getNPts()] = gFuncY.green_integral('c', ord);
  row[1][getIdx() + 5 * getNPts()] = gFuncY.green_integral_t('c', ord);
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFuncY, yp, 0, 1, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFuncY, ym, 0, 1, 'l', 'L');

  for (int i = 0; i < 2; ++i) {
    if (!row[i].empty()) {
      while (row[i].back().idx >= 4 * getNPts()) {
        partMatrixRow[i][row[i].back().idx - 4 * getNPts()] = row[i].back().value;
        row[i].pop_back();
      }
      while (row[i].back().idx >= 2 * getNPts()) {
        rhsMatrixRow[i][row[i].back().idx - 2 * getNPts()] = row[i].back().value;
        row[i].pop_back();
      }
    }
  }
  solMatrixRow[0] = row[0] + row[1];
  if (isInterface) {
    solMatrixRow[1][getIdx() + getNPts()] = UNITVALUE;
  } else {
    solMatrixRow[1] = row[0] - row[1];
  }
}

void AGM::pointHeat::makePhiCoefficient(std::vector<pointHeat> *vector) {
  for (auto item : solMatrixRow[1]) {
    if (item.idx < getNPts()) {
      rb[1] -= item.value * vector->at(item.idx).getValue()["sol"];
    }
  }
}

void AGM::pointHeat::updateRightHandSidePhiPressure(
    const std::function<double(int)> &f,
    const std::function<double(int)> &g,
    std::vector<pointHeat> *points) {
  switch (condition) {
    case 'C':
      updateRightHandSidePhiPressureCross(f, g, points);
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      updateRightHandSidePhiPressureDirichlet(f, g, points);
    case 'I':
      updateRightHandSidePhiPressureInterface(f, g, points);
  }
}

void AGM::pointHeat::updateRightHandSidePhiPressureCross(
    const std::function<double(int)> &f,
    const std::function<double(int)> &g,
    std::vector<pointHeat> *points) {
  rb[0] = rb[1] = ZEROVALUE;
  for (const auto &item : solMatrixRow[0]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * points->at(item.idx)["sol"];
      rb[1] -= item.value * points->at(item.idx)["sol"];
    }
  }
  for (const auto &item : solMatrixRow[1]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * points->at(item.idx)["sol"];
      rb[1] += item.value * points->at(item.idx)["sol"];
    }
  }
  for (const auto &item : partMatrixRow[0]) {
    rb[0] -= item.value * f(item.idx);
    rb[1] -= item.value * f(item.idx);
  }
  for (const auto &item : partMatrixRow[1]) {
    rb[0] -= item.value * g(item.idx - getNPts());
    rb[1] += item.value * g(item.idx - getNPts());
  }
}

void AGM::pointHeat::updateRightHandSidePhiPressureDirichlet(
    const std::function<double(int)> &f,
    const std::function<double(int)> &g,
    std::vector<pointHeat> *points) {
}

void AGM::pointHeat::updateRightHandSidePhiPressureInterface(
    const std::function<double(int)> &f,
    const std::function<double(int)> &g,
    std::vector<pointHeat> *points) {
  bool isinterface{solMatrixRow[1].size() == 1};
  rb[0] = rb[1] = ZEROVALUE;
  for (const auto &item : solMatrixRow[0]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * points->at(item.idx)["sol"];
      rb[1] -= item.value * points->at(item.idx)["sol"];
    }
  }
  for (const auto &item : solMatrixRow[1]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * points->at(item.idx)["sol"];
      rb[1] += item.value * points->at(item.idx)["sol"];
    }
  }
  for (const auto &item : partMatrixRow[0]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * f(item.idx);
      if (!isinterface) rb[1] -= item.value * f(item.idx);
    } else if (item.idx < 2 * getNPts()) {
      rb[0] -= item.value * g(item.idx - getNPts());
      if (!isinterface) rb[1] -= item.value * g(item.idx - getNPts());
    } else {
      printError("AGM::point::updateRightHandSidePhiPressureInterface");
    }
  }
  for (const auto &item : partMatrixRow[1]) {
    if (item.idx < getNPts()) {
      rb[0] -= item.value * f(item.idx);
      if (!isinterface) rb[1] += item.value * f(item.idx);
    } else if (item.idx < 2 * getNPts()) {
      rb[0] -= item.value * g(item.idx - getNPts());
      if (!isinterface) rb[1] += item.value * g(item.idx - getNPts());
    } else {
      printError("AGM::point::updateRightHandSidePhiPressureInterface");
    }
  }
}

void AGM::pointHeat::makeDerivativesCross() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gfuncX{GreenfunctionReactionDiffusion(xm, xb, xp, mp, mp, UNITVALUE / delta)};
  auto gfuncY{GreenfunctionReactionDiffusion(ym, yb, yp, mp, mp, UNITVALUE / delta)};
  deriMatrixRow[0][element[W]->getIdx()] = mp * gfuncX.green_function_ttau(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -mp * gfuncX.green_function_ttau(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral_tau('l', ord);
  deriMatrixRow[0][getIdx() + getNPts()] = gfuncX.green_integral_tau('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral_tau('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('l', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('l', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('c', ord) + UNITVALUE / mp;
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('r', ord);

  deriMatrixRow[1][element[S]->getIdx()] = mp * gfuncY.green_function_ttau(ym);
  deriMatrixRow[1][element[N]->getIdx()] = -mp * gfuncY.green_function_ttau(yp);

  deriMatrixRow[1][element[S]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('l', ord);
  deriMatrixRow[1][getIdx() + getNPts()] = -gfuncY.green_integral_tau('c', ord);
  deriMatrixRow[1][element[N]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('r', ord);

  deriMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('l', ord);
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('c', ord);
  deriMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('r', ord);

  deriMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('l', ord);
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('c', ord) + UNITVALUE / mp;
  deriMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('r', ord);
}

void AGM::pointHeat::makeDerivativesInterface() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::point::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  double mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  double mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFuncX{GreenfunctionReactionDiffusion(xm, xb, xp, mpw, mpe, UNITVALUE / delta)};
  auto gFuncY{GreenfunctionReactionDiffusion(ym, yb, yp, mps, mpn, UNITVALUE / delta)};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    double mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    double mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    double mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    double mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  bool isInterface = checkInterface(this);
  auto checkMatrixRow = [&checkInterface](matrixRow *row, point *ptr, point *ptl) -> void {
    if (checkInterface(ptr)) {
      (*row)[ptl->getIdx() + getNPts()] += (*row)[ptr->getIdx() + getNPts()];
      row->remove(ptr->getIdx() + getNPts());
    } else if (checkInterface(ptl)) {
      (*row)[ptr->getIdx() + getNPts()] += (*row)[ptl->getIdx() + getNPts()];
      row->remove(ptl->getIdx() + getNPts());
    }
  };
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i,
                                                double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{GreenfunctionReactionDiffusion(m, b, p, d, d, UNITVALUE / delta)};
    auto mRow{matrixRow()};
    double sign = i ? -UNITVALUE : UNITVALUE;
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int i, int plus) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [this, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl,
                                                                                  double mp0,
                                                                                  GreenfunctionReactionDiffusion *func,
                                                                                  double d, int i, int i0, char c,
                                                                                  char C) -> void {
    double sign = i0 ? -UNITVALUE : UNITVALUE;
    if (pt) {
      deriMatrixRow[i0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      deriMatrixRow[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral_tau(C, ord)
                                                                 : sign * func->green_integral_tau(c, ord);
      deriMatrixRow[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c, ord);
      deriMatrixRow[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      deriMatrixRow[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
      deriMatrixRow[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral_tau(C, ord), i,
                                                             1)
                                       : linearApproximation(ptr, ptl, sign * func->green_integral_tau(c, ord), i,
                                                             1);
      deriMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), i, i0 + 2);
      deriMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), i, i0 + 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFuncX.green_integral_tau('c', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFuncX.green_integral_tau('c', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFuncX.green_integral_ttau('c', ord) + UNITVALUE / mp;
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

  if (!isInterface) deriMatrixRow[1][getIdx() + getNPts()] = -gFuncY.green_integral_tau('c', ord);
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFuncY.green_integral_tau('c', ord);
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFuncY.green_integral_ttau('c', ord) + UNITVALUE / mp;
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFuncY, yp, 0, 1, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFuncY, ym, 0, 1, 'l', 'L');
}

void AGM::pointHeat::calculateDerivatives(
    const std::vector<pointHeat> *points,
    const std::function<double(int)> &f,
    const std::function<double(int)> &g,
    const std::function<double(int)> &fp,
    const std::function<double(int)> &gp) {
  auto assignDerivatives = [&](int i) -> double {
    auto d{ZEROVALUE};
    for (const auto &item : deriMatrixRow[i]) {
      if (item.idx < getNPts()) {
        d += item.value * points->at(item.idx)["sol"];
      } else if (item.idx < 2 * getNPts()) {
        d += item.value * points->at(item.idx - getNPts())["phi"];
      } else if (item.idx < 3 * getNPts()) {
        d += item.value * f(item.idx - 2 * getNPts());
      } else if (item.idx < 4 * getNPts()) {
        d += item.value * g(item.idx - 3 * getNPts());
      } else if (item.idx < 5 * getNPts()) {
        d += item.value * fp(item.idx - 4 * getNPts());
      } else if (item.idx < 6 * getNPts()) {
        d += item.value * gp(item.idx - 5 * getNPts());
      } else {
        printError(
            "AGM::pointHeat::calculateDerivatives",
            "item.idx (which is %d) is too large",
            item.idx);
      }
    }
    return d;
  };
  values["dx"] = assignDerivatives(0);
  values["dy"] = assignDerivatives(1);
}

void AGM::pointHeat::approximateNaNDerivatives(std::vector<AGM::pointHeat> *points) {
  auto findInnerPointOfBoundary = [this]() -> point * {
    if (getCondition() == 'd' || getCondition() == 'n') {
      for (const auto &item : {E, W, N, S}) {
        if (getElement()[item] && getIdx() != getElement()[item]->getIdx()) {
          return getElement()[item];
        }
      }
    }

    for (const auto &item : {'x', 'y'}) {
      if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
        return getAxialLine(item)->at(1);
      }
      if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
        return *std::prev(getAxialLine(item)->end() - 1);
      }
    }
    printf("condition = %c\n", getCondition());
    printInformation();
    printError("AGM::pointHeat::approximateNaNDerivatives", "findInnerPointOfBoundary");
    return nullptr;
  };
  if (std::isnan(values["dx"])) values["dx"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dx"];
  if (std::isnan(values["dy"])) values["dy"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dy"];
}

void AGM::pointHeat::approximateDiff(std::vector<pointHeat> *points) {
  for (auto &axis : {'x', 'y'}) {
    if (getAxialLine(axis)) {
      if (getAxialLine(axis)->front()->getIdx() == getIdx()) {
        values["dx"] = points->at(getAxialLine(axis)->at(1)->getIdx()).getValue()["dx"];
        values["dy"] = points->at(getAxialLine(axis)->at(1)->getIdx()).getValue()["dy"];
        return;
      } else if (getAxialLine(axis)->back()->getIdx() == getIdx()) {
        values["dx"] = points->at((*std::prev(getAxialLine(axis)->end() - 1))->getIdx()).getValue()["dx"];
        values["dy"] = points->at((*std::prev(getAxialLine(axis)->end() - 1))->getIdx()).getValue()["dy"];
        return;
      }
    }
  }
  for (const auto &item : {E, W, N, S}) {
    if (getElement()[item] && getIdx() != getElement()[item]->getIdx()) {
      values["dx"] = points->at(getElement()[item]->getIdx())["dx"];
      values["dy"] = points->at(getElement()[item]->getIdx())["dy"];
      return;
    }
  }
}
