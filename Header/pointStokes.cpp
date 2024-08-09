//
// Created by 조준홍 on 10/27/23.
//

#include "pointStokes.h"

void AGM::pointStokes::findStencil(const AGM::axialElement *axialElement1, std::vector<pointStokes> *vector) {
  int n{};
  for (const auto &item : *axialElement1) {
    if (item) {
      element.at(n) = &(vector->at(item->getIdx()));
    }
    ++n;
  }
}

void AGM::pointStokes::makeRepresentationFormulaStokes() {
  switch (condition) {
    case 'C':
      makeRepresentationFormulaStokesCross();
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      makeRepresentationFormulaStokesBoundary();
      break;
    case 'I':
      makeRepresentationFormulaStokesInterface();
      break;
    default:
      printError("AGM::pointStokes::makeRepresentationFormulaStokes", "condition (which is %c) is wrong", condition);
  }
}

void AGM::pointStokes::makeRepresentationFormulaStokesCross() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gfuncX{Greenfunction(xm, xb, xp, getMp(), getMp())};
  auto gfuncY{Greenfunction(ym, yb, yp, getMp(), getMp())};
  pMatrixRow[0][element[W]->getIdx()] = getMp() * gfuncX.green_function_ttau(xm);
  pMatrixRow[0][element[E]->getIdx()] = -getMp() * gfuncX.green_function_ttau(xp);

  pMatrixRow[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + getNPts()] = gfuncX.green_integral_tau('c', getOrd());
  pMatrixRow[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral_tau('r', getOrd());

  pMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('c', getOrd());
  pMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('r', getOrd());

  pMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('l', getOrd());
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('c', getOrd());
  pMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('r', getOrd());

  pMatrixRow[1][element[S]->getIdx()] = getMp() * gfuncY.green_function_ttau(ym);
  pMatrixRow[1][element[N]->getIdx()] = -getMp() * gfuncY.green_function_ttau(yp);

  pMatrixRow[1][element[S]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + getNPts()] = -gfuncY.green_integral_tau('c', getOrd());
  pMatrixRow[1][element[N]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('r', getOrd());

  pMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('c', getOrd());
  pMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('r', getOrd());

  pMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('l', getOrd());
  pMatrixRow[1][getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('c', getOrd());
  pMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('r', getOrd());
}

auto AGM::pointStokes::calculateRepresentationFormulaStokesBoundaryAxial(char axis, int axisInt) -> AGM::matrixRow {
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
  auto gFunc{Greenfunction(tm, tb, tp, mp, mp)};
  auto approximateSol = [this, &realAxis](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{Greenfunction(m, b, p, d, d)};
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
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, int i, int i0, char c) -> void {
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
  row[getIdx() + (axisInt + 4) * getNPts()] += gFunc.green_integral_ttau('c', ord);

  if (axis == 'x') {
    assignMatrix(element[E], element[EN], element[ES], -mp, &gFunc, tp, 1, 0, 'r');
    assignMatrix(element[W], element[WN], element[WS], mp, &gFunc, tm, 1, 0, 'l');
  } else if (axis == 'y') {
    assignMatrix(element[N], element[NE], element[NW], -mp, &gFunc, tp, 0, 1, 'r');
    assignMatrix(element[S], element[SE], element[SW], mp, &gFunc, tm, 0, 1, 'l');
  }
  return row;
}

void AGM::pointStokes::makeRepresentationFormulaStokesBoundary() {
  pMatrixRow[0] = calculateRepresentationFormulaStokesBoundaryAxial('x', 0);
  pMatrixRow[1] = calculateRepresentationFormulaStokesBoundaryAxial('y', 1);
}

void AGM::pointStokes::makeRepresentationFormulaStokesInterface() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointStokes::makeRepresentationFormulaStokesInterface", "getEachMp");
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
  auto gFuncX{Greenfunction(xm, xb, xp, mpw, mpe)};
  auto gFuncY{Greenfunction(ym, yb, yp, mps, mpn)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{Greenfunction(m, b, p, d, d)};
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
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, int i, int i0, char c, char C) -> void {
    double sign = i0 ? -UNITVALUE : UNITVALUE;
    if (pt) {
      pMatrixRow[i0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      pMatrixRow[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral_tau(C, ord) : sign * func->green_integral_tau(c, ord);
      pMatrixRow[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c, ord);
      pMatrixRow[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      pMatrixRow[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
      pMatrixRow[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral_tau(C, ord), i, 1) : linearApproximation(ptr, ptl, sign * func->green_integral_tau(c, ord), i, 1);
      pMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), i, i0 + 2);
      pMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), i, i0 + 4);
    }
  };
  if (!isInterface) pMatrixRow[0][getIdx() + getNPts()] += gFuncX.green_integral_tau('c', ord);
  pMatrixRow[0][getIdx() + 2 * getNPts()] += gFuncX.green_integral_tau('c', ord);
  pMatrixRow[0][getIdx() + 4 * getNPts()] += gFuncX.green_integral_ttau('c', ord);
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

  if (!isInterface) pMatrixRow[1][getIdx() + getNPts()] += -gFuncY.green_integral_tau('c', ord);
  pMatrixRow[1][getIdx() + 3 * getNPts()] += gFuncY.green_integral_tau('c', ord);
  pMatrixRow[1][getIdx() + 5 * getNPts()] += gFuncY.green_integral_ttau('c', ord);
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFuncY, yp, 0, 1, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFuncY, ym, 0, 1, 'l', 'L');
}

auto AGM::pointStokes::calculateRepresentationFormulaStokes(const std::vector<pointStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp, int comp) -> double {
  auto assignPressureValue = [&](int i) -> double {
    double d{};
    for (const auto &item : pMatrixRow[i]) {
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
        printError("AGM::pointStokes::calculateRepresentationFormulaStokes", "item.idx (which is %d) is too large", item.idx);
      }
    }
    return d;
  };
  return assignPressureValue(comp);
}

void AGM::pointStokes::calculateDerivatives(const std::vector<pointStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp) {
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
        printError("AGM::pointStokes::calculateDerivatives", "item.idx (which is %d) is too large", item.idx);
      }
    }
    return d;
  };
  values["dx"] = assignDerivatives(0);
  values["dy"] = assignDerivatives(1);
}

void AGM::pointStokes::approximateNaNDerivatives(std::vector<pointStokes> *points) {
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
    printError("AGM::pointStokes::approximateNaNDerivatives", "findInnerPointOfBoundary");
    return nullptr;
  };
  if (std::isnan(values["dx"])) values["dx"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dx"];
  if (std::isnan(values["dy"])) values["dy"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dy"];
}

void AGM::pointStokes::approximateNaNPressure(std::vector<pointStokes> *points) {
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
    printf("values = %f\n", values["sol"]);
    printError("AGM::pointStokes::approximateNaNPressure", "findInnerPointOfBoundary");
    return nullptr;
  };
  if (std::isnan(values["sol"])) values["sol"] = points->at(findInnerPointOfBoundary()->getIdx())["sol"];
}

auto AGM::pointStokes::getPMatrixRow() const -> const std::array<matrixRow, 2> & {
  return pMatrixRow;
}

void AGM::pointStokes::setPMatrixRow(const std::array<matrixRow, 2> &matrix_row) {
  pointStokes::pMatrixRow = matrix_row;
}

auto AGM::pointStokes::getPartMatrixRow() const -> const std::array<matrixRow, 2> & {
  return partMatrixRow;
}

void AGM::pointStokes::makeRepresentationFormulaStokesFull() {
  switch (condition) {
    case 'C':
      makeRepresentationFormulaStokesCrossFull();
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      makeRepresentationFormulaStokesBoundaryFull();
      break;
    case 'I':
      makeRepresentationFormulaStokesInterfaceFull();
      break;
    default:
      printError("AGM::pointStokes::makeRepresentationFormulaStokesFull", "condition (which is %c) is wrong", condition);
  }
}

void AGM::pointStokes::makeRepresentationFormulaStokesCrossFull() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gfuncX{Greenfunction(xm, xb, xp, mp, mp)};
  auto gfuncY{Greenfunction(ym, yb, yp, mp, mp)};
  pMatrixRow[0][element[W]->getIdx()] = mp * gfuncX.green_function_ttau(xm);
  pMatrixRow[0][element[E]->getIdx()] = -mp * gfuncX.green_function_ttau(xp);

  pMatrixRow[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral_tau('l', ord);
  pMatrixRow[0][getIdx() + getNPts()] = gfuncX.green_integral_tau('c', ord);
  pMatrixRow[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral_tau('r', ord);

  prhsMatrixRow[0][element[W]->getIdx()] = gfuncX.green_integral_tau('l', ord);
  prhsMatrixRow[0][getIdx()] = gfuncX.green_integral_tau('c', ord);
  prhsMatrixRow[0][element[E]->getIdx()] = gfuncX.green_integral_tau('r', ord);

  pMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('l', ord);
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('c', ord) + UNITVALUE / mp;
  pMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('r', ord);

  pMatrixRow[1][element[S]->getIdx() + 2 * getNPts()] = mp * gfuncY.green_function_ttau(ym);
  pMatrixRow[1][element[N]->getIdx() + 2 * getNPts()] = -mp * gfuncY.green_function_ttau(yp);

  pMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = -gfuncY.green_integral_tau('l', ord);
  pMatrixRow[1][getIdx() + 3 * getNPts()] = -gfuncY.green_integral_tau('c', ord);
  pMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = -gfuncY.green_integral_tau('r', ord);

  prhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gfuncY.green_integral_tau('l', ord);
  prhsMatrixRow[1][getIdx() + getNPts()] = gfuncY.green_integral_tau('c', ord);
  prhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gfuncY.green_integral_tau('r', ord);

  pMatrixRow[1][element[S]->getIdx() + 4 * getNPts()] = gfuncY.green_integral_ttau('l', ord);
  pMatrixRow[1][getIdx() + 4 * getNPts()] = gfuncY.green_integral_ttau('c', ord) + UNITVALUE / mp;
  pMatrixRow[1][element[N]->getIdx() + 4 * getNPts()] = gfuncY.green_integral_ttau('r', ord);
}

void AGM::pointStokes::makeRepresentationFormulaStokesBoundaryFull() {
  for (const auto &item : getSolMatrixRow()[1]) {
    if (getNPts() <= item.idx && item.idx < 2 * getNPts()) {
      pMatrixRow[0][item.idx + 3 * getNPts()] = item.value;
    } else {
      printError("AGM::pointStokes::makeRepresentationFormulaStokesBoundaryFull", "phi matrix is wrong.");
    }
  }
}

void AGM::pointStokes::makeRepresentationFormulaStokesInterfaceFull() {
  auto row{std::array<matrixRow, 2>{}};
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointStokes::makeRepresentationFormulaStokesInterfaceFull", "getEachMp");
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
  auto gFuncX{Greenfunction(xm, xb, xp, mpw, mpe)};
  auto gFuncY{Greenfunction(ym, yb, yp, mps, mpn)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i, double d) -> matrixRow {
    double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
    auto func{Greenfunction(m, b, p, d, d)};
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
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, int i, int i0, char c, char C) -> void {
    double sign = i0 ? -UNITVALUE : UNITVALUE;
    if (pt) {
      row[i0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral_tau(C, ord) : sign * func->green_integral_tau(c, ord);
      row[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c, ord);
      row[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      row[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
      row[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral_tau(C, ord), i, 1) : linearApproximation(ptr, ptl, sign * func->green_integral_tau(c, ord), i, 1);
      row[i0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), i, i0 + 2);
      row[i0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), i, i0 + 4);
    }
  };
  if (!isInterface) row[0][getIdx() + getNPts()] += gFuncX.green_integral_tau('c', ord);
  row[0][getIdx() + 2 * getNPts()] += gFuncX.green_integral_tau('c', ord);
  row[0][getIdx() + 4 * getNPts()] += gFuncX.green_integral_ttau('c', ord) + UNITVALUE / mp;
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

  if (!isInterface) row[1][getIdx() + getNPts()] += -gFuncY.green_integral_tau('c', ord);
  row[1][getIdx() + 3 * getNPts()] += gFuncY.green_integral_tau('c', ord);
  row[1][getIdx() + 5 * getNPts()] += gFuncY.green_integral_ttau('c', ord) + UNITVALUE / mp;
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFuncY, yp, 0, 1, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFuncY, ym, 0, 1, 'l', 'L');

  for (int i = 0; i < 2; ++i) {
    if (!row[i].empty()) {
      for (auto &item : row[i]) {
        if (item.idx < 2 * getNPts()) {
          pMatrixRow[i][item.idx] = item.value;
        } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
          prhsMatrixRow[i][item.idx - 2 * getNPts()] = item.value;
        } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
          pMatrixRow[i][item.idx] = item.value;
        } else if (5 * getNPts() <= item.idx) {
          pMatrixRow[i][item.idx - getNPts()] += item.value;
        }
      }
    }
  }
}

void AGM::pointStokes::updateRightHandSideFull(const std::function<double(int)> &f, const std::function<double(int)> &g) {
  switch (condition) {
    case 'C':
      updateRightHandSideCrossFull(f, g);
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      break;
    case 'I':
      updateRightHandSideInterfaceFull(f, g);
      break;
    default:
      printError("AGM::pointStokes::updateRightHandSideFull", "condition (which is %d) is wrong", condition);
  }
}

void AGM::pointStokes::updateRightHandSideCrossFull(const std::function<double(int)> &f, const std::function<double(int)> &g) {
  srb[0] = srb[1] = ZEROVALUE;
  for (const auto &item : prhsMatrixRow[0]) {
    srb[0] -= item.value * f(item.idx);
  }
  for (const auto &item : prhsMatrixRow[1]) {
    srb[1] -= item.value * g(item.idx - getNPts());
  }
}

void AGM::pointStokes::updateRightHandSideInterfaceFull(const std::function<double(int)> &f, const std::function<double(int)> &g) {
  bool isinterface{solMatrixRow[1].size() == 1};
  srb[0] = srb[1] = ZEROVALUE;
  for (const auto &item : prhsMatrixRow[0]) {
    if (item.idx < getNPts()) {
      srb[0] -= item.value * f(item.idx);
    } else if (item.idx < 2 * getNPts()) {
      srb[0] -= item.value * g(item.idx - getNPts());
    } else {
      printError("AGM::pointStokes::updateRightHandSideInterfaceFull");
    }
  }
  for (const auto &item : prhsMatrixRow[1]) {
    if (item.idx < getNPts()) {
      srb[1] -= item.value * f(item.idx);
    } else if (item.idx < 2 * getNPts()) {
      srb[1] -= item.value * g(item.idx - getNPts());
    } else {
      printError("AGM::pointStokes::updateRightHandSideInterfaceFull");
    }
  }
}

auto AGM::pointStokes::getSrb() const -> const std::array<double, 2> & {
  return srb;
}

void AGM::pointStokes::setSrb(const std::array<double, 2> &array) {
  pointStokes::srb = array;
}
