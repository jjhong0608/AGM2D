//
// Created by NIMS-JUNHONG on 2022/07/20.
//

#include "pointAxisymmetric.h"

auto AGM::pointAxisymmetric::isOnAxis() const -> bool {
  return is_on_axis;
}

void AGM::pointAxisymmetric::setIsOnAxis(bool isOnAxis) {
  is_on_axis = isOnAxis;
}

void AGM::pointAxisymmetric::findStencil(const AGM::axialElement *axialElement1, std::vector<pointAxisymmetric> *vector) {
  int n{};
  for (const auto &item : *axialElement1) {
    if (item) {
      element.at(n) = &(vector->at(item->getIdx()));
    }
    ++n;
  }
}

void AGM::pointAxisymmetric::checkOnAxis() {
  if (iszero(getXy()[0])) {
    setIsOnAxis(true);
  }
}

void AGM::pointAxisymmetric::EquationOnAxis() {
  if (!isOnAxis()) return;
  std::array<matrixRow, 2> rows{};
  matrixRow row{};

  row[getIdx()] = -UNITVALUE;
  row[getElement()[E]->getIdx()] = UNITVALUE;

  setSolMatrixRow(std::array<matrixRow, 2>{row, getSolMatrixRow()[1]});
  rhsMatrixRow = std::array<matrixRow, 2>{};
  partMatrixRow = rhsMatrixRow;
  rb[0] = ZEROVALUE;
  rb[1] = ZEROVALUE;

  //    auto Error = []() -> double {
  //        printError("AGM::pointAxisymmetric::EquationOnAxis", "nullptr");
  //        return ZEROVALUE;
  //    };
  //    point *ptr{getAxialLine('x')->at(1)};
  //    double tp{ptr->getXy()[0]};
  //    auto gFunc{GreenfunctionAxisymmetric(ZEROVALUE, ZEROVALUE, tp, mp, mp)};
  //
  //    matrixRow row{};
  //    row[getIdx()] = -UNITVALUE;
  //    row[ptr->getIdx()] = UNITVALUE;
  //
  //    row[getIdx() + getNPts()] = gFunc.green_integral_ND('l') + gFunc.green_integral_ND('c');
  //    row[ptr->getIdx() + getNPts()] = gFunc.green_integral_ND('r');
  //
  //    row[getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('l') + gFunc.green_integral_ND('c');
  //    row[ptr->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('r');
  //
  //    row[getIdx() + 4 * getNPts()] =
  //            gFunc.green_integral_t_ND('l') + gFunc.green_function_ND(ZEROVALUE) + gFunc.green_integral_t_ND('c');
  //    row[ptr->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('r');
  //
  //    for (const auto &item: row) {
  //        printf("row value = %f\n", item.value);
  //    }
  //    printf("\n");
  //
  //    setSolMatrixRow(std::array<matrixRow, 2>{row, getSolMatrixRow()[1]});
  //    rhsMatrixRow = std::array<matrixRow, 2>{};
  //    partMatrixRow = rhsMatrixRow;
  //    rb[0] = ZEROVALUE;
  //    rb[1] = ZEROVALUE;
}

void AGM::pointAxisymmetric::calculateRepresentationFormulaCross() {
  std::array<matrixRow, 2> row{};
  row[0] = iszero(getElement()[W]->getXy()[0]) ? calculateRepresentationFormulaCrossSymmetricNearAxis()
                                               : calculateRepresentationFormulaCrossSymmetric();
  row[1] = calculateRepresentationFormulaCrossNonSymmetric();

  solMatrixRow[0] = row[0] + row[1];
  solMatrixRow[1] = row[0] - row[1];
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaCrossSymmetric() -> AGM::matrixRow {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mp, mp)};
  matrixRow row{};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetric::calculateRepresentationFormulaCross", "getEachMp");
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
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[W]->getIdx()] = mp * xm * gFunc.green_function_t(xm);
  row[element[E]->getIdx()] = -mp * xp * gFunc.green_function_t(xp);

  row[element[W]->getIdx() + getNPts()] = gFunc.green_integral('l', ord);
  row[getIdx() + getNPts()] = gFunc.green_integral('c', ord);
  row[element[E]->getIdx() + getNPts()] = gFunc.green_integral('r', ord);

  rhsMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral('l', ord);
  rhsMatrixRow[0][getIdx()] = gFunc.green_integral('c', ord);
  rhsMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral('r', ord);

  partMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_t('l', ord);
  partMatrixRow[0][getIdx()] = gFunc.green_integral_t('c', ord);
  partMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_t('r', ord);

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);

  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaCrossSymmetricNearAxis() -> AGM::matrixRow {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mp, mp)};
  matrixRow row{};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetric::calculateRepresentationFormulaCrossSymmetricNearAxis", "getEachMp");
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
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };
  row[getIdx()] = -UNITVALUE;
  row[element[E]->getIdx()] = UNITVALUE;

  row[element[W]->getIdx() + getNPts()] = gFunc.green_integral_ND('l', ord);
  row[getIdx() + getNPts()] = gFunc.green_integral_ND('c', ord);
  row[element[E]->getIdx() + getNPts()] = gFunc.green_integral_ND('r', ord);

  rhsMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_ND('l', ord);
  rhsMatrixRow[0][getIdx()] = gFunc.green_integral_ND('c', ord);
  rhsMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_ND('r', ord);

  partMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_t_ND('l', ord) + gFunc.green_function_ND(xm);
  partMatrixRow[0][getIdx()] = gFunc.green_integral_t_ND('c', ord);
  partMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_t_ND('r', ord);

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);
  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaCrossNonSymmetric() -> AGM::matrixRow {
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gFunc{Greenfunction(ym, yb, yp, mp * getXy()[0], mp * getXy()[0])};
  matrixRow row{};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetric::calculateRepresentationFormulaCross", "getEachMp");
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
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[S]->getIdx()] = mp * getXy()[0] * gFunc.green_function_t(ym);
  row[element[N]->getIdx()] = -mp * getXy()[0] * gFunc.green_function_t(yp);

  row[element[S]->getIdx() + getNPts()] = -gFunc.green_integral('l', ord);
  row[getIdx() + getNPts()] = -gFunc.green_integral('c', ord);
  row[element[N]->getIdx() + getNPts()] = -gFunc.green_integral('r', ord);

  rhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gFunc.green_integral('l', ord);
  rhsMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral('c', ord);
  rhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gFunc.green_integral('r', ord);

  partMatrixRow[1][element[S]->getIdx() + getNPts()] = gFunc.green_integral_t('l', ord);
  partMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral_t('c', ord);
  partMatrixRow[1][element[N]->getIdx() + getNPts()] = gFunc.green_integral_t('r', ord);

  eraseInterface(element[N], 1);
  eraseInterface(element[S], 1);
  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> AGM::matrixRow {
  if (axis == 'x') return calculateRepresentationFormulaNeumannOnAxialSymmetric();
  else if (axis == 'y')
    return calculateRepresentationFormulaNeumannOnAxialNonSymmetric();
  else
    printError("AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOnAxial", "axis (which is %c) error",
               axis);
  return {};
}

auto AGM::pointAxisymmetric::
    calculateRepresentationFormulaNeumannOnAxialSymmetric() -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOnAxialSymmetric", "nullptr");
    return ZEROVALUE;
  };
  point *ptc = getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(1) : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 1)
                                                                                                                                             : nullptr;
  point *ptl =
      getAxialLine('x')->front()->getIdx() == getIdx() ? this : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 2)
                                                                                                                : nullptr;
  point *ptr = getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(2) : getAxialLine('x')->back()->getIdx() == getIdx() ? this
                                                                                                                                             : nullptr;
  std::string string =
      getAxialLine('x')->front()->getIdx() == getIdx() ? "ND" : getAxialLine('x')->back()->getIdx() == getIdx() ? "DN"
                                                                                                                : "";
  double tm = ptl ? ptl->getXy()[0] : Error();
  double tb = ptc ? ptc->getXy()[0] : Error();
  double tp = ptr ? ptr->getXy()[0] : Error();
  auto gFunc{GreenfunctionAxisymmetric(tm, tb, tp, mp, mp)};

  matrixRow row{};
  if (string == "ND") {
    row[ptl->getIdx()] = -mp * tm * gFunc.green_function_ND(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = UNITVALUE;

    row[ptl->getIdx() + getNPts()] = gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + getNPts()] = gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + getNPts()] = gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('l', ord) + gFunc.green_function_ND(tm);
    row[ptc->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('c', ord);
    row[ptr->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('r', ord);
  } else if (string == "DN") {
    row[ptl->getIdx()] = UNITVALUE;
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = mp * tp * gFunc.green_function_DN(tp);

    row[ptl->getIdx() + getNPts()] = gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + getNPts()] = gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + getNPts()] = gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('l', ord);
    row[ptc->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('c', ord);
    row[ptr->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('r', ord) - gFunc.green_function_DN(tp);
  }
  auto c = -row[getIdx()];
  for (auto &item : row) {
    item.value /= c;
  }
  row.remove(getIdx());
  return row;
}

auto AGM::pointAxisymmetric::
    calculateRepresentationFormulaNeumannOnAxialNonSymmetric() -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOnAxialNonSymmetric", "nullptr");
    return ZEROVALUE;
  };
  point *ptc = getAxialLine('y')->front()->getIdx() == getIdx() ? getAxialLine('y')->at(1) : getAxialLine('y')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('y')->end() - 1)
                                                                                                                                             : nullptr;
  point *ptl =
      getAxialLine('y')->front()->getIdx() == getIdx() ? this : getAxialLine('y')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('y')->end() - 2)
                                                                                                                : nullptr;
  point *ptr = getAxialLine('y')->front()->getIdx() == getIdx() ? getAxialLine('y')->at(2) : getAxialLine('y')->back()->getIdx() == getIdx() ? this
                                                                                                                                             : nullptr;
  std::string string =
      getAxialLine('y')->front()->getIdx() == getIdx() ? "ND" : getAxialLine('y')->back()->getIdx() == getIdx() ? "DN"
                                                                                                                : "";
  double tm = ptl ? ptl->getXy()[1] : Error();
  double tb = ptc ? ptc->getXy()[1] : Error();
  double tp = ptr ? ptr->getXy()[1] : Error();
  auto gFunc{Greenfunction(tm, tb, tp, mp * getXy()[0], mp * getXy()[0])};

  matrixRow row{};
  if (string == "ND") {
    row[ptl->getIdx()] = -mp * tb * gFunc.green_function_ND(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = UNITVALUE;

    row[ptl->getIdx() + getNPts()] = -gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + getNPts()] = -gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + getNPts()] = -gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('l', ord);
    row[ptc->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('c', ord);
    row[ptr->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('r', ord);

    row[ptl->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('l', ord) + gFunc.green_function_ND(tm);
    row[ptc->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('c', ord);
    row[ptr->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('r', ord);
  } else if (string == "DN") {
    row[ptl->getIdx()] = UNITVALUE;
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = mp * tb * gFunc.green_function_DN(tp);

    row[ptl->getIdx() + getNPts()] = -gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + getNPts()] = -gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + getNPts()] = -gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('l', ord);
    row[ptc->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('c', ord);
    row[ptr->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('r', ord);

    row[ptl->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('l', ord);
    row[ptc->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('c', ord);
    row[ptr->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('r', ord) - gFunc.green_function_DN(tp);
  }
  auto c = -row[getIdx()];
  for (auto &item : row) {
    item.value /= c;
  }
  row.remove(getIdx());
  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> AGM::matrixRow {
  if (axis == 'x') return calculateRepresentationFormulaNeumannOffAxialSymmetric();
  else if (axis == 'y')
    return calculateRepresentationFormulaNeumannOffAxialNonSymmetric();
  else
    printError("AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOffAxial", "axis (which is %c) error",
               axis);
  return {};
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOffAxialSymmetric() -> AGM::matrixRow {
  double tm{}, tb{}, tp{};
  tm = element[W] ? element[W]->getXy()[0] : element[WN]->getXy()[0];
  tb = getXy()[0];
  tp = element[E] ? element[E]->getXy()[0] : element[EN]->getXy()[0];
  char realAxis = getAxialLine('y') ? 'y' : '\0';
  auto gFunc{GreenfunctionAxisymmetric(tm, tb, tp, mp, mp)};
  auto approximateSol = [this](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};

    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', ord);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0,
                          GreenfunctionAxisymmetric *func, double d,
                          char c) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += func->green_integral_tau(c, ord);
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, ord);
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), 4);
    }
  };
  row[getIdx() + getNPts()] += gFunc.green_integral_tau('c', ord);
  row[getIdx() + 2 * getNPts()] += gFunc.green_integral_tau('c', ord);
  row[getIdx() + 4 * getNPts()] += gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp / getXy()[0];

  assignMatrix(element[E], element[EN], element[ES], -mp * tp, &gFunc, tp, 'r');
  assignMatrix(element[W], element[WN], element[WS], mp * tm, &gFunc, tm, 'l');
  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaNeumannOffAxialNonSymmetric() -> AGM::matrixRow {
  double tm{}, tb{}, tp{};
  tm = element[S] ? element[S]->getXy()[1] : element[SE]->getXy()[1];
  tb = getXy()[1];
  tp = element[N] ? element[N]->getXy()[1] : element[NE]->getXy()[1];
  char realAxis = getAxialLine('x') ? 'x' : '\0';
  auto gFunc{Greenfunction(tm, tb, tp, mp * getXy()[0], mp * getXy()[0])};
  auto approximateSol = [this](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d * b, d * b)};
    auto mRow{matrixRow()};

    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [this, &row, &approximateSol, &linearApproximation](point *pt, point *ptr, point *ptl,
                                                                          double mp0, Greenfunction *func,
                                                                          double d,
                                                                          char c) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += -func->green_integral_tau(c, ord);
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral_tau(c, ord);
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_ttau(d), std::abs(mp0));
      row += linearApproximation(ptr, ptl, -func->green_integral_tau(c, ord), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), 5);
    }
  };
  row[getIdx() + getNPts()] += -gFunc.green_integral_tau('c', ord);
  row[getIdx() + 3 * getNPts()] += gFunc.green_integral_tau('c', ord);
  row[getIdx() + 5 * getNPts()] += gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp / getXy()[0];

  assignMatrix(element[N], element[NE], element[NW], -mp, &gFunc, tp, 'r');
  assignMatrix(element[S], element[SE], element[SW], mp, &gFunc, tm, 'l');
  return row;
}

void AGM::pointAxisymmetric::calculateRepresentationFormulaInterface() {
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    double mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    double mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    double mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    double mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  bool isInterface = checkInterface(this);

  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  std::array<matrixRow, 2> row{};
  row[0] = iszero(xm) ? calculateRepresentationFormulaInterfaceSymmetricNearAxis()
                      : calculateRepresentationFormulaInterfaceSymmetric();
  row[1] = calculateRepresentationFormulaInterfaceNonSymmetric();

  solMatrixRow[0] = row[0] + row[1];
  if (isInterface) {
    solMatrixRow[1][getIdx() + getNPts()] = UNITVALUE;
  } else {
    solMatrixRow[1] = row[0] - row[1];
  }
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaInterfaceSymmetric() -> AGM::matrixRow {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient,
                                                double d) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl,
                          double mp0,
                          GreenfunctionAxisymmetric *func,
                          double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral(C, ord)
                                                   : func->green_integral(c, ord);
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral(c, ord);
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral(C, ord), 1)
                         : linearApproximation(ptr, ptl, func->green_integral(c, ord), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, ord), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, ord), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral('c', ord);
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral('c', ord);
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t('c', ord);
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  if (!row.empty()) {
    while (row.back().idx >= 4 * getNPts()) {
      partMatrixRow[0][row.back().idx - 4 * getNPts()] = row.back().value;
      row.pop_back();
    }
    while (row.back().idx >= 2 * getNPts()) {
      rhsMatrixRow[0][row.back().idx - 2 * getNPts()] = row.back().value;
      row.pop_back();
    }
  }

  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaInterfaceSymmetricNearAxis() -> AGM::matrixRow {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient,
                                                double d) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl,
                          double mp0,
                          GreenfunctionAxisymmetric *func,
                          double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t_ND(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_ND(C, ord)
                                                   : func->green_integral_ND(c, ord);
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_ND(c, ord);
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t_ND(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t_ND(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_ND(C, ord), 1)
                         : linearApproximation(ptr, ptl, func->green_integral_ND(c, ord), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_ND(c, ord), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t_ND(c, ord), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_ND('c', ord);
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('c', ord);
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('c', ord);
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  if (!row.empty()) {
    while (row.back().idx >= 4 * getNPts()) {
      partMatrixRow[0][row.back().idx - 4 * getNPts()] = row.back().value;
      row.pop_back();
    }
    while (row.back().idx >= 2 * getNPts()) {
      rhsMatrixRow[0][row.back().idx - 2 * getNPts()] = row.back().value;
      row.pop_back();
    }
  }
  return row;
}

auto AGM::pointAxisymmetric::calculateRepresentationFormulaInterfaceNonSymmetric() -> AGM::matrixRow {
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  double mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', ord);

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', ord);
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', ord);

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', ord);
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', ord);

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', ord);
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', ord);
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  matrixRow row{};
  auto assignMatrix = [this, &row, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr,
                                                                                        point *ptl,
                                                                                        double mp0,
                                                                                        Greenfunction *func,
                                                                                        double d, char c,
                                                                                        char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, ord)
                                                   : -func->green_integral(c, ord);
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, ord);
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, ord);
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, ord), 1)
                         : linearApproximation(ptr, ptl, -func->green_integral(c, ord), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, ord), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, ord), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', ord);
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', ord);
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', ord);
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  if (!row.empty()) {
    while (row.back().idx >= 4 * getNPts()) {
      partMatrixRow[1][row.back().idx - 4 * getNPts()] = row.back().value;
      row.pop_back();
    }
    while (row.back().idx >= 2 * getNPts()) {
      rhsMatrixRow[1][row.back().idx - 2 * getNPts()] = row.back().value;
      row.pop_back();
    }
  }
  return row;
}

void AGM::pointAxisymmetric::makeDerivativesCross() {
  if (iszero(getElement()[W]->getXy()[0])) {
    makeDerivativesCrossSymmetricNearAxis();
  } else {
    makeDerivativesCrossSymmetric();
  }
  makeDerivativesCrossNonSymmetric();
}

void AGM::pointAxisymmetric::makeDerivativesCrossSymmetric() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mp, mp)};
  deriMatrixRow[0][element[W]->getIdx()] = mp * xm * gFunc.green_function_ttau(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -mp * xp * gFunc.green_function_ttau(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', ord);
  deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('l', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp / getXy()[0];
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', ord);
}

void AGM::pointAxisymmetric::makeDerivativesCrossSymmetricNearAxis() {
  double xm = element[W]->getXy()[0];
  double xb = getXy()[0];
  double xp = element[E]->getXy()[0];
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mp, mp)};
  deriMatrixRow[0][element[W]->getIdx()] = mp * xm * gFunc.green_function_ttau_ND(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -mp * xp * gFunc.green_function_ttau_ND(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('l', ord);
  deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau_ND('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('l', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('c', ord);
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('r', ord);

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('l', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('c', ord) + UNITVALUE / mp / getXy()[0];
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('r', ord);
}

void AGM::pointAxisymmetric::makeDerivativesCrossNonSymmetric() {
  double ym = element[S]->getXy()[1];
  double yb = getXy()[1];
  double yp = element[N]->getXy()[1];
  auto gFunc{Greenfunction(ym, yb, yp, mp * getXy()[0], mp * getXy()[0])};

  deriMatrixRow[1][element[S]->getIdx()] = mp * getXy()[0] * gFunc.green_function_ttau(ym);
  deriMatrixRow[1][element[N]->getIdx()] = -mp * getXy()[0] * gFunc.green_function_ttau(yp);

  deriMatrixRow[1][element[S]->getIdx() + getNPts()] = -gFunc.green_integral_tau('l', ord);
  deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', ord);
  deriMatrixRow[1][element[N]->getIdx() + getNPts()] = -gFunc.green_integral_tau('r', ord);

  deriMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('l', ord);
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('r', ord);

  deriMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('l', ord);
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp / getXy()[0];
  deriMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('r', ord);
}

void AGM::pointAxisymmetric::calculateDerivatives(const std::vector<pointAxisymmetric> *points,
                                                  const std::function<double(int)> &f,
                                                  const std::function<double(int)> &g,
                                                  const std::function<double(int)> &fp,
                                                  const std::function<double(int)> &gp) {
  auto assignDerivatives = [&](int i) -> double {
    double d{};
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
        printError("AGM::pointAxisymmetric::calculateDerivatives", "item.idx (which is %d) is too large",
                   item.idx);
      }
    }
    return d;
  };
  values["dx"] = assignDerivatives(0);
  values["dy"] = assignDerivatives(1);
}

void AGM::pointAxisymmetric::approximateNaNDerivatives(std::vector<pointAxisymmetric> *points) {
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
    printError("AGM::pointAxisymmetric::approximateNaNDerivatives", "findInnerPointOfBoundary");
    return nullptr;
  };
  if (std::isnan(values["dx"])) values["dx"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dx"];
  if (std::isnan(values["dy"])) values["dy"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dy"];
}

void AGM::pointAxisymmetric::makeDerivativesInterface() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  if (iszero(xm)) {
    makeDerivativesInterfaceSymmetricNearAxis();
  } else {
    makeDerivativesInterfaceSymmetric();
  }
  makeDerivativesInterfaceNonSymmetric();
}

void AGM::pointAxisymmetric::makeDerivativesInterfaceSymmetric() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [this, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl,
                                                                                  double mp0,
                                                                                  GreenfunctionAxisymmetric *func,
                                                                                  double d, char c,
                                                                                  char C) -> void {
    if (pt) {
      deriMatrixRow[0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      deriMatrixRow[0][pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, ord)
                                                                : func->green_integral_tau(c, ord);
      deriMatrixRow[0][pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, ord);
      deriMatrixRow[0][pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      deriMatrixRow[0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      deriMatrixRow[0] += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, ord), 1)
                                      : linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 1);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 2);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp / getXy()[0];
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');
}

void AGM::pointAxisymmetric::makeDerivativesInterfaceSymmetricNearAxis() {
  double xm = getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0];
  double xb = getXy()[0];
  double xp = getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  double mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe)};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', ord);

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [this, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl,
                                                                                  double mp0,
                                                                                  GreenfunctionAxisymmetric *func,
                                                                                  double d, char c,
                                                                                  char C) -> void {
    if (pt) {
      deriMatrixRow[0][pt->getIdx()] += mp0 * func->green_function_ttau_ND(d);
      deriMatrixRow[0][pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau_ND(C, ord)
                                                                : func->green_integral_tau_ND(c, ord);
      deriMatrixRow[0][pt->getIdx() + 2 * getNPts()] += func->green_integral_tau_ND(c, ord);
      deriMatrixRow[0][pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau_ND(c, ord);
    } else {
      deriMatrixRow[0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau_ND(d), std::abs(mp0));
      deriMatrixRow[0] += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau_ND(C, ord), 1)
                                      : linearApproximation(ptr, ptl, func->green_integral_tau_ND(c, ord), 1);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_tau_ND(c, ord), 2);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_ttau_ND(c, ord), 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau_ND('c', ord);
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('c', ord);
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('c', ord) + UNITVALUE / mp / getXy()[0];
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');
}

void AGM::pointAxisymmetric::makeDerivativesInterfaceNonSymmetric() {
  double ym = getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1];
  double yb = getXy()[1];
  double yp = getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1];
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetric::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
        : ptl->getCondition() == 'C'                           ? ptl->getMp()
                                                               : Error();
    return rtv;
  };
  double mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  double mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
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
  auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', ord);

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', ord);
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', ord);

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', ord);
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', ord);

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', ord);
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', ord);

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', ord);
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', ord);

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', ord);
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', ord);
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    double m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [this, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl,
                                                                                  double mp0, Greenfunction *func,
                                                                                  double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[1][pt->getIdx()] += mp0 * getXy()[0] * func->green_function_ttau(d);
      deriMatrixRow[1][pt->getIdx() + getNPts()] += isInterface ? -func->green_integral_tau(C, ord)
                                                                : -func->green_integral_tau(c, ord);
      deriMatrixRow[1][pt->getIdx() + 3 * getNPts()] += func->green_integral_tau(c, ord);
      deriMatrixRow[1][pt->getIdx() + 5 * getNPts()] += func->green_integral_ttau(c, ord);
    } else {
      deriMatrixRow[1] += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_ttau(d),
                                         std::abs(mp0));
      deriMatrixRow[1] += isInterface ? linearApproximation(ptr, ptl, -func->green_integral_tau(C, ord), 1)
                                      : linearApproximation(ptr, ptl, -func->green_integral_tau(c, ord), 1);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_tau(c, ord), 3);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, ord), 5);
    }
  };
  if (!isInterface) deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', ord);
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', ord);
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('c', ord) + UNITVALUE / mp;
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');
}
