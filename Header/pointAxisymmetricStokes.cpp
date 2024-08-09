//
// Created by 조준홍 on 2/19/24.
//

#include "pointAxisymmetricStokes.h"

auto AGM::pointAxisymmetricStokes::isOnAxis() const -> bool {
  return is_on_axis;
}

void AGM::pointAxisymmetricStokes::setIsOnAxis(bool isOnAxis) {
  is_on_axis = isOnAxis;
}

void AGM::pointAxisymmetricStokes::findStencil(const AGM::axialElement *axialElement1, std::vector<pointAxisymmetricStokes> *vector) {
  auto n{0};
  for (const auto &item : *axialElement1) {
    if (item) {
      element.at(n) = &(vector->at(item->getIdx()));
    }
    ++n;
  }
}

void AGM::pointAxisymmetricStokes::checkOnAxis() {
  if (iszero(getXy()[0]) && getCondition() == 'N') {
    setIsOnAxis(true);
  }
}

void AGM::pointAxisymmetricStokes::EquationOnAxis() {
  if (!isOnAxis()) return;
  matrixRow row{};

  row[getIdx()] = -UNITVALUE;
  row[getElement()[E]->getIdx()] = UNITVALUE;

  setSolMatrixRow(std::array<matrixRow, 2>{row, getSolMatrixRow()[1]});
  rhsMatrixRow = std::array<matrixRow, 2>{};
  partMatrixRow = rhsMatrixRow;
  rb[0] = ZEROVALUE;
  rb[1] = ZEROVALUE;
}

void AGM::pointAxisymmetricStokes::calculateRepresentationFormula(char rz, int order) {
  switch (condition) {
    case 'C':
      calculateRepresentationFormulaCross(rz);
      break;
    case 'D':
    case 'd':
      calculateRepresentationFormulaDirichlet(order);
      break;
    case 'N':
      calculateRepresentationFormulaNeumann(order, rz);
      break;
    case 'n':
      calculateRepresentationFormulaInterfaceNeumann();
      break;
    case 'I':
      calculateRepresentationFormulaInterface(rz);
      break;
    default:
      printError("AGM::pointAxisymmetricStokes::calcRepresentationFormula", "boundary condition (which is %c) is wrong", condition);
  }
}

void AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCross(char rz) {
  auto row{std::array<matrixRow, 2>()};
  if (rz == 'r') {
    row[0] = calculateRepresentationFormulaCrossSymmetricReactionDiffusion();
    row[1] = calculateRepresentationFormulaCrossNonSymmetricR();
  } else if (rz == 'z') {
    row[0] = iszero(getElement()[W]->getXy()[0]) ? calculateRepresentationFormulaCrossSymmetricNearAxis() : calculateRepresentationFormulaCrossSymmetricElliptic();
    row[1] = calculateRepresentationFormulaCrossNonSymmetricZ();
  } else {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCross", "Argument `rz` must be `r` or `z`, got %c.", rz);
  }
  solMatrixRow[0] = row[0] + row[1];
  solMatrixRow[1] = row[0] - row[1];
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricElliptic() -> AGM::matrixRow {
  auto xm{element[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{element[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  auto row{matrixRow()};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricElliptic", "getEachMp");
          return ZEROVALUE;
        };
        auto rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                         ? ptl->getMp()
                                                                 : Error();
        return rtv;
      };
      auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[W]->getIdx()] = getMp() * xm * gFunc.green_function_t(xm);
  row[element[E]->getIdx()] = -getMp() * xp * gFunc.green_function_t(xp);

  row[element[W]->getIdx() + getNPts()] = gFunc.green_integral('l', getOrd());
  row[getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  row[element[E]->getIdx() + getNPts()] = gFunc.green_integral('r', getOrd());

  rhsMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral('l', getOrd());
  rhsMatrixRow[0][getIdx()] = gFunc.green_integral('c', getOrd());
  rhsMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral('r', getOrd());

  partMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_t('l', getOrd());
  partMatrixRow[0][getIdx()] = gFunc.green_integral_t('c', getOrd());
  partMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_t('r', getOrd());

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);

  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricReactionDiffusion() -> AGM::matrixRow {
  auto xm{element[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{element[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetricStokes(xm, xb, xp, getMp(), getMp())};
  auto row{matrixRow()};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricReactionDiffusion", "getEachMp");
          return ZEROVALUE;
        };
        auto rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                         ? ptl->getMp()
                                                                 : Error();
        return rtv;
      };
      auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[W]->getIdx()] = getMp() * xm * gFunc.green_function_t(xm);
  row[element[E]->getIdx()] = -getMp() * xp * gFunc.green_function_t(xp);

  row[element[W]->getIdx() + getNPts()] = gFunc.green_integral('l', getOrd());
  row[getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  row[element[E]->getIdx() + getNPts()] = gFunc.green_integral('r', getOrd());

  rhsMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral('l', getOrd());
  rhsMatrixRow[0][getIdx()] = gFunc.green_integral('c', getOrd());
  rhsMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral('r', getOrd());

  partMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_t('l', getOrd()) + gFunc.green_integral('l', getOrd());
  partMatrixRow[0][getIdx()] = gFunc.green_integral_t('c', getOrd()) + gFunc.green_integral('c', getOrd());
  partMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_t('r', getOrd()) + gFunc.green_integral('r', getOrd());

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);

  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricNearAxis() -> AGM::matrixRow {
  auto xm{element[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{element[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  auto row{matrixRow()};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossSymmetricNearAxis", "getEachMp");
          return ZEROVALUE;
        };
        auto rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                         ? ptl->getMp()
                                                                 : Error();
        return rtv;
      };
      auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };
  row[getIdx()] = -UNITVALUE;
  row[element[E]->getIdx()] = UNITVALUE;

  row[element[W]->getIdx() + getNPts()] = gFunc.green_integral_ND('l', getOrd());
  row[getIdx() + getNPts()] = gFunc.green_integral_ND('c', getOrd());
  row[element[E]->getIdx() + getNPts()] = gFunc.green_integral_ND('r', getOrd());

  rhsMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_ND('l', getOrd());
  rhsMatrixRow[0][getIdx()] = gFunc.green_integral_ND('c', getOrd());
  rhsMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_ND('r', getOrd());

  partMatrixRow[0][element[W]->getIdx()] = gFunc.green_integral_t_ND('l', getOrd()) + gFunc.green_function_ND(xm);
  partMatrixRow[0][getIdx()] = gFunc.green_integral_t_ND('c', getOrd());
  partMatrixRow[0][element[E]->getIdx()] = gFunc.green_integral_t_ND('r', getOrd());

  eraseInterface(element[E], 0);
  eraseInterface(element[W], 0);
  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossNonSymmetricR() -> AGM::matrixRow {
  auto ym{element[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{element[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  auto row{matrixRow()};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossNonSymmetricR", "getEachMp");
          return ZEROVALUE;
        };
        auto rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                         ? ptl->getMp()
                                                                 : Error();
        return rtv;
      };
      auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[S]->getIdx()] = getMp() * getXy()[0] * gFunc.green_function_t(ym);
  row[element[N]->getIdx()] = -getMp() * getXy()[0] * gFunc.green_function_t(yp);

  row[element[S]->getIdx() + getNPts()] = -gFunc.green_integral('l', getOrd());
  row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[element[N]->getIdx() + getNPts()] = -gFunc.green_integral('r', getOrd());

  rhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gFunc.green_integral('l', getOrd());
  rhsMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  rhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gFunc.green_integral('r', getOrd());

  partMatrixRow[1][element[S]->getIdx() + getNPts()] = gFunc.green_integral_t('l', getOrd());
  partMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral_t('c', getOrd());
  partMatrixRow[1][element[N]->getIdx() + getNPts()] = gFunc.green_integral_t('r', getOrd());

  eraseInterface(element[N], 1);
  eraseInterface(element[S], 1);

  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossNonSymmetricZ() -> AGM::matrixRow {
  auto ym{element[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{element[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  auto row{matrixRow()};
  auto eraseInterface = [this, &row](point *pt, int i) -> void {
    auto checkInterface = [](point *pt) -> bool {
      auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
        auto Error = []() -> double {
          printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaCrossNonSymmetricZ", "getEachMp");
          return ZEROVALUE;
        };
        auto rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
            : ptl->getCondition() == 'C'                         ? ptl->getMp()
                                                                 : Error();
        return rtv;
      };
      auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
      auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
      auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
      auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
      return pt->getCondition() == 'I' && !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
    };
    if (checkInterface(pt)) {
      row[getIdx() + getNPts()] += row[pt->getIdx() + getNPts()];
      row.remove(pt->getIdx() + getNPts());
    }
  };

  row[getIdx()] = -UNITVALUE;
  row[element[S]->getIdx()] = getMp() * getXy()[0] * gFunc.green_function_t(ym);
  row[element[N]->getIdx()] = -getMp() * getXy()[0] * gFunc.green_function_t(yp);

  row[element[S]->getIdx() + getNPts()] = -gFunc.green_integral('l', getOrd());
  row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[element[N]->getIdx() + getNPts()] = -gFunc.green_integral('r', getOrd());

  rhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gFunc.green_integral('l', getOrd());
  rhsMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  rhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gFunc.green_integral('r', getOrd());

  partMatrixRow[1][element[S]->getIdx() + getNPts()] = getXy()[0] * gFunc.green_integral_t('l', getOrd());
  partMatrixRow[1][getIdx() + getNPts()] = getXy()[0] * gFunc.green_integral_t('c', getOrd());
  partMatrixRow[1][element[N]->getIdx() + getNPts()] = getXy()[0] * gFunc.green_integral_t('r', getOrd());

  eraseInterface(element[N], 1);
  eraseInterface(element[S], 1);

  return row;
}

void AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumann(int order, char rz) {
  auto row{std::array<matrixRow, 2>()};
  row[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial(rz, 'x', 0) : calculateRepresentationFormulaNeumannOffAxial(rz, 'x', 0);
  row[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial(rz, 'y', 1) : calculateRepresentationFormulaNeumannOffAxial(rz, 'y', 1);
  for (int i = 0; i < 2; ++i) {
    if (!row[i].empty() && !iszero(normal[i])) {
      while (row[i].back().idx >= 4 * getNPts()) {
        partMatrixRow[i][row[i].back().idx - 4 * getNPts()] = row[i].back().value * normal[i];
        row[i].pop_back();
      }
      if (i == 1) {
        for (auto &item : partMatrixRow[1]) {
          item.value *= getXy()[0];
        }
      }
      while (row[i].back().idx >= 2 * getNPts()) {
        rhsMatrixRow[i][row[i].back().idx - 2 * getNPts()] = row[i].back().value * normal[i];
        row[i].pop_back();
      }
      solMatrixRow[0] += row[i] * normal[i];
    }
  }
  auto findInnerPointOfBoundary = [this](int i) -> point * {
    for (const auto &item : {'x', 'y'}) {
      if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
        return getAxialLine(item)->at(i);
      }
      if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
        return *std::prev(getAxialLine(item)->end() - i);
      }
    }
    printInformation();
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumann", "findInnerPointOfBoundary");
    return {};
  };
  auto pt1{*findInnerPointOfBoundary(1)};
  auto pt2{order > 0 ? *findInnerPointOfBoundary(2) : point()};
  if (order == 2) {
    if (min(*this - pt1, pt1 - pt2) / max(*this - pt1, pt1 - pt2) < 0.1) {
      order = 0;
    } else if (min(*this - pt1, pt1 - pt2) / max(*this - pt1, pt1 - pt2) < 0.2) {
      order = 1;
    }
  }
  approximatePhiAtBoundary1(order);
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxial(char rz, char axis, int axisInt) -> AGM::matrixRow {
  if (rz == 'r') {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxial", "Neumann condition is not implemented");
  } else if (rz == 'z') {
    if (axis == 'x') {
      return calculateRepresentationFormulaNeumannOnAxialSymmetricElliptic();
    } else if (axis == 'y') {
      return calculateRepresentationFormulaNeumannOnAxialNonSymmetricElliptic();
    } else {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxial", "axis (which is %c) error", axis);
    }
  } else {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxial", "input argument rz (which is %c) must be 'r' or 'z'.", rz);
  }
  return {};
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxialSymmetricElliptic() -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxialSymmetric", "nullptr");
    return ZEROVALUE;
  };
  point *ptc{getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(1) : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 1)
                                                                                                                                           : nullptr};
  point *ptl{getAxialLine('x')->front()->getIdx() == getIdx() ? this : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 2)
                                                                                                                       : nullptr};
  point *ptr{getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(2) : getAxialLine('x')->back()->getIdx() == getIdx() ? this
                                                                                                                                           : nullptr};
  std::string string{getAxialLine('x')->front()->getIdx() == getIdx() ? "ND" : getAxialLine('y')->back()->getIdx() == getIdx() ? "DN"
                                                                                                                               : ""};
  auto tm{ptl ? ptl->getXy()[0] : Error()};
  auto tb{ptc ? ptc->getXy()[0] : Error()};
  auto tp{ptr ? ptr->getXy()[0] : Error()};
  auto gFunc{GreenfunctionAxisymmetric(tm, tb, tp, getMp(), getMp())};
  auto row{matrixRow()};

  if (string == "ND") {
    row[ptl->getIdx()] = -getMp() * tm * gFunc.green_function_ND(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = UNITVALUE;

    row[ptl->getIdx() + getNPts()] = gFunc.green_integral_ND('l', getOrd());
    row[ptc->getIdx() + getNPts()] = gFunc.green_integral_ND('c', getOrd());
    row[ptr->getIdx() + getNPts()] = gFunc.green_integral_ND('r', getOrd());

    row[ptl->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('l', getOrd());
    row[ptc->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('c', getOrd());
    row[ptr->getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('r', getOrd());

    row[ptl->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('l', getOrd()) + gFunc.green_function_ND(tm);
    row[ptc->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('c', getOrd());
    row[ptr->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('r', getOrd());
  } else if (string == "DN") {
    row[ptl->getIdx()] = UNITVALUE;
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = getMp() * tp * gFunc.green_function_DN(tp);

    row[ptl->getIdx() + getNPts()] = gFunc.green_integral_DN('l', getOrd());
    row[ptc->getIdx() + getNPts()] = gFunc.green_integral_DN('c', getOrd());
    row[ptr->getIdx() + getNPts()] = gFunc.green_integral_DN('r', getOrd());

    row[ptl->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('l', getOrd());
    row[ptc->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('c', getOrd());
    row[ptr->getIdx() + 2 * getNPts()] = gFunc.green_integral_DN('r', getOrd());

    row[ptl->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('l', getOrd());
    row[ptc->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('c', getOrd());
    row[ptr->getIdx() + 4 * getNPts()] = gFunc.green_integral_t_DN('r', getOrd()) - gFunc.green_function_DN(tp);
  }
  auto c = -row[getIdx()];
  for (auto &item : row) {
    item.value /= c;
  }
  row.remove(getIdx());
  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxialNonSymmetricElliptic() -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOnAxialSymmetric", "nullptr");
    return ZEROVALUE;
  };
  point *ptc{getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(1) : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 1)
                                                                                                                                           : nullptr};
  point *ptl{getAxialLine('x')->front()->getIdx() == getIdx() ? this : getAxialLine('x')->back()->getIdx() == getIdx() ? *std::prev(getAxialLine('x')->end() - 2)
                                                                                                                       : nullptr};
  point *ptr{getAxialLine('x')->front()->getIdx() == getIdx() ? getAxialLine('x')->at(2) : getAxialLine('x')->back()->getIdx() == getIdx() ? this
                                                                                                                                           : nullptr};
  std::string string{getAxialLine('x')->front()->getIdx() == getIdx() ? "ND" : getAxialLine('y')->back()->getIdx() == getIdx() ? "DN"
                                                                                                                               : ""};
  auto tm{ptl ? ptl->getXy()[1] : Error()};
  auto tb{ptc ? ptc->getXy()[1] : Error()};
  auto tp{ptr ? ptr->getXy()[1] : Error()};
  auto gFunc{Greenfunction(tm, tb, tp, getMp() * getXy()[0], getMp() * getXy()[0])};
  auto row{matrixRow()};

  if (string == "ND") {
    row[ptl->getIdx()] = -getMp() * tb * gFunc.green_function_ND(tm);
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = UNITVALUE;

    row[ptl->getIdx() + getNPts()] = -gFunc.green_integral_ND('l', getOrd());
    row[ptc->getIdx() + getNPts()] = -gFunc.green_integral_ND('c', getOrd());
    row[ptr->getIdx() + getNPts()] = -gFunc.green_integral_ND('r', getOrd());

    row[ptl->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('l', getOrd());
    row[ptc->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('c', getOrd());
    row[ptr->getIdx() + 3 * getNPts()] = gFunc.green_integral_ND('r', getOrd());

    row[ptl->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('l', getOrd()) + gFunc.green_function_ND(tm);
    row[ptc->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('c', getOrd());
    row[ptr->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_ND('r', getOrd());
  } else if (string == "DN") {
    row[ptl->getIdx()] = UNITVALUE;
    row[ptc->getIdx()] = -UNITVALUE;
    row[ptr->getIdx()] = getMp() * tb * gFunc.green_function_DN(tp);

    row[ptl->getIdx() + getNPts()] = -gFunc.green_integral_DN('l', getOrd());
    row[ptc->getIdx() + getNPts()] = -gFunc.green_integral_DN('c', getOrd());
    row[ptr->getIdx() + getNPts()] = -gFunc.green_integral_DN('r', getOrd());

    row[ptl->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('l', getOrd());
    row[ptc->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('c', getOrd());
    row[ptr->getIdx() + 3 * getNPts()] = gFunc.green_integral_DN('r', getOrd());

    row[ptl->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('l', getOrd());
    row[ptc->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('c', getOrd());
    row[ptr->getIdx() + 5 * getNPts()] = gFunc.green_integral_t_DN('r', getOrd()) - gFunc.green_function_DN(tp);
  }
  auto c = -row[getIdx()];
  for (auto &item : row) {
    item.value /= c;
  }
  row.remove(getIdx());
  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxial(char rz,
                                                                                 char axis,
                                                                                 int axisInt) -> AGM::matrixRow {
  if (rz == 'r') {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxial", "Neumann condition is not constructed");
  } else if (rz == 'z') {
    if (axis == 'x') {
      return calculateRepresentationFormulaNeumannOffAxialSymmetricElliptic();
    } else if (axis == 'y') {
      return calculateRepresentationFormulaNeumannOffAxialNonSymmetricElliptic();
    } else {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxial", "axis (which is %c) error", axis);
    }
  } else {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxial", "input argument rz (which is %c) must be 'r' or 'z'.", rz);
  }
  return {};
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxialSymmetricElliptic() -> AGM::matrixRow {
  auto Error = []() -> double {
    printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxialSymmetricElliptic", "nullptr");
    return ZEROVALUE;
  };
  auto tm{element[W] ? element[W]->getXy()[0] : element[WN]->getXy()[0]};
  auto tb{getXy()[0]};
  auto tp{element[E] ? element[E]->getXy()[0] : element[EN]->getXy()[0]};
  auto realAxis{getAxialLine('y') ? 'y' : '\0'};
  auto gFunc{GreenfunctionAxisymmetric(tm, tb, tp, getMp(), getMp())};

  auto approximateSol = [this](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};

    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 4);
    }
  };
  row[getIdx() + getNPts()] += gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 2 * getNPts()] += gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 4 * getNPts()] += gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];

  assignMatrix(element[E], element[EN], element[ES], -getMp() * tp, &gFunc, tp, 'r');
  assignMatrix(element[W], element[WN], element[WS], getMp() * tm, &gFunc, tm, 'l');
  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaNeumannOffAxialNonSymmetricElliptic() -> AGM::matrixRow {
  auto tm{element[S] ? element[S]->getXy()[1] : element[SE]->getXy()[1]};
  auto tb{getXy()[1]};
  auto tp{element[N] ? element[N]->getXy()[1] : element[NE]->getXy()[1]};
  auto realAxis{getAxialLine('x') ? 'x' : '\0'};
  auto gFunc{Greenfunction(tm, tb, tp, getMp() * getXy()[0], getMp() * getXy()[0])};

  auto approximateSol = [this](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d * b, d * b)};
    auto mRow{matrixRow()};

    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += -func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_ttau(d), std::abs(mp0));
      row += linearApproximation(ptr, ptl, -func->green_integral_tau(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 5);
    }
  };
  row[getIdx() + getNPts()] += -gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 3 * getNPts()] += gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 5 * getNPts()] += gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];

  assignMatrix(element[N], element[NE], element[NW], -getMp(), &gFunc, tp, 'r');
  assignMatrix(element[S], element[SE], element[SW], getMp(), &gFunc, tm, 'l');
  return row;
}

void AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterface(char rz) {
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterface", "getEachMp");
      return ZEROVALUE;
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto checkInterface = [&](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto row{std::array<matrixRow, 2>()};
  row[0] = rz == 'r' ? calculateRepresentationFormulaInterfaceSymmetricReactionDiffusion() : iszero(xm) ? calculateRepresentationFormulaInterfaceSymmetricNearAxis()
                                                                                                        : calculateRepresentationFormulaInterfaceSymmetricElliptic();
  row[1] = rz == 'r' ? calculateRepresentationFormulaInterfaceNonSymmetricR() : calculateRepresentationFormulaInterfaceNonSymmetricZ();
  solMatrixRow[0] = row[0] + row[1];
  if (isInterface) {
    solMatrixRow[1][getIdx() + getNPts()] = UNITVALUE;
  } else {
    solMatrixRow[1] = row[0] - row[1];
  }
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricElliptic() -> AGM::matrixRow {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricElliptic", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral(C, getOrd()) : func->green_integral(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t('c', getOrd());
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

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricReactionDiffusion() -> AGM::matrixRow {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricReactionDiffusion", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetricStokes(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetricStokes *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral(C, getOrd()) : func->green_integral(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t(c, getOrd()) + func->green_integral(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()) + func->green_integral(c, getOrd()), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t('c', getOrd()) + gFunc.green_integral('c', getOrd());
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

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricNearAxis() -> AGM::matrixRow {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceSymmetricNearAxis", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t_ND(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_ND(C, getOrd()) : func->green_integral_ND(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_ND(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t_ND(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t_ND(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_ND(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_ND(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_ND(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t_ND(c, getOrd()), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_ND('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('c', getOrd());
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

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceNonSymmetricR() -> AGM::matrixRow {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceNonSymmetricR", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetricStokes(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    //        if (iszero(m)) {
    //            auto mRow0{matrixRow()};
    //            mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
    //            mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);
    //
    //            mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());
    //
    //            checkMatrixRow(&mRow0, ptr, ptl);
    //
    //            mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());
    //
    //            mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
    //            return mRow0 * coefficient;
    //        }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, getOrd()) : -func->green_integral(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, getOrd()), 1) : linearApproximation(ptr, ptl, -func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', getOrd());
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

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceNonSymmetricZ() -> AGM::matrixRow {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaInterfaceNonSymmetricZ", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, getOrd()) : -func->green_integral(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, getOrd()), 1) : linearApproximation(ptr, ptl, -func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', getOrd());
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  if (!row.empty()) {
    while (row.back().idx >= 4 * getNPts()) {
      partMatrixRow[1][row.back().idx - 4 * getNPts()] = getXy()[0] * row.back().value;
      row.pop_back();
    }
    while (row.back().idx >= 2 * getNPts()) {
      rhsMatrixRow[1][row.back().idx - 2 * getNPts()] = row.back().value;
      row.pop_back();
    }
  }
  return row;
}

auto AGM::pointAxisymmetricStokes::calculateRepresentationFormulaStokes(const std::vector<pointAxisymmetricStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp, int comp) -> double {
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
        printError("AGM::pointAxisymmetricStokes::calculateRepresentationFormulaStokes", "item.idx (which is %d) is too large", item.idx);
      }
    }
    return d;
  };
  return assignPressureValue(comp);
}

void AGM::pointAxisymmetricStokes::approximateNaNPressure(std::vector<pointAxisymmetricStokes> *points) {
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
    printf("values = %f\n", getValue()["sol"]);
    printError("AGM::pointAxisymmetricStokes::approximateNaNPressure", "findInnerPointOfBoundary");
    return {};
  };
  if (std::isnan(getValue()["sol"])) values["sol"] = points->at(findInnerPointOfBoundary()->getIdx())["sol"];
}

void AGM::pointAxisymmetricStokes::constantExtensionPressure(std::vector<pointAxisymmetricStokes> *points) {
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
    printf("value = %f\n", getValue()["so"]);
    printError("AGM::pointAxisymmetricStokes::constantExtensionPressure", "findInnerPointOfBoundary");
    return {};
  };
  if (getCondition() == 'D' || getCondition() == 'N') values["sol"] = points->at(findInnerPointOfBoundary()->getIdx())["sol"];
}

void AGM::pointAxisymmetricStokes::makeDerivatives(char rz) {
  switch (condition) {
    case 'C':
      makeDerivativesCross(rz);
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      makeDerivativesBoundary(rz);
      break;
    case 'I':
      makeDerivativesInterface(rz);
      break;
    default: printError("AGM::pointAxisymmetricStokes::makeDerivatives", "condition (which is %c) is wrong", condition);
  }
}

void AGM::pointAxisymmetricStokes::makeDerivativesCross(char rz) {
  if (rz == 'r') {
    makeDerivativesCrossSymmetricReactionDiffusion();
    makeDerivativesCrossNonSymmetricR();
  } else if (rz == 'z') {
    if (iszero(getElement()[W]->getXy()[0])) {
      makeDerivativesCrossSymmetricNearAxis();
    } else {
      makeDerivativesCrossSymmetricElliptic();
    }
    makeDerivativesCrossNonSymmetricZ();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeDerivativesCross", "input argument rz (which is %c) must be 'r' or 'z'.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeDerivativesCrossSymmetricElliptic() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  deriMatrixRow[0][element[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeDerivativesCrossSymmetricReactionDiffusion() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetricStokes(xm, xb, xp, getMp(), getMp())};
  deriMatrixRow[0][element[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd()) + gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd()) + gFunc.green_integral_tau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeDerivativesCrossSymmetricNearAxis() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  deriMatrixRow[0][element[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau_ND(xm);
  deriMatrixRow[0][element[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau_ND(xp);

  deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('l', getOrd());
  deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau_ND('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('l', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('c', getOrd());
  deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('r', getOrd());

  deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('l', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeDerivativesCrossNonSymmetricR() {
  auto ym{getElement()[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  deriMatrixRow[1][element[S]->getIdx()] = getMp() * getXy()[0] * gFunc.green_function_ttau(ym);
  deriMatrixRow[1][element[N]->getIdx()] = -getMp() * getXy()[0] * gFunc.green_function_ttau(yp);

  deriMatrixRow[1][element[S]->getIdx() + getNPts()] = -gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][element[N]->getIdx() + getNPts()] = -gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('l', getOrd());
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  deriMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeDerivativesCrossNonSymmetricZ() {
  auto ym{getElement()[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  deriMatrixRow[1][element[S]->getIdx()] = getMp() * getXy()[0] * gFunc.green_function_ttau(ym);
  deriMatrixRow[1][element[N]->getIdx()] = -getMp() * getXy()[0] * gFunc.green_function_ttau(yp);

  deriMatrixRow[1][element[S]->getIdx() + getNPts()] = -gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][element[N]->getIdx() + getNPts()] = -gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  deriMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('l', getOrd());
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp();
  deriMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeDerivativesBoundary(char rz) {
  deriMatrixRow[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial(rz, 'x', 0)
                                       : calculateRepresentationFormulaNeumannOffAxial(rz, 'x', 0);
  deriMatrixRow[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial(rz, 'y', 1)
                                       : calculateRepresentationFormulaNeumannOffAxial(rz, 'y', 1);
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterface(char rz) {
  if (rz == 'r') {
    makeDerivativesInterfaceSymmetricReactionDiffusion();
    makeDerivativesInterfaceNonSymmetricR();
  } else if (rz == 'z') {
    auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
    if (iszero(xm)) {
      makeDerivativesInterfaceSymmetricNearAxis();
    } else {
      makeDerivativesInterfaceSymmetricElliptic();
    }
    makeDerivativesInterfaceNonSymmetricZ();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeDerivativesInterface", "input argument rz (which is %c) must be `r` or `z`.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricElliptic() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricElliptic", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [this, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      deriMatrixRow[0][pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, getOrd())
                                                                : func->green_integral_tau(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      deriMatrixRow[0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      deriMatrixRow[0] += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, getOrd()), 1)
                                      : linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricReactionDiffusion() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricReactionDiffusion", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetricStokes(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetricStokes *func, double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
      deriMatrixRow[0][pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, getOrd()) : func->green_integral_tau(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 4 * getNPts()] +=
          func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd());
    } else {
      deriMatrixRow[0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      deriMatrixRow[0] += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd()), 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricNearAxis() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeDerivativesInterfaceSymmetricNearAxis", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[0][pt->getIdx()] += mp0 * func->green_function_ttau_ND(d);
      deriMatrixRow[0][pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau_ND(C, getOrd()) : func->green_integral_tau_ND(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 2 * getNPts()] += func->green_integral_tau_ND(c, getOrd());
      deriMatrixRow[0][pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau_ND(c, getOrd());
    } else {
      deriMatrixRow[0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau_ND(d), std::abs(mp0));
      deriMatrixRow[0] += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau_ND(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_tau_ND(c, getOrd()), 1);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_tau_ND(c, getOrd()), 2);
      deriMatrixRow[0] += linearApproximation(ptr, ptl, func->green_integral_ttau_ND(c, getOrd()), 4);
    }
  };
  if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau_ND('c', getOrd());
  deriMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau_ND('c', getOrd());
  deriMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterfaceNonSymmetricR() {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeDerivativesInterfaceNonSymmetricR", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetricStokes(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    //        if (iszero(m)) {
    //            auto mRow0{matrixRow()};
    //            mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
    //            mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);
    //
    //            mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());
    //
    //            checkMatrixRow(&mRow0, ptr, ptl);
    //
    //            mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());
    //
    //            mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
    //            mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
    //            return mRow0 * coefficient;
    //        }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[1][pt->getIdx()] += mp0 * getXy()[0] * func->green_function_ttau(d);
      deriMatrixRow[1][pt->getIdx() + getNPts()] += isInterface ? -func->green_integral_tau(C, getOrd()) : -func->green_integral_tau(c, getOrd());
      deriMatrixRow[1][pt->getIdx() + 3 * getNPts()] += func->green_integral_tau(c, getOrd());
      deriMatrixRow[1][pt->getIdx() + 5 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      deriMatrixRow[1] += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_ttau(d), std::abs(mp0));
      deriMatrixRow[1] += isInterface ? linearApproximation(ptr, ptl, -func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, -func->green_integral_tau(c, getOrd()), 1);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 3);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 5);
    }
  };
  if (!isInterface) deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp();
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');
}

void AGM::pointAxisymmetricStokes::makeDerivativesInterfaceNonSymmetricZ() {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeDerivativesInterfaceNonSymmetricZ", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      deriMatrixRow[1][pt->getIdx()] += mp0 * getXy()[0] * func->green_function_ttau(d);
      deriMatrixRow[1][pt->getIdx() + getNPts()] += isInterface ? -func->green_integral_tau(C, getOrd()) : -func->green_integral_tau(c, getOrd());
      deriMatrixRow[1][pt->getIdx() + 3 * getNPts()] += func->green_integral_tau(c, getOrd());
      deriMatrixRow[1][pt->getIdx() + 5 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      deriMatrixRow[1] += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_ttau(d), std::abs(mp0));
      deriMatrixRow[1] += isInterface ? linearApproximation(ptr, ptl, -func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, -func->green_integral_tau(c, getOrd()), 1);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 3);
      deriMatrixRow[1] += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 5);
    }
  };
  if (!isInterface) deriMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  deriMatrixRow[1][getIdx() + 5 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  for (auto &item : deriMatrixRow[1]) {
    if (item.idx > 4 * getNPts()) {
      item.value *= getXy()[0];
    }
  }
}

void AGM::pointAxisymmetricStokes::calculateDerivatives(const std::vector<pointAxisymmetricStokes> *points, const std::function<double(int)> &f, const std::function<double(int)> &g, const std::function<double(int)> &fp, const std::function<double(int)> &gp) {
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
        printError("AGM::pointAxisymmetricStokes::calculateDerivatives", "item.idx (which is %d) is too large", item.idx);
      }
    }
    return d;
  };
  values["dx"] = assignDerivatives(0);
  values["dy"] = assignDerivatives(1);
}

void AGM::pointAxisymmetricStokes::approximateNaNDerivatives(std::vector<pointAxisymmetricStokes> *points) {
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
    printError("AGM::pointAxisymmetricStokes::approximateNaNDerivatives", "findInnerPointOfBoundary");
    return {};
  };
  if (std::isnan(values["dx"])) values["dx"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dx"];
  if (std::isnan(values["dy"])) values["dy"] = points->at(findInnerPointOfBoundary()->getIdx()).getValue()["dy"];
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokes(char rz) {
  switch (condition) {
    case 'C':
      makeRepresentationFormulaStokesCross(rz);
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      makeRepresentationFormulaStokesBoundary(rz);
      break;
    case 'I':
      makeRepresentationFormulaStokesInterface(rz);
      break;
    default:
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokes", "Argument `condition` is wrong, got %c.", condition);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCross(char rz) {
  if (rz == 'r') {
    makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion0();
  } else if (rz == 'z') {
    makeRepresentationFormulaStokesCrossNonSymmetricZ0();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCross", "Argument `rz` must be `r` or `z`, got %c.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion0() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetricStokes(xm, xb, xp, getMp(), getMp())};

  pMatrixRow[0][getElement()[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau(xm);
  pMatrixRow[0][getIdx()] = UNITVALUE / xb;
  pMatrixRow[0][getElement()[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau(xp);

  pMatrixRow[0][getElement()[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[0][getElement()[W]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[0][getElement()[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd()) + gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + gFunc.green_integral_tau('c', getOrd());
  //  pMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd()) + gFunc.green_integral_tau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossNonSymmetricZ0() {
  auto ym{getElement()[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};

  pMatrixRow[1][getElement()[S]->getIdx()] = getMp() * getXy()[0] * gFunc.green_function_ttau(ym);
  pMatrixRow[1][getElement()[N]->getIdx()] = -getMp() * getXy()[0] * gFunc.green_function_ttau(yp);

  pMatrixRow[1][getElement()[S]->getIdx() + getNPts()] = -gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[1][getElement()[N]->getIdx() + getNPts()] = -gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[1][getElement()[S]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[1][getElement()[N]->getIdx() + 3 * getNPts()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[1][getElement()[S]->getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('l', getOrd());
  pMatrixRow[1][getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('c', getOrd());
  pMatrixRow[1][getElement()[N]->getIdx() + 5 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesBoundary(char rz) {
  for (const auto &item : getSolMatrixRow()[1]) {
    if (getNPts() <= item.idx && item.idx < 2 * getNPts()) {
      pMatrixRow[0][item.idx + 3 * getNPts()] = item.value;
    } else {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesBoundary", "Solution matrix index must be between nPts and 2 * nPts.");
    }
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterface(char rz) {
  if (rz == 'r') {
    makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion0();
  } else if (rz == 'z') {
    makeRepresentationFormulaStokesInterfaceNonSymmetricZ0();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterface", "Argument `rz` must be `r` or `z`, got %c.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion0() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion0", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc{GreenfunctionAxisymmetricStokes(xm, xb, xp, mpw, mpe)};
  auto checkInterface = [&](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
  auto checkMatrixRow = [&checkInterface](matrixRow *row, point *ptr, point *ptl) -> void {
    if (checkInterface(ptr)) {
      (*row)[ptl->getIdx() + getNPts()] += (*row)[ptr->getIdx() + getNPts()];
      row->remove(ptr->getIdx() + getNPts());
    } else if (checkInterface(ptl)) {
      (*row)[ptr->getIdx() + getNPts()] += (*row)[ptl->getIdx() + getNPts()];
      row->remove(ptl->getIdx() + getNPts());
    }
  };
  auto approximateSol = [&](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetricStokes *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, getOrd()) : func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd()), 4);
    }
  };
  row[getIdx()] = UNITVALUE / getXy()[0];
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  pMatrixRow[0] = row;
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricZ0() {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricZ0", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
  auto checkMatrixRow = [&](matrixRow *row, point *ptr, point *ptl) -> void {
    if (checkInterface(ptr)) {
      (*row)[ptl->getIdx() + getNPts()] += (*row)[ptr->getIdx() + getNPts()];
      row->remove(ptr->getIdx() + getNPts());
    } else if (checkInterface(ptl)) {
      (*row)[ptr->getIdx() + getNPts()] += (*row)[ptl->getIdx() + getNPts()];
      row->remove(ptl->getIdx() + getNPts());
    }
  };
  auto approximateSol = [&](point *ptr, point *ptl, double coefficient, double d) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, getOrd())
                                                   : -func->green_integral(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, getOrd()), 1)
                         : linearApproximation(ptr, ptl, -func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', getOrd());
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  pMatrixRow[1] = row;
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesFull(char rz) {
  switch (condition) {
    case 'C':
      makeRepresentationFormulaStokesCrossFull(rz);
      break;
    case 'D':
    case 'd':
    case 'N':
    case 'n':
      makeRepresentationFormulaStokesBoundaryFull();
      break;
    case 'I':
      makeRepresentationFormulaStokesInterfaceFull(rz);
      break;
    default:
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesFull", "Argument `condition` is wrong, got %c.", condition);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossFull(char rz) {
  if (rz == 'r') {
    makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion();
    //    makeRepresentationFormulaStokesCrossNonSymmetricR();
  } else if (rz == 'z') {
    auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
    //    if (iszero(xm)) {
    //      makeRepresentationFormulaStokesCrossSymmetricNearAxis();
    //    } else {
    //      makeRepresentationFormulaStokesCrossSymmetricElliptic();
    //    }
    makeRepresentationFormulaStokesCrossNonSymmetricZ();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossFull", "Argument `rz` must be `r` or `z`, got %c.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossSymmetricElliptic() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  pMatrixRow[0][getElement()[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau(xm);
  pMatrixRow[0][getElement()[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau(xp);

  pMatrixRow[0][getElement()[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  prhsMatrixRow[0][getElement()[W]->getIdx()] = gFunc.green_integral_tau('l', getOrd());
  prhsMatrixRow[0][getIdx()] = gFunc.green_integral_tau('c', getOrd());
  prhsMatrixRow[0][getElement()[E]->getIdx()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[0][getElement()[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd());
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / xb;
  pMatrixRow[0][getElement()[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossSymmetricReactionDiffusion() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetricStokes(xm, xb, xp, getMp(), getMp())};
  pMatrixRow[0][getElement()[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau(xm);
  pMatrixRow[0][getIdx()] = UNITVALUE / xb;
  pMatrixRow[0][getElement()[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau(xp);

  pMatrixRow[0][getElement()[W]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  prhsMatrixRow[0][getElement()[W]->getIdx()] = gFunc.green_integral_tau('l', getOrd());
  prhsMatrixRow[0][getIdx()] = gFunc.green_integral_tau('c', getOrd());
  prhsMatrixRow[0][getElement()[E]->getIdx()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[0][getElement()[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd()) + gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd()) + gFunc.green_integral_tau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossSymmetricNearAxis() {
  auto xm{getElement()[W]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E]->getXy()[0]};
  auto gFunc{GreenfunctionAxisymmetric(xm, xb, xp, getMp(), getMp())};
  pMatrixRow[0][getElement()[W]->getIdx()] = getMp() * xm * gFunc.green_function_ttau_ND(xm);
  pMatrixRow[0][getElement()[E]->getIdx()] = -getMp() * xp * gFunc.green_function_ttau_ND(xp);

  pMatrixRow[0][getElement()[W]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('l', getOrd());
  pMatrixRow[0][getIdx() + getNPts()] = gFunc.green_integral_tau_ND('c', getOrd());
  pMatrixRow[0][getElement()[E]->getIdx() + getNPts()] = gFunc.green_integral_tau_ND('r', getOrd());

  prhsMatrixRow[0][getElement()[W]->getIdx()] = gFunc.green_integral_tau_ND('l', getOrd());
  prhsMatrixRow[0][getIdx()] = gFunc.green_integral_tau_ND('c', getOrd());
  prhsMatrixRow[0][getElement()[E]->getIdx()] = gFunc.green_integral_tau_ND('r', getOrd());

  pMatrixRow[0][getElement()[W]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('l', getOrd());
  pMatrixRow[0][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('c', getOrd()) + UNITVALUE / getMp() / xb;
  pMatrixRow[0][getElement()[E]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau_ND('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossNonSymmetricR() {
  auto ym{getElement()[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  pMatrixRow[1][getElement()[S]->getIdx() + 2 * getNPts()] = getMp() * getXy()[0] * gFunc.green_function_ttau(ym);
  pMatrixRow[1][getElement()[N]->getIdx() + 2 * getNPts()] = -getMp() * getXy()[0] * gFunc.green_function_ttau(yp);

  pMatrixRow[1][getElement()[S]->getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[1][getElement()[N]->getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('r', getOrd());

  prhsMatrixRow[1][getElement()[S]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  prhsMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  prhsMatrixRow[1][getElement()[N]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[1][getElement()[S]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('l', getOrd());
  pMatrixRow[1][getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  pMatrixRow[1][getElement()[N]->getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesCrossNonSymmetricZ() {
  auto ym{getElement()[S]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N]->getXy()[1]};
  auto gFunc{Greenfunction(ym, yb, yp, getMp() * getXy()[0], getMp() * getXy()[0])};
  pMatrixRow[1][getElement()[S]->getIdx() + 2 * getNPts()] = getMp() * getXy()[0] * gFunc.green_function_ttau(ym);
  pMatrixRow[1][getElement()[N]->getIdx() + 2 * getNPts()] = -getMp() * getXy()[0] * gFunc.green_function_ttau(yp);

  pMatrixRow[1][getElement()[S]->getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('l', getOrd());
  pMatrixRow[1][getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('c', getOrd());
  pMatrixRow[1][getElement()[N]->getIdx() + 3 * getNPts()] = -gFunc.green_integral_tau('r', getOrd());

  prhsMatrixRow[1][getElement()[S]->getIdx() + getNPts()] = gFunc.green_integral_tau('l', getOrd());
  prhsMatrixRow[1][getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  prhsMatrixRow[1][getElement()[N]->getIdx() + getNPts()] = gFunc.green_integral_tau('r', getOrd());

  pMatrixRow[1][getElement()[S]->getIdx() + 4 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('l', getOrd());
  pMatrixRow[1][getIdx() + 4 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp();
  pMatrixRow[1][getElement()[N]->getIdx() + 4 * getNPts()] = getXy()[0] * gFunc.green_integral_ttau('r', getOrd());
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceFull(char rz) {
  if (rz == 'r') {
    makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion();
    makeRepresentationFormulaStokesInterfaceNonSymmetricR();
  } else if (rz == 'z') {
    auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
    if (iszero(xm)) {
      makeRepresentationFormulaStokesInterfaceSymmetricNearAxis();
    } else {
      makeRepresentationFormulaStokesInterfaceSymmetricElliptic();
    }
    makeRepresentationFormulaStokesInterfaceNonSymmetricZ();
  } else {
    printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceFull", "input argument rz (which is %c) must be 'r' or 'z'.", rz);
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricElliptic() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricElliptic", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, getOrd()) : func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()), 4);
    }
  };
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() / getXy()[0];
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  if (!row.empty()) {
    for (auto &item : row) {
      if (item.idx < 2 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
        prhsMatrixRow[0][item.idx - 2 * getNPts()] = item.value;
      } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (5 * getNPts() <= item.idx) {
        pMatrixRow[0][item.idx - getNPts()] += item.value;
      }
    }
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricReactionDiffusion", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetricStokes(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetricStokes *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_tau(C, getOrd()) : func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_tau(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_tau(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_tau(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_ttau(c, getOrd()) + func->green_integral_tau(c, getOrd()), 4);
    }
  };
  row[getIdx()] = UNITVALUE / getXy()[0];
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_tau('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_ttau('c', getOrd()) + UNITVALUE / getMp() + gFunc.green_integral_tau('c', getOrd());
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  if (!row.empty()) {
    for (auto &item : row) {
      if (item.idx < 2 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
        prhsMatrixRow[0][item.idx - 2 * getNPts()] = item.value;
      } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (5 * getNPts() <= item.idx) {
        pMatrixRow[0][item.idx - getNPts()] += item.value;
      }
    }
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricNearAxis() {
  auto xm{getElement()[W] ? getElement()[W]->getXy()[0] : getElement()[WN]->getXy()[0]};
  auto xb{getXy()[0]};
  auto xp{getElement()[E] ? getElement()[E]->getXy()[0] : getElement()[EN]->getXy()[0]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceSymmetricNearAxis", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpe{getEachMp(getElement()[E], getElement()[EN], getElement()[ES])};
  auto mpw{getEachMp(getElement()[W], getElement()[WN], getElement()[WS])};
  auto gFunc(GreenfunctionAxisymmetric(xm, xb, xp, mpw, mpe));
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};
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
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto func{Greenfunction(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = -func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = -func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 3 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 3 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 5 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 5 * getNPts()] = func.green_integral_t('R', getOrd());

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[1]}, b{getXy()[1]}, p{ptr->getXy()[1]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, GreenfunctionAxisymmetric *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * func->green_function_t_ND(d);
      row[pt->getIdx() + getNPts()] += isInterface ? func->green_integral_ND(C, getOrd()) : func->green_integral_ND(c, getOrd());
      row[pt->getIdx() + 2 * getNPts()] += func->green_integral_ND(c, getOrd());
      row[pt->getIdx() + 4 * getNPts()] += func->green_integral_t_ND(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * func->green_function_t_ND(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, func->green_integral_ND(C, getOrd()), 1) : linearApproximation(ptr, ptl, func->green_integral_ND(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral_ND(c, getOrd()), 2);
      row += linearApproximation(ptr, ptl, func->green_integral_t_ND(c, getOrd()), 4);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = gFunc.green_integral_ND('c', getOrd());
  row[getIdx() + 2 * getNPts()] = gFunc.green_integral_ND('c', getOrd());
  row[getIdx() + 4 * getNPts()] = gFunc.green_integral_t_ND('c', getOrd());
  assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe * xp, &gFunc, xp, 'r', 'R');
  assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw * xm, &gFunc, xm, 'l', 'L');

  if (!row.empty()) {
    for (auto &item : row) {
      if (item.idx < 2 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
        prhsMatrixRow[0][item.idx - 2 * getNPts()] = item.value;
      } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
        pMatrixRow[0][item.idx] = item.value;
      } else if (5 * getNPts() <= item.idx) {
        pMatrixRow[0][item.idx - getNPts()] += item.value;
      }
    }
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricR() {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricR", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, getOrd()) : -func->green_integral(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, getOrd()), 1) : linearApproximation(ptr, ptl, -func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', getOrd());
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  if (!row.empty()) {
    for (auto &item : row) {
      if (item.idx < 2 * getNPts()) {
        pMatrixRow[1][item.idx] = item.value;
      } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
        prhsMatrixRow[1][item.idx - 2 * getNPts()] = item.value;
      } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
        pMatrixRow[1][item.idx] = item.value;
      } else if (5 * getNPts() <= item.idx) {
        pMatrixRow[1][item.idx - getNPts()] += item.value;
      }
    }
  }
}

void AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricZ() {
  auto ym{getElement()[S] ? getElement()[S]->getXy()[1] : getElement()[SE]->getXy()[1]};
  auto yb{getXy()[1]};
  auto yp{getElement()[N] ? getElement()[N]->getXy()[1] : getElement()[NE]->getXy()[1]};
  auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
    auto Error = []() -> double {
      printError("AGM::pointAxisymmetricStokes::makeRepresentationFormulaStokesInterfaceNonSymmetricZ", "getEachMp");
      return {};
    };
    auto rtv{pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp()
                 : ptl->getCondition() == 'C'              ? ptl->getMp()
                                                           : Error()};
    return rtv;
  };
  auto mpn{getEachMp(getElement()[N], getElement()[NE], getElement()[NW])};
  auto mps{getEachMp(getElement()[S], getElement()[SE], getElement()[SW])};
  auto gFunc{Greenfunction(ym, yb, yp, mps * getXy()[0], mpn * getXy()[0])};
  auto checkInterface = [&getEachMp](point *pt) -> bool {
    auto mpe{getEachMp(pt->getElement()[E], pt->getElement()[EN], pt->getElement()[ES])};
    auto mpw{getEachMp(pt->getElement()[W], pt->getElement()[WN], pt->getElement()[WS])};
    auto mpn{getEachMp(pt->getElement()[N], pt->getElement()[NE], pt->getElement()[NW])};
    auto mps{getEachMp(pt->getElement()[S], pt->getElement()[SE], pt->getElement()[SW])};
    return !(isclose(mpe, mpw) && isclose(mpn, mps) && isclose(mpe, mps));
  };
  auto isInterface{checkInterface(this)};

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
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto func{GreenfunctionAxisymmetric(m, b, p, d, d)};
    auto mRow{matrixRow()};
    mRow[ptl->getIdx()] = d * m * func.green_function_t(m);
    mRow[ptr->getIdx()] = -d * p * func.green_function_t(p);

    mRow[ptl->getIdx() + getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + getNPts()] = func.green_integral('R', getOrd());

    checkMatrixRow(&mRow, ptr, ptl);

    mRow[ptl->getIdx() + 2 * getNPts()] = func.green_integral('L', getOrd());
    mRow[ptr->getIdx() + 2 * getNPts()] = func.green_integral('R', getOrd());

    mRow[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t('L', getOrd());
    mRow[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t('R', getOrd());

    if (iszero(m)) {
      auto mRow0{matrixRow()};
      mRow0[ptl->getIdx()] = d * m * func.green_function_t_ND(m);
      mRow0[ptr->getIdx()] = -d * p * func.green_function_t_ND(p);

      mRow0[ptl->getIdx() + getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + getNPts()] = func.green_integral_ND('R', getOrd());

      checkMatrixRow(&mRow0, ptr, ptl);

      mRow0[ptl->getIdx() + 2 * getNPts()] = func.green_integral_ND('L', getOrd());
      mRow0[ptr->getIdx() + 2 * getNPts()] = func.green_integral_ND('R', getOrd());

      mRow0[ptl->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('L', getOrd());
      mRow0[ptr->getIdx() + 4 * getNPts()] = func.green_integral_t_ND('R', getOrd());
      return mRow0 * coefficient;
    }

    return mRow * coefficient;
  };
  auto linearApproximation = [this](point *ptr, point *ptl, double coefficient, int plus) -> matrixRow {
    auto m{ptl->getXy()[0]}, b{getXy()[0]}, p{ptr->getXy()[0]};
    auto mRow = matrixRow();
    mRow[ptl->getIdx() + plus * getNPts()] = (p - b) / (p - m);
    mRow[ptr->getIdx() + plus * getNPts()] = (b - m) / (p - m);

    return mRow * coefficient;
  };
  auto row{matrixRow()};
  auto assignMatrix = [&](point *pt, point *ptr, point *ptl, double mp0, Greenfunction *func, double d, char c, char C) -> void {
    if (pt) {
      row[pt->getIdx()] += mp0 * getXy()[0] * func->green_function_t(d);
      row[pt->getIdx() + getNPts()] += isInterface ? -func->green_integral(C, getOrd())
                                                   : -func->green_integral(c, getOrd());
      row[pt->getIdx() + 3 * getNPts()] += func->green_integral(c, getOrd());
      row[pt->getIdx() + 5 * getNPts()] += func->green_integral_t(c, getOrd());
    } else {
      row += approximateSol(ptr, ptl, mp0 * getXy()[0] * func->green_function_t(d), std::abs(mp0));
      row += isInterface ? linearApproximation(ptr, ptl, -func->green_integral(C, getOrd()), 1)
                         : linearApproximation(ptr, ptl, -func->green_integral(c, getOrd()), 1);
      row += linearApproximation(ptr, ptl, func->green_integral(c, getOrd()), 3);
      row += linearApproximation(ptr, ptl, func->green_integral_t(c, getOrd()), 5);
    }
  };
  row[getIdx()] = -UNITVALUE;
  if (!isInterface) row[getIdx() + getNPts()] = -gFunc.green_integral('c', getOrd());
  row[getIdx() + 3 * getNPts()] = gFunc.green_integral('c', getOrd());
  row[getIdx() + 5 * getNPts()] = gFunc.green_integral_t('c', getOrd());
  assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFunc, yp, 'r', 'R');
  assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFunc, ym, 'l', 'L');

  if (!row.empty()) {
    for (auto &item : row) {
      if (item.idx < 2 * getNPts()) {
        pMatrixRow[1][item.idx] = item.value;
      } else if (2 * getNPts() <= item.idx && item.idx < 4 * getNPts()) {
        prhsMatrixRow[1][item.idx - 2 * getNPts()] = item.value;
      } else if (4 * getNPts() <= item.idx && item.idx < 5 * getNPts()) {
        pMatrixRow[1][item.idx] = item.value;
      } else if (5 * getNPts() <= item.idx) {
        pMatrixRow[1][item.idx - getNPts()] += item.value * getXy()[0];
      }
    }
  }
}
