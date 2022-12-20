//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include <vector>
#include "point.h"

int AGM::point::nPts;
std::vector<AGM::axialLine> *AGM::point::xline;
std::vector<AGM::axialLine> *AGM::point::yline;
std::vector<AGM::boundaryLine2D> *AGM::point::bdLine;
std::vector<AGM::point> *AGM::point::pts;

AGM::point::point() = default;

AGM::point::point(int idx) : idx(idx) {}

AGM::point::point(const AGM::coordinate &xy) : xy(xy) {}

AGM::point::point(int idx, const AGM::coordinate &xy) : idx(idx), xy(xy) {}

AGM::point::point(const AGM::coordinate &xy, double mp) : xy(xy), mp(mp) {}

AGM::point::point(int idx, const AGM::coordinate &xy, double mp) : idx(idx), xy(xy), mp(mp) {}

auto AGM::point::getIdx() const -> int {
    return idx;
}

void AGM::point::setIdx(int i) {
    point::idx = i;
}

auto AGM::point::getXy() const -> const AGM::coordinate & {
    return xy;
}

void AGM::point::setXy(const AGM::coordinate &coordinate) {
    point::xy = coordinate;
}

auto AGM::point::getNormal() const -> const AGM::coordinate & {
    return normal;
}

void AGM::point::setNormal(const AGM::coordinate &coordinate) {
    point::normal = coordinate;
}

auto AGM::point::getMp() const -> double {
    return mp;
}

void AGM::point::setMp(double d) {
    point::mp = d;
}

auto AGM::point::getCondition() const -> char {
    return condition;
}

void AGM::point::setCondition(char i) {
    point::condition = i;
}

auto AGM::point::getElement() const -> const std::array<AGM::point *, 12> & {
    return element;
}

void AGM::point::setElement(const std::array<AGM::point *, 12> &array) {
    point::element = array;
}

auto AGM::point::getValue() const -> const AGM::value & {
    return values;
}

void AGM::point::setValue(const AGM::value &value) {
    point::values = value;
}

auto AGM::point::getSolMatrixRow() const -> const std::array<AGM::matrixRow, 2> & {
    return solMatrixRow;
}

void AGM::point::setSolMatrixRow(const std::array<AGM::matrixRow, 2> &row) {
    point::solMatrixRow = row;
}

auto AGM::point::getDeriMatrixRow() const -> const std::array<AGM::matrixRow, 2> & {
    return deriMatrixRow;
}

void AGM::point::setDeriMatrixRow(const std::array<AGM::matrixRow, 2> &row) {
    point::deriMatrixRow = row;
}

auto AGM::point::getRb() const -> const std::array<double, 2> & {
    return rb;
}

void AGM::point::setRb(const std::array<double, 2> &array) {
    point::rb = array;
}

auto AGM::point::getAxialLine() const -> const std::array<AGM::axialLine *, 2> & {
    return aline;
}

auto AGM::point::getAxialLine(char i) -> AGM::axialLine *& {
    if (i == 'x') return aline[0];
    if (i == 'y') return aline[1];
    if (i == 'z') return aline[2];
    printError("AGM::axialLine *&AGM::point::getAxialLine", "index (which = %c) is wrong", i);
    return aline[0];
}

void AGM::point::setAxialLine(const std::array<AGM::axialLine *, 2> &array) {
    point::aline = array;
}

void AGM::point::setAxialLine(AGM::axialLine *line, char i) {
    if (i == 'x') aline[0] = line;
    if (i == 'y') aline[1] = line;
    if (i == 'z') aline[2] = line;
}

auto AGM::point::getNPts() -> int {
    return nPts;
}

void AGM::point::setNPts(int i) {
    point::nPts = i;
}

auto AGM::point::getAxialLines(char i) -> std::vector<AGM::axialLine> *& {
    if (i == 'x') return xline;
    if (i == 'y') return yline;
    printError("std::vector<AGM::axialLine> *&AGM::point::getAxialLines", "input character (which is %c) is wrong", i);
    return xline;
}

void AGM::point::setAxialLines(std::vector<AGM::axialLine> *line, char i) {
    if (i == 'x') xline = line;
    if (i == 'y') yline = line;
}

auto AGM::point::getBdLine() -> std::vector<AGM::boundaryLine2D> * {
    return bdLine;
}

void AGM::point::setBdLine(std::vector<AGM::boundaryLine2D> *line) {
    point::bdLine = line;
}

auto AGM::point::getPts() -> std::vector<AGM::point> * {
    return pts;
}

void AGM::point::setPts(std::vector<AGM::point> *vector) {
    point::pts = vector;
}

auto AGM::point::operator[](int i) -> double & {
    return xy[i];
}

auto AGM::point::operator[](int i) const -> const double & {
    return xy[i];
}

auto AGM::point::operator[](AGM::EWNS ewns) -> AGM::point *& {
    return element[ewns];
}

auto AGM::point::operator[](const std::string &string) -> double & {
    return values[string];
}

auto AGM::point::operator[](const std::string &string) const -> const double & {
    return values[string];
}

auto AGM::point::operator-(const AGM::point &src) -> double {
    return (xy - src.xy).norm();
}

auto AGM::point::operator-(axialLine &src) -> double {
    if (src.getMark() == 'x') {
        return xy[1] - src[2];
    } else if (src.getMark() == 'y') {
        return xy[0] - src[0];
    } else {
        printError("AGM::point::operator-", "axialLine mark (which is %c) error", src.getMark());
    }
    return ZEROVALUE;
}

auto AGM::point::operator=(const AGM::point &src) -> AGM::point & {
    if (this != &src) {
        idx = src.idx;
        xy = src.xy;
        normal = src.normal;
        mp = src.mp;
        condition = src.condition;
        element = src.element;
        values = src.values;
        solMatrixRow = src.solMatrixRow;
        deriMatrixRow = src.deriMatrixRow;
        rb = src.rb;
        dv = src.dv;
        aline = src.aline;
    }
    return *this;
}

void AGM::point::findStencil() {
    switch (condition) {
        case 'D':
        case 'N':
            findStencilBoundary();
            break;
        case 'd':
        case 'n':
            findStencilAppendBoundary();
            break;
        case 'I':
            findStencilInterface();
            break;
        default:
            break;
    }
}

void AGM::point::findStencilBoundary() {
    auto isContainLine = [](axialLine *aln, int lineIdx, double pt) -> bool {
        return isnegativezero((*aln)[lineIdx] - pt) && ispositivezero((*aln)[lineIdx + 1] - pt);
    };
    auto findLeftLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt)) {
                vec.emplace_back(&item);
            }
        }
        auto pline = std::find_if(vec.begin(), vec.end(), [this, &lineChar](axialLine *&aln) -> bool {
            return *std::next(&aln) == getAxialLine(lineChar);
        });
        auto vec0 = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (ispositive(*getAxialLine(lineChar) - item)) {
                vec0.emplace_back(&item);
            }
        }
        auto pline0 = vec0.empty() ? vec0.end() : std::prev(vec0.end());
        if (pline == vec.end() || pline0 == vec0.end() || !iszero(**pline - **pline0)) return nullptr;
        return *pline;
    };
    auto findRightLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt)) {
                vec.emplace_back(&item);
            }
        }
        auto pline = std::find_if(vec.begin(), vec.end(), [this, &lineChar](axialLine *&aln) -> bool {
            return *std::prev(&aln) == getAxialLine(lineChar);
        });
        auto vec0 = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isnegative(*getAxialLine(lineChar) - item)) {
                vec0.emplace_back(&item);
            }
        }
        auto pline0 = vec0.empty() ? vec0.end() : vec0.begin();
        if (pline == vec.end() || pline0 == vec0.end() || !iszero(**pline - **pline0)) return nullptr;
        return *pline;
    };
    auto findLeftPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return isnegativezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.back();
    };
    auto findRightPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return ispositivezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.front();
    };
    auto assignStencil = [&](EWNS ewns, EWNS ewns0, EWNS ewns1, auto func, char lineChar, int lineIdx, double pt,
                             double sign, double &n) -> void {
        if (element[ewns] == this && ispositive(sign * n)) return;
        auto alin = func(lineChar, lineIdx, pt);
        if (alin) {
            auto leftPt = findLeftPt(alin, lineIdx);
            auto rightPt = findRightPt(alin, lineIdx);
            if (leftPt && rightPt && leftPt == rightPt) {
                element[ewns] = leftPt;
                leftPt = nullptr;
                rightPt = nullptr;
                return;
            }
            if (rightPt) {
                element[ewns0] = rightPt;
                element[ewns] = nullptr;
            }
            if (leftPt) {
                element[ewns1] = leftPt;
                element[ewns] = nullptr;
            }

            if (!(rightPt || leftPt)) {
                printError("assignStencil in findStencil", "rightPt & leftPt do not exist");
            }
        } else {
            n = ZEROVALUE;
            double norm{std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2))};
            for (int i = 0; i < 2; ++i) {
                normal[i] /= norm;
            }
            element[ewns] = this;
        }
    };
    if (!element[E]) assignStencil(E, EN, ES, findRightLine, 'y', 1, xy[1], UNITVALUE, normal[0]);
    if (!element[W]) assignStencil(W, WN, WS, findLeftLine, 'y', 1, xy[1], -UNITVALUE, normal[0]);

    if (!element[N]) assignStencil(N, NE, NW, findRightLine, 'x', 0, xy[0], UNITVALUE, normal[1]);
    if (!element[S]) assignStencil(S, SE, SW, findLeftLine, 'x', 0, xy[0], -UNITVALUE, normal[1]);
}

void AGM::point::findStencilAppendBoundary() {
    auto isContainLine = [](axialLine *aln, int lineIdx, double pt) -> bool {
        return isnegativezero((*aln)[lineIdx] - pt) && ispositivezero((*aln)[lineIdx + 1] - pt);
    };
    auto findLeftLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt) && ispositive(*this - item)) {
                vec.emplace_back(&item);
            }
        }
        if (vec.empty()) return nullptr;
        std::sort(vec.begin(), vec.end(), [](axialLine *a, axialLine *b) -> bool {
            return ispositivezero(*b - *a);
        });
        return vec.back();
    };
    auto findRightLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt) && isnegative(*this - item)) {
                vec.emplace_back(&item);
            }
        }
        if (vec.empty()) return nullptr;
        std::sort(vec.begin(), vec.end(), [](axialLine *a, axialLine *b) -> bool {
            return ispositivezero(*b - *a);
        });
        return vec.front();
    };
    auto findLeftPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return isnegativezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.back();
    };
    auto findRightPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return ispositivezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.front();
    };
    auto assignStencil = [&](EWNS ewns, EWNS ewns0, EWNS ewns1, auto func, char lineChar, int lineIdx, double pt,
                             double sign, double &n) -> void {
        if (element[ewns] == this && ispositive(sign * n)) return;
        auto alin = func(lineChar, lineIdx, pt);
        if (alin) {
            auto leftPt = findLeftPt(alin, lineIdx);
            auto rightPt = findRightPt(alin, lineIdx);
            if (leftPt && rightPt && leftPt == rightPt) {
                element[ewns] = leftPt;
                leftPt = nullptr;
                rightPt = nullptr;
                return;
            }
            if (rightPt) {
                element[ewns0] = rightPt;
                element[ewns] = nullptr;
            }
            if (leftPt) {
                element[ewns1] = leftPt;
                element[ewns] = nullptr;
            }

            if (!(rightPt || leftPt)) {
                printError("assignStencil in findStencil", "rightPt & leftPt do not exist");
            }
        } else {
            n = ZEROVALUE;
            double norm{std::sqrt(std::pow(normal[0], 2) + std::pow(normal[1], 2))};
            for (int i = 0; i < 2; ++i) {
                normal[i] /= norm;
            }
            element[ewns] = this;
        }
    };
    if (!element[E]) assignStencil(E, EN, ES, findRightLine, 'y', 1, xy[1], UNITVALUE, normal[0]);
    if (!element[W]) assignStencil(W, WN, WS, findLeftLine, 'y', 1, xy[1], -UNITVALUE, normal[0]);

    if (!element[N]) assignStencil(N, NE, NW, findRightLine, 'x', 0, xy[0], UNITVALUE, normal[1]);
    if (!element[S]) assignStencil(S, SE, SW, findLeftLine, 'x', 0, xy[0], -UNITVALUE, normal[1]);
}

void AGM::point::findStencilInterface() {
    auto assignDist = [this]() -> double {
        double vv{};
        for (const auto &item: {E, W, N, S}) {
            if (getElement()[item]) {
                if (5e0 * (*getElement()[item] - *this) > vv) {
                    vv = 5e0 * (*getElement()[item] - *this);
                }
            }
        }
        if (!iszero(vv)) {
            return vv;
        }
        printError("AGM::point::findStencilInterface, assignDist");
        return ZEROVALUE;
    };
    double dist{assignDist()};
    auto isContainLine = [](axialLine *aln, int lineIdx, double pt) -> bool {
        return isnegativezero((*aln)[lineIdx] - pt) && ispositivezero((*aln)[lineIdx + 1] - pt);
    };
    auto findLeftLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt)) {
                vec.emplace_back(&item);
            }
        }
        std::sort(vec.begin(), vec.end(), [](axialLine *&a, axialLine *&b) -> bool {
            return isnegative(*a - *b);
        });
        auto vec0 = std::vector<axialLine *>{};
        for (auto &item: vec) {
            if (isnegative((*item)[(2 * lineIdx + 2) % 4] - xy[(lineIdx + 1) % 2])) {
                vec0.emplace_back(item);
            }
        }
        if (vec0.empty()) return nullptr;
        return vec0.back();
    };
    auto findRightLine = [this, &isContainLine](char lineChar, int lineIdx, double pt) -> axialLine * {
        auto vec = std::vector<axialLine *>{};
        for (auto &item: *getAxialLines(lineChar)) {
            if (isContainLine(&item, 2 * lineIdx, pt)) {
                vec.emplace_back(&item);
            }
        }
        std::sort(vec.begin(), vec.end(), [](axialLine *&a, axialLine *&b) -> bool {
            return isnegative(*a - *b);
        });
        auto vec0 = std::vector<axialLine *>{};
        for (auto &item: vec) {
            if (ispositive((*item)[(2 * lineIdx + 2) % 4] - xy[(lineIdx + 1) % 2])) {
                vec0.emplace_back(item);
            }
        }
        if (vec0.empty()) return nullptr;
        return vec0.front();
    };
    auto findLeftPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return isnegativezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.back();
    };
    auto findRightPt = [this](axialLine *aln, int ptIdx) -> point * {
        auto vec = std::vector<point *>{};
        std::copy_if(aln->begin(), aln->end(), std::back_inserter(vec), [this, &ptIdx](point *pt) -> bool {
            return ispositivezero((*pt)[ptIdx] - xy[ptIdx]);
        });
        return vec.front();
    };
    auto lineEndPt = [this, &dist](EWNS ewns) -> std::vector<double> {
        if (ewns == E) {
            return std::vector<double>{getXy()[0] + dist, getXy()[1]};
        } else if (ewns == W) {
            return std::vector<double>{getXy()[0] - dist, getXy()[1]};
        } else if (ewns == N) {
            return std::vector<double>{getXy()[0], getXy()[1] + dist};
        } else if (ewns == S) {
            return std::vector<double>{getXy()[0], getXy()[1] - dist};
        }
        printError("AGM::point::findStencilInterface, lineEndPt");
        return std::vector<double>{};
    };
    auto assignAppendElement = [this](point &pt, EWNS ewns, const coordinate &coordinate) -> void {
        auto func = [this, &pt, &coordinate](EWNS ewns0, EWNS ewns1, EWNS ewns2, EWNS ewns3, int i) -> void {
            pt[ewns1] = this;
            pt[ewns0] = &pt;
            if (ispositive(coordinate[i])) pt[ewns2] = &pt;
            else if (isnegative(coordinate[i])) pt[ewns3] = &pt;
            else pt[ewns2] = pt[ewns3] = &pt;
        };
        if (ewns == E) func(ewns, W, N, S, 1);
        else if (ewns == W) func(ewns, E, N, S, 1);
        else if (ewns == N) func(ewns, S, E, W, 0);
        else if (ewns == S) func(ewns, N, E, W, 0);
        else printError("AGM::point::findStencilInterface", "assignAppendElement");
    };
    auto assignStencil = [&](EWNS ewns, EWNS ewns0, EWNS ewns1, auto func, char lineChar, int lineIdx,
                             double pt) -> void {
        auto alin = func(lineChar, lineIdx, pt);
        if (alin && std::fabs((*this) - *alin) > dist) {
            std::cout << "(xx, yy) = " << xy[0] << ", " << xy[1] << "\n";
            alin = nullptr;
        }
        if (alin) {
            auto leftPt = findLeftPt(alin, lineIdx);
            auto rightPt = findRightPt(alin, lineIdx);
            if (leftPt && rightPt && leftPt == rightPt) {
                element[ewns] = leftPt;
                leftPt = nullptr;
                rightPt = nullptr;
                return;
            }
            if (rightPt) {
                element[ewns0] = rightPt;
                element[ewns] = nullptr;
            }
            if (leftPt) {
                element[ewns1] = leftPt;
                element[ewns] = nullptr;
            }

            if (!(rightPt || leftPt)) {
                printError("assignStencil in findStencil", "rightPt & leftPt do not exist");
            }
        } else {
            auto line{line2D{}};
            auto start{vector{std::vector<double>{getXy()[0], getXy()[1]}}};
            auto end{vector{lineEndPt(ewns)}};
            auto vec{vector{}};
            vec.resize(2);
            auto crossLine{std::vector<point>{}};
            line.setStart(start);
            line.setAnEnd(end);
            line.calcProperties();
            for (const auto &item: *(getBdLine())) {
                if (line.iscross(item, vec) && item.getCondition() != 'I') {
                    crossLine.emplace_back(point(coordinate(vec[0], vec[1]), getMp()));
                    crossLine.back().setCondition(char(std::tolower(item.getCondition())));
                    crossLine.back().setNormal(coordinate(item.getNormal()[0], item.getNormal()[1]));
                    crossLine.back()["bdv"] = item.getBoundaryValue();
                }
            }
            std::sort(crossLine.begin(), crossLine.end(), [this](point &a, point &b) -> bool {
                return (a - *this) < (b - *this);
            });
            if (crossLine.size() > 1) {
                for (const auto &item: crossLine) {
                    std::cout << "(x, y) = " << item[0] << ", " << item[1] << "\n";
                }
//                printError("AGM::point::findStencilInterface",
//                           "The size of crossLine (which is %d) more than 1, should be corrected later. (using sort algorithm)",
//                           crossLine.size());
            } else if (crossLine.empty()) {
                std::cout << "(x, y) = (" << getXy()[0] << ", " << getXy()[1] << ")\n";
                printError("AGM::point::findStencilInterface",
                           "The size of crossLine is empty");
            }
            getPts()->emplace_back(crossLine[0]);
            getPts()->back().setIdx(int(getPts()->size()) - 1);
            assignAppendElement(getPts()->back(), ewns, getPts()->back().getNormal());
            element[ewns] = &(getPts()->back());

            std::cout << "\nAppend point\n";
            std::cout << "idx = " << getPts()->back().getIdx() << "\n";
            std::cout << "xy = (" << getPts()->back().getXy()[0] << ", " << getPts()->back().getXy()[1] << ")\n";
            std::cout << "condition = " << getPts()->back().getCondition() << "\n";
            std::cout << "boundary value = " << getPts()->back()["bdv"] << "\n";
            std::cout << "normal = (" << getPts()->back().getNormal()[0] << ", " << getPts()->back().getNormal()[1]
                      << ")\n\n";
//            std::cout << "ewns = " << ewns << "\n";
//            std::cout << "E = " << getPts()->back()[E]->getIdx() << "\n";
//            std::cout << "W = " << getPts()->back()[W]->getIdx() << "\n";
//            std::cout << "N = " << getPts()->back()[N]->getIdx() << "\n";
//            std::cout << "S = " << getPts()->back()[S]->getIdx() << "\n";
//            printError("AGM::point::findStencilInterface", "there is no axial line");
        }
    };
    if (!element[E]) assignStencil(E, EN, ES, findRightLine, 'y', 1, xy[1]);
    if (!element[W]) assignStencil(W, WN, WS, findLeftLine, 'y', 1, xy[1]);

    if (!element[N]) assignStencil(N, NE, NW, findRightLine, 'x', 0, xy[0]);
    if (!element[S]) assignStencil(S, SE, SW, findLeftLine, 'x', 0, xy[0]);
}

void AGM::point::findStencil(std::vector<AGM::point> *src, std::vector<AGM::point> *tgt) {
    int index{};
    auto elt = axialElement{};
    for (int i = 0; i < 13; ++i) {
        if (src->at(getIdx()).getElement()[i]) {
            element[i] = &(tgt->at(src->at(getIdx()).getElement()[i]->getIdx()));
        }
    }
}

void AGM::point::calculateRepresentationFormula() {
    switch (condition) {
        case 'C':
            calculateRepresentationFormulaCross();
            break;
        case 'D':
        case 'd':
            calculateRepresentationFormulaDirichlet();
            break;
        case 'N':
        case 'n':
            calculateRepresentationFormulaNeumann();
            break;
        case 'I':
            calculateRepresentationFormulaInterface();
            break;
        default:
            printError("AGM::Point::calcRepresentationFormula", "boundary condition (which is %c) is wrong", condition);
    }
}

void AGM::point::calculateRepresentationFormulaCross() {
    double xm = element[W]->getXy()[0];
    double xb = getXy()[0];
    double xp = element[E]->getXy()[0];
    double ym = element[S]->getXy()[1];
    double yb = getXy()[1];
    double yp = element[N]->getXy()[1];
    auto gfuncX{Greenfunction(xm, xb, xp, mp, mp)};
    auto gfuncY{Greenfunction(ym, yb, yp, mp, mp)};
    std::array<matrixRow, 2> row{};
    auto eraseInterface = [this, &row](point *pt, int i) -> void {
        auto checkInterface = [](point *pt) -> bool {
            auto getEachMp = [](point *pt, point *ptr, point *ptl) -> double {
                auto Error = []() -> double {
                    printError("AGM::point::calculateRepresentationFormulaCross", "getEachMp");
                    return ZEROVALUE;
                };
                double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ? ptr->getMp() : ptl->getCondition() == 'C'
                                                                                            ? ptl->getMp() : Error();
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

    row[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral('l');
    row[0][getIdx() + getNPts()] = gfuncX.green_integral('c');
    row[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral('r');

    rhsMatrixRow[0][element[W]->getIdx()] = gfuncX.green_integral('l');
    rhsMatrixRow[0][getIdx()] = gfuncX.green_integral('c');
    rhsMatrixRow[0][element[E]->getIdx()] = gfuncX.green_integral('r');

    partMatrixRow[0][element[W]->getIdx()] = gfuncX.green_integral_t('l');
    partMatrixRow[0][getIdx()] = gfuncX.green_integral_t('c');
    partMatrixRow[0][element[E]->getIdx()] = gfuncX.green_integral_t('r');

    row[1][getIdx()] = -UNITVALUE;
    row[1][element[S]->getIdx()] = mp * gfuncY.green_function_t(ym);
    row[1][element[N]->getIdx()] = -mp * gfuncY.green_function_t(yp);

    row[1][element[S]->getIdx() + getNPts()] = -gfuncY.green_integral('l');
    row[1][getIdx() + getNPts()] = -gfuncY.green_integral('c');
    row[1][element[N]->getIdx() + getNPts()] = -gfuncY.green_integral('r');

    rhsMatrixRow[1][element[S]->getIdx() + getNPts()] = gfuncY.green_integral('l');
    rhsMatrixRow[1][getIdx() + getNPts()] = gfuncY.green_integral('c');
    rhsMatrixRow[1][element[N]->getIdx() + getNPts()] = gfuncY.green_integral('r');

    partMatrixRow[1][element[S]->getIdx() + getNPts()] = gfuncY.green_integral_t('l');
    partMatrixRow[1][getIdx() + getNPts()] = gfuncY.green_integral_t('c');
    partMatrixRow[1][element[N]->getIdx() + getNPts()] = gfuncY.green_integral_t('r');

    eraseInterface(element[E], 0);
    eraseInterface(element[W], 0);
    eraseInterface(element[N], 1);
    eraseInterface(element[S], 1);

    solMatrixRow[0] = row[0] + row[1];
    solMatrixRow[1] = row[0] - row[1];
}

void AGM::point::calculateRepresentationFormulaDirichlet() {
    solMatrixRow[0][getIdx()] = UNITVALUE;
    if (getCondition() == 'D') approximatePhiAtBoundary(1);
    else if (getCondition() == 'd') approximatePhiAtAppend();
}

void AGM::point::calculateRepresentationFormulaNeumann() {
    std::array<matrixRow, 2> row{};
    row[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial('x', 0)
                               : calculateRepresentationFormulaNeumannOffAxial('x', 0);
    row[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial('y', 1)
                               : calculateRepresentationFormulaNeumannOffAxial('y', 1);
    for (int i = 0; i < 2; ++i) {
        if (!row[i].empty() && !iszero(normal[i])) {
            while (row[i].back().idx >= 4 * getNPts()) {
                partMatrixRow[i][row[i].back().idx - 4 * getNPts()] = row[i].back().value * normal[i];
                row[i].pop_back();
            }
            while (row[i].back().idx >= 2 * getNPts()) {
                rhsMatrixRow[i][row[i].back().idx - 2 * getNPts()] = row[i].back().value * normal[i];
                row[i].pop_back();
            }
            solMatrixRow[0] += row[i] * normal[i];
        }
    }
    if (getCondition() == 'N') approximatePhiAtBoundary(1);
    else if (getCondition() == 'n') approximatePhiAtAppend();
}

auto AGM::point::calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> AGM::matrixRow {
    auto Error = []() -> double {
        printError("AGM::point::calculateRepresentationFormulaNeumannOnAxial", "nullptr");
        return ZEROVALUE;
    };
    point *ptc = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(1) :
                 getAxialLine(axis)->back()->getIdx() == getIdx() ? *std::prev(getAxialLine(axis)->end() - 1) : nullptr;
    point *ptl =
            getAxialLine(axis)->front()->getIdx() == getIdx() ? this : getAxialLine(axis)->back()->getIdx() == getIdx()
                                                                       ? *std::prev(getAxialLine(axis)->end() - 2)
                                                                       : nullptr;
    point *ptr = getAxialLine(axis)->front()->getIdx() == getIdx() ? getAxialLine(axis)->at(2) :
                 getAxialLine(axis)->back()->getIdx() == getIdx() ? this : nullptr;
    std::string string =
            getAxialLine(axis)->front()->getIdx() == getIdx() ? "ND" : getAxialLine(axis)->back()->getIdx() == getIdx()
                                                                       ? "DN" : "";
    double tm = ptl ? ptl->getXy()[axisInt] : Error();
    double tb = ptc ? ptc->getXy()[axisInt] : Error();
    double tp = ptr ? ptr->getXy()[axisInt] : Error();
    double signPhi0 = axis == 'x' ? UNITVALUE : -UNITVALUE;
    auto gFunc{Greenfunction(tm, tb, tp, mp, mp)};

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
        row[ptr->getIdx()] = UNITVALUE;

        row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('l');
        row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('c');
        row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_ND('r');

        row[ptl->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('l');
        row[ptc->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('c');
        row[ptr->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_ND('r');

        row[ptl->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_ND('l') + gFunc.green_function_ND(tm);
        row[ptc->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_ND('c');
        row[ptr->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_ND('r');
    } else if (string == "DN") {
        row[ptl->getIdx()] = UNITVALUE;
        row[ptc->getIdx()] = -UNITVALUE;
        row[ptr->getIdx()] = mp * gFunc.green_function_DN(tp);

        row[ptl->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('l');
        row[ptc->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('c');
        row[ptr->getIdx() + getNPts()] = signPhi0 * gFunc.green_integral_DN('r');

        row[ptl->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('l');
        row[ptc->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('c');
        row[ptr->getIdx() + (axisInt + 2) * getNPts()] = gFunc.green_integral_DN('r');

        row[ptl->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_DN('l');
        row[ptc->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_DN('c');
        row[ptr->getIdx() + (axisInt + 4) * getNPts()] = gFunc.green_integral_t_DN('r') - gFunc.green_function_DN(tp);
    }
    auto c = -row[getIdx()];
    for (auto &item: row) {
        item.value /= c;
    }
    row.remove(getIdx());
    return row;
}

auto AGM::point::calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> AGM::matrixRow {
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
    for (const auto &item: {'y', 'x'}) {
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

        mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R');

        mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L');
        mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R');

        mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L');
        mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R');

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
    auto assignMatrix = [&row, &approximateSol, &linearApproximation, &sign](point *pt, point *ptr, point *ptl,
                                                                             double mp0, Greenfunction *func, double d,
                                                                             int i, int i0, char c) -> void {
        if (pt) {
            row[pt->getIdx()] += mp0 * func->green_function_ttau(d);
            row[pt->getIdx() + getNPts()] += sign * func->green_integral_tau(c);
            row[pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c);
            row[pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c);
        } else {
            row += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
            row += linearApproximation(ptr, ptl, sign * func->green_integral_tau(c), i, 1);
            row += linearApproximation(ptr, ptl, func->green_integral_tau(c), i, i0 + 2);
            row += linearApproximation(ptr, ptl, func->green_integral_ttau(c), i, i0 + 4);
        }
    };
    row[getIdx() + getNPts()] += sign * gFunc.green_integral_tau('c');
    row[getIdx() + (axisInt + 2) * getNPts()] += gFunc.green_integral_tau('c');
    row[getIdx() + (axisInt + 4) * getNPts()] += gFunc.green_integral_ttau('c') + UNITVALUE / mp;

    if (axis == 'x') {
        assignMatrix(element[E], element[EN], element[ES], -mp, &gFunc, tp, 1, 0, 'r');
        assignMatrix(element[W], element[WN], element[WS], mp, &gFunc, tm, 1, 0, 'l');
    } else if (axis == 'y') {
        assignMatrix(element[N], element[NE], element[NW], -mp, &gFunc, tp, 0, 1, 'r');
        assignMatrix(element[S], element[SE], element[SW], mp, &gFunc, tm, 0, 1, 'l');
    }
    return row;
}

void AGM::point::calculateRepresentationFormulaInterface() {
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
        double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ?
                                        ptr->getMp() : ptl->getCondition() == 'C' ?
                                                       ptl->getMp() : Error();
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
    auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i,
                                                  double d) -> matrixRow {
        double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
        auto func{Greenfunction(m, b, p, d, d)};
        auto mRow{matrixRow()};
        double sign = i ? -UNITVALUE : UNITVALUE;
        mRow[ptl->getIdx()] = d * func.green_function_t(m);
        mRow[ptr->getIdx()] = -d * func.green_function_t(p);

        mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptr, ptl);

        mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L');
        mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R');

        mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L');
        mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R');

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
    auto assignMatrix = [&row, &approximateSol, &linearApproximation, &isInterface](point *pt, point *ptr, point *ptl,
                                                                                    double mp0, Greenfunction *func,
                                                                                    double d, int i, int i0, char c,
                                                                                    char C) -> void {
        double sign = i0 ? -UNITVALUE : UNITVALUE;
        if (pt) {
            row[i0][pt->getIdx()] += mp0 * func->green_function_t(d);
            row[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral(C)
                                                             : sign * func->green_integral(c);
            row[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral(c);
            row[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_t(c);
        } else {
            row[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_t(d), i, std::abs(mp0));
            row[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral(C), i, 1)
                                   : linearApproximation(ptr, ptl, sign * func->green_integral(c), i, 1);
            row[i0] += linearApproximation(ptr, ptl, func->green_integral(c), i, i0 + 2);
            row[i0] += linearApproximation(ptr, ptl, func->green_integral_t(c), i, i0 + 4);
        }
    };
    row[0][getIdx()] = -UNITVALUE;
    if (!isInterface) row[0][getIdx() + getNPts()] = gFuncX.green_integral('c');
    row[0][getIdx() + 2 * getNPts()] = gFuncX.green_integral('c');
    row[0][getIdx() + 4 * getNPts()] = gFuncX.green_integral_t('c');
    assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
    assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

    row[1][getIdx()] = -UNITVALUE;
    if (!isInterface) row[1][getIdx() + getNPts()] = -gFuncY.green_integral('c');
    row[1][getIdx() + 3 * getNPts()] = gFuncY.green_integral('c');
    row[1][getIdx() + 5 * getNPts()] = gFuncY.green_integral_t('c');
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

void AGM::point::approximatePhiAtBoundary(int order) {
    auto findInnerPointOfBoundary = [this]() -> point * {
        for (const auto &item: {'x', 'y'}) {
            if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
                return getAxialLine(item)->at(2);
            }
            if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
                return *std::prev(getAxialLine(item)->end() - 2);
            }
        }
        std::cout << "(x, y) = (" << xy[0] << ", " << xy[1] << ")\n";
        std::cout << "boundary condition = " << getCondition() << "\n";
        printError("AGM::point::approximatePhiAtBoundary", "findInnerPointOfBoundary");
        return nullptr;
    };
    auto findStencil = [this]() -> EWNS {
        if (getAxialLine('x') && getAxialLine('x')->front()->getIdx() == getIdx()) return W;
        if (getAxialLine('x') && getAxialLine('x')->back()->getIdx() == getIdx()) return E;
        if (getAxialLine('y') && getAxialLine('y')->front()->getIdx() == getIdx()) return S;
        if (getAxialLine('y') && getAxialLine('y')->back()->getIdx() == getIdx()) return N;
        printError("AGM::point::approximatePhiAtBoundary", "findStencil");
        return E;
    };
    auto pt = *findInnerPointOfBoundary();
    auto ewns = findStencil();
    double tm{}, tb{}, tp{};
    if (getAxialLine('x')) {
        tm = pt[W]->getXy()[0];
        tb = pt.getXy()[0];
        tp = pt[E]->getXy()[0];
    } else if (getAxialLine('y')) {
        tm = pt[S]->getXy()[1];
        tb = pt.getXy()[1];
        tp = pt[N]->getXy()[1];
    }
    auto secondOrderExtrapolation = [this, &pt](int i, double m, double b, double p, EWNS ewns1, EWNS ewns2) -> void {
        double d{(p - m) * (p - b) * (b - m)}, t0{getXy()[i]};
        auto firstTerm = [&t0, &d](double w0, double w1) -> double { return t0 * t0 * (w0 - w1) / d; };
        auto secondTerm = [&t0, &d](double w0, double w1) -> double { return t0 * (w0 * w0 - w1 * w1) / d; };
        auto thirdTerm = [&d](double w0, double w1) -> double { return w0 * w1 * (w0 - w1) / d; };
        auto secondOrder = [&firstTerm, &secondTerm, &thirdTerm](double w0, double w1) -> double {
            return firstTerm(w0, w1) - secondTerm(w0, w1) + thirdTerm(w0, w1);
        };
        solMatrixRow[1][pt[ewns1]->getIdx() + getNPts()] = secondOrder(p, b);
        solMatrixRow[1][getIdx() + getNPts()] = -secondOrder(p, m);
        solMatrixRow[1][pt[ewns2]->getIdx() + getNPts()] = secondOrder(b, m);
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;
    };
    auto firstOrderExtrapolation = [this, &pt](int i, double t1, double t2, EWNS ewns0) -> void {
        double d{t2 - t1}, t0{getXy()[i]};
        auto firstTerm = [&t0, &d]() -> double { return t0 / d; };
        auto secondTerm = [&d](double t) -> double { return t / d; };
        auto firstOrder = [&firstTerm, &secondTerm](double t) -> double { return -firstTerm() + secondTerm(t); };
        solMatrixRow[1][pt[ewns0]->getIdx() + getNPts()] = -firstOrder(t1);
        solMatrixRow[1][pt.getIdx() + getNPts()] = firstOrder(t2);
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;
    };
    auto zeroOrderExtrapolation = [this, &pt](EWNS ewns1) -> void {
        solMatrixRow[1][pt[ewns1]->getIdx() + getNPts()] = UNITVALUE;
        solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;
    };
    if (getAxialLine('x')) {
        if (order == 2) {
            secondOrderExtrapolation(0, tm, tb, tp, W, E);
        } else if (order == 1) {
            if (ewns == E) {
                firstOrderExtrapolation(0, tb, tp, ewns);
            } else if (ewns == W) {
                firstOrderExtrapolation(0, tb, tm, ewns);
            } else {
                printError("AGM::point::approximatePhiAtBoundary", "ewns (which is %d) is wrong", ewns);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::point::approximatePhiAtBoundary", "order (which is %d) is wrong", order);
        }
    } else if (getAxialLine('y')) {
        if (order == 2) {
            secondOrderExtrapolation(1, tm, tb, tp, S, N);
        } else if (order == 1) {
            if (ewns == N) {
                firstOrderExtrapolation(1, tb, tp, ewns);
            } else if (ewns == S) {
                firstOrderExtrapolation(1, tb, tm, ewns);
            } else {
                printError("AGM::point::approximatePhiAtBoundary", "ewns (which is %d) is wrong", ewns);
            }
        } else if (order == 0) {
            zeroOrderExtrapolation(ewns);
        } else {
            printError("AGM::point::approximatePhiAtBoundary", "order (which is %d) is wrong", order);
        }
    } else {
        printError("AGM::point::approximatePhiAtBoundary", "getAxialLine error");
    }
}

void AGM::point::approximatePhiAtAppend() {
    for (const auto &item: {E, W, N, S}) {
        if (getElement()[item] && getIdx() != getElement()[item]->getIdx()) {
            solMatrixRow[1][getElement()[item]->getIdx() + getNPts()] = UNITVALUE;
            solMatrixRow[1][getIdx() + getNPts()] = -UNITVALUE;

            std::cout << "this point = " << getIdx() << "\n";
            std::cout << "near point = " << getElement()[item]->getIdx() << "\n\n";

            return;
        }
    }
}

void AGM::point::updateRightHandSide(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    switch (condition) {
        case 'C':
            updateRightHandSideCross(f, g);
            break;
        case 'D':
        case 'd':
            updateRightHandSideDirichlet(f, g);
            break;
        case 'N':
        case 'n':
            updateRightHandSideNeumann(f, g);
            break;
        case 'I':
            updateRightHandSideInterface(f, g);
            break;
        default:
            printError("AGM::point::updateRightHandSide", "condition (which is %d) is wrong", condition);
    }
}

void AGM::point::updateRightHandSideCross(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &item: rhsMatrixRow[0]) {
        rb[0] -= item.value * f(item.idx);
        rb[1] -= item.value * f(item.idx);
    }
    for (const auto &item: rhsMatrixRow[1]) {
        rb[0] -= item.value * g(item.idx - getNPts());
        rb[1] += item.value * g(item.idx - getNPts());
    }
}

void
AGM::point::updateRightHandSideDirichlet(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    rb[0] = values["bdv"];
    rb[1] = ZEROVALUE;
}

void AGM::point::updateRightHandSideNeumann(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    rb[0] = values["bdv"];
    rb[1] = ZEROVALUE;
    for (const auto &item: rhsMatrixRow) {
        for (const auto &item0: item) {
            if (item0.idx < getNPts()) {
                rb[0] -= item0.value * f(item0.idx);
            } else if (item0.idx < 2 * getNPts()) {
                rb[0] -= item0.value * g(item0.idx - getNPts());
            }
        }
    }
}

void
AGM::point::updateRightHandSideInterface(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    bool isinterface{solMatrixRow[1].size() == 1};
    rb[0] = rb[1] = ZEROVALUE;
    for (const auto &item: rhsMatrixRow[0]) {
        if (item.idx < getNPts()) {
            rb[0] -= item.value * f(item.idx);
            if (!isinterface) rb[1] -= item.value * f(item.idx);
        } else if (item.idx < 2 * getNPts()) {
            rb[0] -= item.value * g(item.idx - getNPts());
            if (!isinterface) rb[1] -= item.value * g(item.idx - getNPts());
        } else {
            printError("AGM::point::updateRightHandSideInterface");
        }
    }
    for (const auto &item: rhsMatrixRow[1]) {
        if (item.idx < getNPts()) {
            rb[0] -= item.value * f(item.idx);
            if (!isinterface) rb[1] += item.value * f(item.idx);
        } else if (item.idx < 2 * getNPts()) {
            rb[0] -= item.value * g(item.idx - getNPts());
            if (!isinterface) rb[1] += item.value * g(item.idx - getNPts());
        } else {
            printError("AGM::point::updateRightHandSideInterface");
        }
    }
}

void AGM::point::updateRightHandSidePart(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    switch (condition) {
        case 'C':
            updateRightHandSideCrossPart(f, g);
            break;
        case 'D':
        case 'd':
            updateRightHandSideDirichletPart(f, g);
            break;
        case 'N':
        case 'n':
            updateRightHandSideNeumannPart(f, g);
            break;
        case 'I':
            updateRightHandSideInterfacePart(f, g);
            break;
        default:
            printError("AGM::point::updateRightHandSidePart", "condition (which is %d) is wrong", condition);
    }
}

void
AGM::point::updateRightHandSideCrossPart(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    for (const auto &item: partMatrixRow[0]) {
        rb[0] -= item.value * f(item.idx);
        rb[1] -= item.value * f(item.idx);
    }
    for (const auto &item: partMatrixRow[1]) {
        rb[0] -= item.value * g(item.idx - getNPts());
        rb[1] += item.value * g(item.idx - getNPts());
    }
}

void
AGM::point::updateRightHandSideDirichletPart(const std::function<double(int)> &f, const std::function<double(int)> &g) {

}

void
AGM::point::updateRightHandSideNeumannPart(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    for (const auto &item: partMatrixRow) {
        for (const auto &item0: item) {
            if (item0.idx < getNPts()) {
                rb[0] -= item0.value * f(item0.idx);
            } else if (item0.idx < 2 * getNPts()) {
                rb[0] -= item0.value * g(item0.idx - getNPts());
            } else {
                printError("AGM::point::updateRightHandSideNeumannPart", "idx (which is %d) is too large", item0.idx);
            }
        }
    }
}

void
AGM::point::updateRightHandSideInterfacePart(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    bool isinterface{solMatrixRow[1].size() == 1};
    for (const auto &item: partMatrixRow[0]) {
        if (item.idx < getNPts()) {
            rb[0] -= item.value * f(item.idx);
            if (!isinterface) rb[1] -= item.value * f(item.idx);
        } else if (item.idx < 2 * getNPts()) {
            rb[0] -= item.value * g(item.idx - getNPts());
            if (!isinterface) rb[1] -= item.value * g(item.idx - getNPts());
        } else {
            printError("AGM::point::updateRightHandSideInterfacePart");
        }
    }
    for (const auto &item: partMatrixRow[1]) {
        if (item.idx < getNPts()) {
            rb[0] -= item.value * f(item.idx);
            if (!isinterface) rb[1] += item.value * f(item.idx);
        } else if (item.idx < 2 * getNPts()) {
            rb[0] -= item.value * g(item.idx - getNPts());
            if (!isinterface) rb[1] += item.value * g(item.idx - getNPts());
        } else {
            printError("AGM::point::updateRightHandSideInterfacePart");
        }
    }
}

void AGM::point::makeDerivatives() {
    switch (condition) {
        case 'C':
            makeDerivativesCross();
            break;
        case 'D':
        case 'd':
        case 'N':
        case 'n':
            makeDerivativesBoundary();
            break;
        case 'I':
            makeDerivativesInterface();
            break;
        default:
            printError("AGM::point::makeDerivatives", "condition (which is %c) is wrong", condition);
    }
}

void AGM::point::makeDerivativesCross() {
    double xm = element[W]->getXy()[0];
    double xb = getXy()[0];
    double xp = element[E]->getXy()[0];
    double ym = element[S]->getXy()[1];
    double yb = getXy()[1];
    double yp = element[N]->getXy()[1];
    auto gfuncX{Greenfunction(xm, xb, xp, mp, mp)};
    auto gfuncY{Greenfunction(ym, yb, yp, mp, mp)};
    deriMatrixRow[0][element[W]->getIdx()] = mp * gfuncX.green_function_ttau(xm);
    deriMatrixRow[0][element[E]->getIdx()] = -mp * gfuncX.green_function_ttau(xp);

    deriMatrixRow[0][element[W]->getIdx() + getNPts()] = gfuncX.green_integral_tau('l');
    deriMatrixRow[0][getIdx() + getNPts()] = gfuncX.green_integral_tau('c');
    deriMatrixRow[0][element[E]->getIdx() + getNPts()] = gfuncX.green_integral_tau('r');

    deriMatrixRow[0][element[W]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('l');
    deriMatrixRow[0][getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('c');
    deriMatrixRow[0][element[E]->getIdx() + 2 * getNPts()] = gfuncX.green_integral_tau('r');

    deriMatrixRow[0][element[W]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('l');
    deriMatrixRow[0][getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('c') + UNITVALUE / mp;
    deriMatrixRow[0][element[E]->getIdx() + 4 * getNPts()] = gfuncX.green_integral_ttau('r');

    deriMatrixRow[1][element[S]->getIdx()] = mp * gfuncY.green_function_ttau(ym);
    deriMatrixRow[1][element[N]->getIdx()] = -mp * gfuncY.green_function_ttau(yp);

    deriMatrixRow[1][element[S]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('l');
    deriMatrixRow[1][getIdx() + getNPts()] = -gfuncY.green_integral_tau('c');
    deriMatrixRow[1][element[N]->getIdx() + getNPts()] = -gfuncY.green_integral_tau('r');

    deriMatrixRow[1][element[S]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('l');
    deriMatrixRow[1][getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('c');
    deriMatrixRow[1][element[N]->getIdx() + 3 * getNPts()] = gfuncY.green_integral_tau('r');

    deriMatrixRow[1][element[S]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('l');
    deriMatrixRow[1][getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('c') + UNITVALUE / mp;
    deriMatrixRow[1][element[N]->getIdx() + 5 * getNPts()] = gfuncY.green_integral_ttau('r');
}

void AGM::point::makeDerivativesBoundary() {
    deriMatrixRow[0] = getAxialLine('x') ? calculateRepresentationFormulaNeumannOnAxial('x', 0)
                                         : calculateRepresentationFormulaNeumannOffAxial('x', 0);
    deriMatrixRow[1] = getAxialLine('y') ? calculateRepresentationFormulaNeumannOnAxial('y', 1)
                                         : calculateRepresentationFormulaNeumannOffAxial('y', 1);
}

void AGM::point::makeDerivativesInterface() {
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
        double rtv = pt ? pt->getMp() : ptr->getCondition() == 'C' ?
                                        ptr->getMp() : ptl->getCondition() == 'C' ?
                                                       ptl->getMp() : Error();
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
    auto approximateSol = [this, &checkMatrixRow](point *ptr, point *ptl, double coefficient, int i,
                                                  double d) -> matrixRow {
        double m{ptl->getXy()[i]}, b{getXy()[i]}, p{ptr->getXy()[i]};
        auto func{Greenfunction(m, b, p, d, d)};
        auto mRow{matrixRow()};
        double sign = i ? -UNITVALUE : UNITVALUE;
        mRow[ptl->getIdx()] = d * func.green_function_t(m);
        mRow[ptr->getIdx()] = -d * func.green_function_t(p);

        mRow[ptl->getIdx() + getNPts()] = sign * func.green_integral('L');
        mRow[ptr->getIdx() + getNPts()] = sign * func.green_integral('R');

        checkMatrixRow(&mRow, ptr, ptl);

        mRow[ptl->getIdx() + (i + 2) * getNPts()] = func.green_integral('L');
        mRow[ptr->getIdx() + (i + 2) * getNPts()] = func.green_integral('R');

        mRow[ptl->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('L');
        mRow[ptr->getIdx() + (i + 4) * getNPts()] = func.green_integral_t('R');

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
                                                                                    double mp0, Greenfunction *func,
                                                                                    double d, int i, int i0, char c,
                                                                                    char C) -> void {
        double sign = i0 ? -UNITVALUE : UNITVALUE;
        if (pt) {
            deriMatrixRow[i0][pt->getIdx()] += mp0 * func->green_function_ttau(d);
            deriMatrixRow[i0][pt->getIdx() + getNPts()] += isInterface ? sign * func->green_integral_tau(C)
                                                                       : sign * func->green_integral_tau(c);
            deriMatrixRow[i0][pt->getIdx() + (i0 + 2) * getNPts()] += func->green_integral_tau(c);
            deriMatrixRow[i0][pt->getIdx() + (i0 + 4) * getNPts()] += func->green_integral_ttau(c);
        } else {
            deriMatrixRow[i0] += approximateSol(ptr, ptl, mp0 * func->green_function_ttau(d), i, std::abs(mp0));
            deriMatrixRow[i0] += isInterface ? linearApproximation(ptr, ptl, sign * func->green_integral_tau(C), i, 1)
                                             : linearApproximation(ptr, ptl, sign * func->green_integral_tau(c), i, 1);
            deriMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_tau(c), i, i0 + 2);
            deriMatrixRow[i0] += linearApproximation(ptr, ptl, func->green_integral_ttau(c), i, i0 + 4);
        }
    };
    if (!isInterface) deriMatrixRow[0][getIdx() + getNPts()] += gFuncX.green_integral_tau('c');
    deriMatrixRow[0][getIdx() + 2 * getNPts()] += gFuncX.green_integral_tau('c');
    deriMatrixRow[0][getIdx() + 4 * getNPts()] += gFuncX.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(getElement()[E], getElement()[EN], getElement()[ES], -mpe, &gFuncX, xp, 1, 0, 'r', 'R');
    assignMatrix(getElement()[W], getElement()[WN], getElement()[WS], mpw, &gFuncX, xm, 1, 0, 'l', 'L');

    if (!isInterface) deriMatrixRow[1][getIdx() + getNPts()] += -gFuncY.green_integral_tau('c');
    deriMatrixRow[1][getIdx() + 3 * getNPts()] += gFuncY.green_integral_tau('c');
    deriMatrixRow[1][getIdx() + 5 * getNPts()] += gFuncY.green_integral_ttau('c') + UNITVALUE / mp;
    assignMatrix(getElement()[N], getElement()[NE], getElement()[NW], -mpn, &gFuncY, yp, 0, 1, 'r', 'R');
    assignMatrix(getElement()[S], getElement()[SE], getElement()[SW], mps, &gFuncY, ym, 0, 1, 'l', 'L');
}

void AGM::point::calculateDerivatives(const std::vector<point> *points, const std::function<double(int)> &f,
                                      const std::function<double(int)> &g, const std::function<double(int)> &fp,
                                      const std::function<double(int)> &gp) {
    auto assignDerivatives = [&](int i) -> double {
        double d{};
        for (const auto &item: deriMatrixRow[i]) {
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
                printError("AGM::point::calculateDerivatives", "item.idx (which is %d) is too large", item.idx);
            }
        }
        return d;
    };
    values["dx"] = assignDerivatives(0);
    values["dy"] = assignDerivatives(1);
}

void AGM::point::approximateNaNDerivatives(std::vector<AGM::point> *points) {
    auto findInnerPointOfBoundary = [this]() -> point * {
        if (getCondition() == 'd' || getCondition() == 'n') {
            for (const auto &item: {E, W, N, S}) {
                if (getElement()[item] && getIdx() != getElement()[item]->getIdx()) {
                    return getElement()[item];
                }
            }
        }

        for (const auto &item: {'x', 'y'}) {
            if (getAxialLine(item) && getAxialLine(item)->front()->getIdx() == getIdx()) {
                return getAxialLine(item)->at(1);
            }
            if (getAxialLine(item) && getAxialLine(item)->back()->getIdx() == getIdx()) {
                return *std::prev(getAxialLine(item)->end() - 1);
            }
        }
        printInformation();
        printError("AGM::point::approximateNaNDerivatives", "findInnerPointOfBoundary");
        return nullptr;
    };
    if (std::isnan(values["dx"])) values["dx"] = findInnerPointOfBoundary()->getValue()["dx"];
    if (std::isnan(values["dy"])) values["dy"] = findInnerPointOfBoundary()->getValue()["dy"];
}

void AGM::point::calculateDerivativesTwice(const std::function<double(int)> &f, const std::function<double(int)> &g) {
    values["dxx"] = -(f(getIdx()) + values["phi"]) / mp;
    values["dyy"] = -(g(getIdx()) - values["phi"]) / mp;
}

void AGM::point::printInformation() {
    char c0[256]{""}, c1[256]{""}, c2[256]{""}, c3[256]{""};
    sprintf(c0, "(%f, %f), with boundary condition %c", getXy()[0], getXy()[1], getCondition());
    if (getCondition() == 'D' || getCondition() == 'N') {
        sprintf(c1, ", boundary value = %f", getValue()["bdv"]);
        if (getCondition() == 'N') {
            sprintf(c2, ", norval vector = (%f, %f)", getNormal()[0], getNormal()[1]);
        }
    }
    sprintf(c3, "\n");
    printf("%s%s%s%s", c0, c1, c2, c3);
}

AGM::point::~point() = default;
