///
// Created by NIMS-JUNHONG on 2022/07/04.
//

#include "boundaryLine.h"

AGM::denseMatrix::denseMatrix(int row, int col) : vector(row),
                                                  row(row),
                                                  col(col) {
    for (auto &i: *this) {
        i = std::vector<double>(col);
    }
}

AGM::denseMatrix::~denseMatrix() = default;

auto AGM::denseMatrix::operator+(const AGM::denseMatrix &src) -> AGM::denseMatrix {
    if (row != src.row || col != src.col)
        printError("AGM::denseMatrix::operator+",
                   "size of denseMatrix (which is (%d, %d)) is not equal size of source denseMatrix (which is (%d, %d))",
                   row,
                   col, src.row, src.col);
    denseMatrix denseMatrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            denseMatrix[i][j] = at(i).at(j) + src.at(i).at(j);
        }
    }
    return denseMatrix;
}

auto AGM::denseMatrix::operator-(const AGM::denseMatrix &src) -> AGM::denseMatrix {
    if (row != src.row || col != src.col)
        printError("AGM::denseMatrix::operator-",
                   "the size of denseMatrix (which is (%d, %d)) is not equal to the size of source denseMatrix (which is (%d, %d))",
                   row, col, src.row, src.col);
    denseMatrix mat = denseMatrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = at(i).at(j) - src.at(i).at(j);
        }
    }
    return mat;
}

auto AGM::denseMatrix::operator*(double d) -> AGM::denseMatrix {
    denseMatrix mat = denseMatrix(row, col);
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            mat[i][j] = at(i).at(j) * d;
        }
    }
    return mat;
}

auto AGM::denseMatrix::operator*(const AGM::denseMatrix &src) -> AGM::denseMatrix {
    if (col != src.row)
        printError("AGM::denseMatrix::operator*",
                   "the column size of denseMatrix (which is %d) is not equal to the row size of source denseMatrix (which is %d)",
                   col, src.row);
    denseMatrix mat = denseMatrix(row, src.col);
    double sum;
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < src.col; ++j) {
            sum = ZEROVALUE;
            for (int k = 0; k < col; ++k) {
                sum += at(i).at(k) * src.at(k).at(j);
            }
            mat[i][j] = sum;
        }
    }
    return mat;
}

auto AGM::denseMatrix::operator*(const AGM::vector &src) -> AGM::vector {
    if (col != src.size())
        printError("AGM::denseMatrix::operator*",
                   "the column size of denseMatrix (which is %d) is not equal to the size of source vector (which is %d)",
                   col, src.size());
    AGM::vector vec = AGM::vector{};
    double sum{};
    for (int i = 0; i < row; ++i) {
        sum = ZEROVALUE;
        for (int j = 0; j < col; ++j) {
            sum += at(i).at(j) * src.at(j);
        }
        vec.emplace_back(sum);
    }
    return vec;
}

auto AGM::denseMatrix::operator=(const AGM::denseMatrix &src) -> AGM::denseMatrix & {
    row = src.row;
    col = src.col;
    resize(row);
    for (auto &i: *this) {
        i.resize(col);
    }
    for (int i = 0; i < row; ++i) {
        for (int j = 0; j < col; ++j) {
            at(i).at(j) = src.at(i).at(j);
        }
    }
    return *this;
}

AGM::line2D::line2D() = default;

AGM::line2D::line2D(const AGM::vector &start, const AGM::vector &anEnd) : start(start), end(anEnd) {
}

AGM::line2D::line2D(double s0, double s1, double e0, double e1) {
    start.emplace_back(s0);
    start.emplace_back(s1);
    end.emplace_back(e0);
    end.emplace_back(e1);
}

void AGM::line2D::calcProperties() {
    auto line = end - start;
    length = line.norm();
    normal = vector(std::vector<double>{-line[1] / length, line[0] / length});
    tangent = line.unitVector();
}

auto AGM::line2D::getStart() -> AGM::vector & {
    return start;
}

auto AGM::line2D::getStart() const -> const AGM::vector & {
    return start;
}

void AGM::line2D::setStart(const AGM::vector &vector) {
    if (vector.size() != 2) {
        printError("AGM::line2D::setStart", "start vector size (which is %d) is not 2", vector.size());
    }
    start = vector;
}

auto AGM::line2D::getAnEnd() -> AGM::vector & {
    return end;
}

auto AGM::line2D::getAnEnd() const -> const AGM::vector & {
    return end;
}

void AGM::line2D::setAnEnd(const AGM::vector &vector) {
    if (vector.size() != 2) {
        printError("AGM::line2D::setAnEnd", "end vector size (which is %d) is not 2", vector.size());
    }
    end = vector;
}

auto AGM::line2D::getNormal() -> AGM::vector & {
    return normal;
}

auto AGM::line2D::getNormal() const -> const AGM::vector & {
    return normal;
}

auto AGM::line2D::getTangent() -> AGM::vector & {
    return tangent;
}

auto AGM::line2D::getTangent() const -> const AGM::vector & {
    return tangent;
}

auto AGM::line2D::getLength() const -> double {
    return length;
}

auto AGM::line2D::iscross(const AGM::line2D &src, AGM::vector &vec) -> bool {
    if (isclose(std::fabs((src.end - src.start).unitVector() * tangent), UNITVALUE)) {
        return false;
    }
    auto projection = [this]() -> denseMatrix {
        denseMatrix mat = denseMatrix(2, 2);
        mat[0][0] = UNITVALUE - tangent[0] * tangent[0];
        mat[0][1] = -tangent[0] * tangent[1];
        mat[0][1] = -tangent[0] * tangent[1];
        mat[1][1] = UNITVALUE - tangent[1] * tangent[1];
        return mat;
    };
    auto P{projection()};
    auto projectedStart{P * src.start};
    auto projectedEnd{P * src.end};
    auto projectedPoint{P * start};

    if ((projectedPoint - projectedStart) * (projectedPoint - projectedEnd) < NEARZERO) {
        double lenStart{(projectedPoint - projectedStart).norm()};
        double lenEnd{(projectedPoint - projectedEnd).norm()};
        double ratioStart{lenStart / (lenStart + lenEnd)};
        double t{((src.start + (src.end - src.start) * ratioStart - start) * (end - start)) /
                 std::pow((end - start).norm(), 2)};
        if (ispositive(t * (t - UNITVALUE)) || std::isnan(ratioStart)) {
            return false;
        } else {
            double len{src.getLength()};
            vec = src.start + ((src.end - src.start) * (ratioStart * len)) / len;
            return true;
        }
    } else {
        return false;
    }
}

AGM::line2D::~line2D() = default;

AGM::boundaryLine2D::boundaryLine2D() = default;

AGM::boundaryLine2D::boundaryLine2D(double s0, double s1, double e0, double e1, char condition, double boundaryValue)
        : line2D(s0, s1, e0, e1), condition(condition), boundary_value(boundaryValue) {}

AGM::boundaryLine2D::boundaryLine2D(const AGM::vector &start, const AGM::vector &anEnd, char condition,
                                    double boundaryValue) : line2D(start, anEnd), condition(condition),
                                                            boundary_value(boundaryValue) {}

auto AGM::boundaryLine2D::getCondition() const -> char {
    return condition;
}

void AGM::boundaryLine2D::setCondition(char i) {
    boundaryLine2D::condition = i;
}

auto AGM::boundaryLine2D::getBoundaryValue() const -> double {
    return boundary_value;
}

void AGM::boundaryLine2D::setBoundaryValue(double boundaryValue) {
    boundary_value = boundaryValue;
}