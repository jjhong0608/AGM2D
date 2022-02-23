//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include "matrixRow.h"

void AGM::matrixRow::remove(int i) {
    for (int j = 0; j < size(); ++j) {
        if (at(j).idx == i) {
            erase(begin() + j);
        }
    }
}

double &AGM::matrixRow::operator[](int i) {
    if (empty()) {
        matrixElement args;
        args.idx = i;
        args.value = ZEROVALUE;
        emplace_back(args);
        return front().value;
    } else {
        for (int j = 0; j < size(); ++j) {
            if (at(j).idx == i) {
                return at(j).value;
            } else if (at(j).idx > i) {
                matrixElement args;
                args.idx = i;
                args.value = ZEROVALUE;
                emplace(begin() + j, args);
                return at(j).value;
            }
        }
    }
    matrixElement args;
    args.idx = i;
    args.value = ZEROVALUE;
    emplace_back(args);
    return back().value;
}

AGM::matrixRow AGM::matrixRow::operator+(const AGM::matrixRow &src) const {
    auto row = matrixRow();
    row = *this;
    for (const auto &i: src) {
        row[i.idx] += i.value;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator-(const AGM::matrixRow &src) const {
    auto row = matrixRow();
    row = *this;
    for (const auto &i: src) {
        row[i.idx] -= i.value;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator*(double d) const {
    auto row = matrixRow();
    row = *this;
    for (auto &i: row) {
        row[i.idx] *= d;
    }
    return row;
}

AGM::matrixRow AGM::matrixRow::operator+=(const AGM::matrixRow &src) {
    for (const auto &i: src) {
        (*this)[i.idx] += i.value;
    }
    return *this;
}

AGM::matrixRow AGM::matrixRow::operator-=(const AGM::matrixRow &src) {
    for (const auto &i: src) {
        (*this)[i.idx] -= i.value;
    }
    return *this;
}
