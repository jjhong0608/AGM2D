//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#include "axialLine.h"

AGM::axialLine::axialLine() = default;

AGM::axialLine::axialLine(char mark) : mark(mark) {}

char AGM::axialLine::getMark() const {
    return mark;
}

void AGM::axialLine::setMark(char i) {
    axialLine::mark = i;
}

double &AGM::axialLine::operator[](int i) {
    return coordinate[i];
}

double AGM::axialLine::operator-(AGM::axialLine &line) {
    if (mark == 'x') return coordinate[2] - line[2];
    else if (mark == 'y') return coordinate[0] - line[0];
    else printError("AGM::axialLine::operator-", "mark (which is %c) is wrong", mark);
    return ZEROVALUE;
}

AGM::axialLine::~axialLine() = default;