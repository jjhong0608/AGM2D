//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_COORDINATE_H
#define AGM_COORDINATE_H

#include "matrixRow.h"

namespace AGM {
    class coordinate : public std::array<double, 2> {
    public:
        coordinate();

        virtual ~coordinate();

        coordinate(double x, double y);

        [[nodiscard]] double norm() const;

        coordinate operator+(const coordinate &src) const;

        coordinate operator-(const coordinate &src) const;

        coordinate operator*(double d) const;

        bool operator==(const coordinate &src) const;

        bool operator!=(const coordinate &rhs) const;
    };
}

#endif //AGM_COORDINATE_H
