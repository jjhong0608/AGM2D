//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_MATRIXROW_H
#define AGM_MATRIXROW_H

#include "axialLine.h"

namespace AGM {
    struct matrixElement {
        int idx{};
        double value{};
    };

    class matrixRow : public std::vector<matrixElement> {
    public:
        void remove(int i);

        double &operator[](int i);

        matrixRow operator+(const matrixRow &src) const;

        matrixRow operator-(const matrixRow &src) const;

        matrixRow operator*(double d) const;

        matrixRow operator+=(const matrixRow &src);

        matrixRow operator-=(const matrixRow &src);
    };
}


#endif //AGM_MATRIXROW_H
