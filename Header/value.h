//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_VALUE_H
#define AGM_VALUE_H

#include "coordinate.h"

namespace AGM {
    class value : public std::array<double, 8> {

    public:
        value();

        virtual ~value();

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;
    };
}

#endif //AGM_VALUE_H
