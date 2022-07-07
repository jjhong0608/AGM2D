//
// Created by NIMS-JUNHONG on 2022/07/04.
//

#ifndef AGM2D_VECTOR_H
#define AGM2D_VECTOR_H

#include "util.h"

namespace AGM {
    class vector : public std::vector<double> {
    public:
        vector();

        explicit vector(const std::vector<double> &x);

        virtual ~vector();

        double norm();

        double dot(const vector &src);

        vector cross(const vector &src);

        vector unitVector();

        vector operator+(const vector &src);

        vector operator-(const vector &src);

        vector operator*(double d);

        vector operator/(double d);

        vector operator+(const vector &src) const;

        vector operator-(const vector &src) const;

        vector operator*(double d) const;

        vector operator/(double d) const;

        double operator*(const vector &src);

        vector &operator+=(const vector &src);

        vector &operator-=(const vector &src);

        vector &operator*=(double d);

        bool operator<(const vector &src);

        bool operator>(const vector &src);
    };

}

#endif //AGM2D_VECTOR_H
