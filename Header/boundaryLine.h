//
// Created by NIMS-JUNHONG on 2022/07/04.
//

#ifndef AGM2D_BOUNDARYLINE_H
#define AGM2D_BOUNDARYLINE_H

#include "axialLine.h"

namespace AGM {
    class denseMatrix : public std::vector<std::vector<double>> {
    private:
        int row{}, col{};
    public:
        denseMatrix(int row, int col);

        virtual ~denseMatrix();

        denseMatrix operator+(const denseMatrix &src);

        denseMatrix operator-(const denseMatrix &src);

        denseMatrix operator*(double d);

        denseMatrix operator*(const denseMatrix &src);

        AGM::vector operator*(const AGM::vector &src);

        denseMatrix &operator=(const denseMatrix &src);
    };

    class line2D {
    private:
        vector start{}, end{}, normal{}, tangent{};
        double length{};

    public:
        line2D();

        line2D(const vector& start, const vector &anEnd);

        line2D(double s0, double s1, double e0, double e1);

        void calcProperties();

        [[nodiscard]] vector &getStart();

        [[nodiscard]] const vector &getStart() const;

        void setStart(const vector &vector);

        [[nodiscard]] vector &getAnEnd();

        [[nodiscard]] const vector &getAnEnd() const;

        void setAnEnd(const vector &vector);

        [[nodiscard]] vector &getNormal();

        [[nodiscard]] const vector &getNormal() const;

        [[nodiscard]] vector &getTangent();

        [[nodiscard]] const vector &getTangent() const;

        [[nodiscard]] double getLength() const;

        bool iscross(const line2D &src, AGM::vector &vec);

        virtual ~line2D();
    };

    class boundaryLine2D : public line2D {
    private:
        char condition{};
        double boundary_value{};
    public:
        boundaryLine2D();

        boundaryLine2D(const vector &start, const vector &anEnd, char condition, double boundaryValue);

        boundaryLine2D(double s0, double s1, double e0, double e1, char condition, double boundaryValue);

        [[nodiscard]] char getCondition() const;

        void setCondition(char i);

        [[nodiscard]] double getBoundaryValue() const;

        void setBoundaryValue(double boundaryValue);
    };
}


#endif //AGM2D_BOUNDARYLINE_H
