//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_POINT_H
#define AGM_POINT_H

#include <vector>
#include "value.h"

namespace AGM {
    class point;

    using axialElement = std::array<point *, 12>;

    class point {
    protected:
        int idx{};
        coordinate xy{}, normal{};
        double mp{};
        char condition{};
        axialElement element{};
        value values{};
        std::array<matrixRow, 2> solMatrixRow{}, deriMatrixRow{}, rhsMatrixRow{}, partMatrixRow{};
        std::array<double, 2> rb{}, dv{};
        std::array<axialLine *, 2> aline{};
        static int nPts;
        static std::vector<axialLine> *xline, *yline;
        static std::vector<boundaryLine2D> *bdLine;
        static std::vector<point> *pts;

    public:
        point();

        explicit point(int idx);

        explicit point(const coordinate &xy);

        point(int idx, const coordinate &xy);

        point(const coordinate &xy, double mp);

        point(int idx, const coordinate &xy, double mp);

        virtual ~point();

        int getIdx() const;

        void setIdx(int i);

        const coordinate &getXy() const;

        void setXy(const coordinate &coordinate);

        const coordinate &getNormal() const;

        void setNormal(const coordinate &coordinate);

        double getMp() const;

        void setMp(double d);

        char getCondition() const;

        void setCondition(char i);

        const std::array<point *, 12> &getElement() const;

        void setElement(const std::array<point *, 12> &array);

        const value &getValue() const;

        void setValue(const value &value);

        const std::array<matrixRow, 2> &getSolMatrixRow() const;

        void setSolMatrixRow(const std::array<matrixRow, 2> &row);

        const std::array<matrixRow, 2> &getDeriMatrixRow() const;

        void setDeriMatrixRow(const std::array<matrixRow, 2> &row);

        const std::array<double, 2> &getRb() const;

        void setRb(const std::array<double, 2> &array);

        const std::array<AGM::axialLine *, 2> &getAxialLine() const;

        axialLine *&getAxialLine(char i);

        void setAxialLine(const std::array<AGM::axialLine *, 2> &array);

        void setAxialLine(AGM::axialLine *line, char i);

        static int getNPts();

        static void setNPts(int i);

        static std::vector<axialLine> *&getAxialLines(char i);

        static void setAxialLines(std::vector<axialLine> *line, char i);

        static std::vector<boundaryLine2D> *getBdLine();

        static void setBdLine(std::vector<boundaryLine2D> *line);

        static std::vector<point> *getPts();

        static void setPts(std::vector<point> *vector);

        double &operator[](int i);

        const double &operator[](int i) const;

        point *&operator[](EWNS ewns);

        double &operator[](const std::string &string);

        const double &operator[](const std::string &string) const;

        double operator-(const point &src);

        double operator-(axialLine &src);

        point &operator=(const point &src);

        void findStencil();

        void findStencilBoundary();

        void findStencilAppendBoundary();

        void findStencilInterface();

        void findStencil(std::vector<AGM::point> *src, std::vector<AGM::point> *tgt);

        void calculateRepresentationFormula();

        virtual void calculateRepresentationFormulaCross();

        void calculateRepresentationFormulaDirichlet();

        void calculateRepresentationFormulaNeumann();

        virtual matrixRow calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt);

        virtual matrixRow calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt);

        virtual void calculateRepresentationFormulaInterface();

        void approximatePhiAtBoundary(int order);

        void approximatePhiAtAppend();

        void updateRightHandSide(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideCross(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideDirichlet(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideNeumann(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideInterface(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSidePart(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideCrossPart(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideDirichletPart(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideNeumannPart(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void updateRightHandSideInterfacePart(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void makeDerivatives();

        virtual void makeDerivativesCross();

        void makeDerivativesBoundary();

        virtual void makeDerivativesInterface();

        virtual void calculateDerivatives(const std::vector<point> *points, const std::function<double(int)> &f,
                                          const std::function<double(int)> &g, const std::function<double(int)> &fp,
                                          const std::function<double(int)> &gp);

        virtual void approximateNaNDerivatives(std::vector<point> *points);

        void calculateDerivativesTwice(const std::function<double(int)> &f, const std::function<double(int)> &g);
    };

}


#endif //AGM_POINT_H
