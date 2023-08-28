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
        std::array<matrixRow, 2> phiPressureMatrixRow{};
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

        [[nodiscard]] auto getIdx() const -> int;

        void setIdx(int i);

        [[nodiscard]] auto getXy() const -> const coordinate &;

        void setXy(const coordinate &coordinate);

        [[nodiscard]] auto getNormal() const -> const coordinate &;

        void setNormal(const coordinate &coordinate);

        [[nodiscard]] auto getMp() const -> double;

        void setMp(double d);

        [[nodiscard]] auto getCondition() const -> char;

        void setCondition(char i);

        [[nodiscard]] auto getElement() const -> const std::array<point *, 12> &;

        void setElement(const std::array<point *, 12> &array);

        [[nodiscard]] auto getValue() const -> const value &;

        void setValue(const value &value);

        [[nodiscard]] auto getSolMatrixRow() const -> const std::array<matrixRow, 2> &;

        void setSolMatrixRow(const std::array<matrixRow, 2> &row);

        [[nodiscard]] auto getDeriMatrixRow() const -> const std::array<matrixRow, 2> &;

        void setDeriMatrixRow(const std::array<matrixRow, 2> &row);

        [[nodiscard]] auto getPhiPressureMatrixRow() const -> const std::array<matrixRow, 2> &;

        void setPhiPressureMatrixRow(const std::array<matrixRow, 2> &row);

        [[nodiscard]] auto getRb() const -> const std::array<double, 2> &;

        void setRb(const std::array<double, 2> &array);

        [[nodiscard]] auto getAxialLine() const -> const std::array<AGM::axialLine *, 2> &;

        auto getAxialLine(char i) -> axialLine *&;

        void setAxialLine(const std::array<AGM::axialLine *, 2> &array);

        void setAxialLine(AGM::axialLine *line, char i);

        static auto getNPts() -> int;

        static void setNPts(int i);

        static auto getAxialLines(char i) -> std::vector<axialLine> *&;

        static void setAxialLines(std::vector<axialLine> *line, char i);

        static auto getBdLine() -> std::vector<boundaryLine2D> *;

        static void setBdLine(std::vector<boundaryLine2D> *line);

        static auto getPts() -> std::vector<point> *;

        static void setPts(std::vector<point> *vector);

        auto operator[](int i) -> double &;

        auto operator[](int i) const -> const double &;

        auto operator[](EWNS ewns) -> point *&;

        auto operator[](const std::string &string) -> double &;

        auto operator[](const std::string &string) const -> const double &;

        auto operator-(const point &src) -> double;

        auto operator-(axialLine &src) -> double;

        auto operator=(const point &src) -> point &;

        void findStencil();

        void findStencilBoundary();

        void findStencilAppendBoundary();

        void findStencilInterface();

        void findStencil(std::vector<AGM::point> *src, std::vector<AGM::point> *tgt);

        void calculateRepresentationFormula();

        virtual void calculateRepresentationFormulaCross();

        void calculateRepresentationFormulaDirichlet();

        void calculateRepresentationFormulaNeumann();

        virtual auto calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> matrixRow;

        virtual auto calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> matrixRow;

        virtual void calculateRepresentationFormulaInterface();

        void calculateRepresentationFormulaPhiPressure(char comp);

        void calculateRepresentationFormulaPhiPressureCross(char comp);

        void calculateRepresentationFormulaPhiPressureDirichlet();

        void calculateRepresentationFormulaPhiPressureInterface(char comp);

        void approximatePhiAtBoundary(int order);

        void approximatePhiAtBoundary1(int order);

        void approximatePhiAtBoundary2();

        void approximatePhiAtAppend();

        void approximateDiff(std::vector<point> *points);

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

        void updateRightHandSidePhiPressure(
                const std::function<double(int)> &f,
                const std::function<double(int)> &g,
                std::vector<point> *points
        );

        void updateRightHandSidePhiPressureCross(
                const std::function<double(int)> &f,
                const std::function<double(int)> &g,
                std::vector<point> *points
        );

        void updateRightHandSidePhiPressureDirichlet(
                const std::function<double(int)> &f,
                const std::function<double(int)> &g,
                std::vector<point> *points
        );

        void updateRightHandSidePhiPressureInterface(
                const std::function<double(int)> &f,
                const std::function<double(int)> &g,
                std::vector<point> *points
        );

        void makeDerivatives();

        virtual void makeDerivativesCross();

        void makeDerivativesBoundary();

        virtual void makeDerivativesInterface();

        void makePhiCoefficient(std::vector<point> *vector);

        virtual void calculateDerivatives(const std::vector<point> *points, const std::function<double(int)> &f,
                                          const std::function<double(int)> &g, const std::function<double(int)> &fp,
                                          const std::function<double(int)> &gp);

        virtual void approximateNaNDerivatives(std::vector<point> *points);

        void calculateDerivativesTwice(const std::function<double(int)> &f, const std::function<double(int)> &g);

        void printInformation();
    };

}


#endif //AGM_POINT_H
