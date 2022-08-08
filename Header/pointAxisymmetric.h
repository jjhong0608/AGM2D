//
// Created by NIMS-JUNHONG on 2022/07/20.
//

#ifndef AGM2D_POINTAXISYMMETRIC_H
#define AGM2D_POINTAXISYMMETRIC_H

#include "point.h"

namespace AGM {
    class pointAxisymmetric : public point {
    private:
        bool is_on_axis{false};

    public:
        auto isOnAxis() const -> bool;

        void setIsOnAxis(bool isOnAxis);

        void findStencil(const axialElement *axialElement1, std::vector<pointAxisymmetric> *vector);

        void checkOnAxis();

        void EquationOnAxis();

        void calculateRepresentationFormulaCross() override;

        auto calculateRepresentationFormulaCrossSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaCrossSymmetricNearAxis() -> matrixRow;

        auto calculateRepresentationFormulaCrossNonSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) -> matrixRow override;

        auto calculateRepresentationFormulaNeumannOnAxialSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaNeumannOnAxialNonSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) -> matrixRow override;

        auto calculateRepresentationFormulaNeumannOffAxialSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaNeumannOffAxialNonSymmetric() -> matrixRow;

        void calculateRepresentationFormulaInterface() override;

        auto calculateRepresentationFormulaInterfaceSymmetric() -> matrixRow;

        auto calculateRepresentationFormulaInterfaceSymmetricNearAxis() -> matrixRow;

        auto calculateRepresentationFormulaInterfaceNonSymmetric() -> matrixRow;

        void makeDerivativesCross() override;

        void makeDerivativesCrossSymmetric();

        void makeDerivativesCrossSymmetricNearAxis();

        void makeDerivativesCrossNonSymmetric();

        void calculateDerivatives(const std::vector<pointAxisymmetric> *points, const std::function<double(int)> &f,
                                  const std::function<double(int)> &g, const std::function<double(int)> &fp,
                                  const std::function<double(int)> &gp);

        void approximateNaNDerivatives(std::vector<pointAxisymmetric> *points);

        void makeDerivativesInterface() override;

        void makeDerivativesInterfaceSymmetric();

        void makeDerivativesInterfaceSymmetricNearAxis();

        void makeDerivativesInterfaceNonSymmetric();
    };

}


#endif //AGM2D_POINTAXISYMMETRIC_H
