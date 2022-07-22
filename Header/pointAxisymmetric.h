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
        bool isOnAxis() const;

        void setIsOnAxis(bool isOnAxis);

        void findStencil(const axialElement *axialElement1, std::vector<pointAxisymmetric> *vector);

        void checkOnAxis();

        void EquationOnAxis();

        void calculateRepresentationFormulaCross() override;

        matrixRow calculateRepresentationFormulaCrossSymmetric();

        matrixRow calculateRepresentationFormulaCrossSymmetricNearAxis();

        matrixRow calculateRepresentationFormulaCrossNonSymmetric();

        matrixRow calculateRepresentationFormulaNeumannOnAxial(char axis, int axisInt) override;

        matrixRow calculateRepresentationFormulaNeumannOnAxialSymmetric();

        matrixRow calculateRepresentationFormulaNeumannOnAxialNonSymmetric();

        matrixRow calculateRepresentationFormulaNeumannOffAxial(char axis, int axisInt) override;

        matrixRow calculateRepresentationFormulaNeumannOffAxialSymmetric();

        matrixRow calculateRepresentationFormulaNeumannOffAxialNonSymmetric();

        void calculateRepresentationFormulaInterface() override;

        matrixRow calculateRepresentationFormulaInterfaceSymmetric();

        matrixRow calculateRepresentationFormulaInterfaceNonSymmetric();

        void makeDerivativesCross() override;

        void makeDerivativesCrossSymmetric();

        void makeDerivativesCrossNonSymmetric();

        void calculateDerivatives(const std::vector<pointAxisymmetric> *points, const std::function<double(int)> &f,
                                  const std::function<double(int)> &g, const std::function<double(int)> &fp,
                                  const std::function<double(int)> &gp);

        void approximateNaNDerivatives(std::vector<pointAxisymmetric> *points);

        void makeDerivativesInterface() override;

        void makeDerivativesInterfaceSymmetric();

        void makeDerivativesInterfaceNonSymmetric();
    };

}


#endif //AGM2D_POINTAXISYMMETRIC_H
