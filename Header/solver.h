//
// Created by 조준홍 on 2022/02/19.
//

#ifndef AGM_SOLVER_H
#define AGM_SOLVER_H

#include "writeFileMultiple.h"

namespace AGM {
    class solver {
    public:
        explicit solver(std::vector<point> *pts);

        virtual ~solver();

        auto getPts() const -> std::vector<point> *;

        void setPts(std::vector<point> *vector);

        void ellipticSolver();

        void streamSolver();

        void axisymmetricEllipticSolver();

        void heatSolver();

        void NavierStokesSolver();

        void FluidStructureInteraction();

    private:
        std::vector<point> *pts{};
    };

}


#endif //AGM_SOLVER_H
