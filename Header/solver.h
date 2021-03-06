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

        std::vector<point> *getPts() const;

        void setPts(std::vector<point> *vector);

        void ellipticSolver();

        void axisymmetricEllipticSolver();

        void heatSolver();

        void NavierStokesSolver();

    private:
        std::vector<point> *pts{};
    };

}


#endif //AGM_SOLVER_H
