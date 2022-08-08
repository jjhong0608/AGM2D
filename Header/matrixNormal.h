//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#ifndef AGM_MATRIXNORMAL_H
#define AGM_MATRIXNORMAL_H

#include "matrixMulti.h"
#include "StdVector"
#include "Eigen/Sparse"

namespace AGM {
    template<typename pt>
    class matrixNormal : public matrix<pt> {
    public:
        matrixNormal();

        explicit matrixNormal(std::vector<pt> *pts);

        matrixNormal(std::vector<pt> *pts, int fixedPointIdx);

        virtual ~matrixNormal();

        auto getFixedPointIdx() const -> int;

        void setFixedPointIdx(int i);

        void makeMatrix() override;

        void calculateMatrix() override;

        void factorizeMatrix() override;


    private:
        Eigen::SparseMatrix<double, Eigen::RowMajor> A{};
        int fixedPointIdx{};
        int *iaT{}, *jaT{};
        double *entT{};
    };

}


#endif //AGM_MATRIXNORMAL_H
