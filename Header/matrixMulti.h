//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#ifndef AGM_MATRIXMULTI_H
#define AGM_MATRIXMULTI_H

#include "matrix.h"

namespace AGM {
    template<typename pt>
    class matrixMulti : public matrix<pt> {
    public:
        matrixMulti();

        matrixMulti(std::vector<pt> *pts, std::vector<pt> *pts0);

        virtual ~matrixMulti();

        void calculateMatrix() override;

        void factorizeMatrix() override;


    private:
        std::vector<pt> *pts0{};
    };

}

#endif //AGM_MATRIXMULTI_H
