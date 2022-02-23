//
// Created by 조준홍 on 2022/02/13.
//

#ifndef AGM_MATRIX_H
#define AGM_MATRIX_H

#include "function.h"
#include "mkl_pardiso.h"

namespace AGM {
    struct pardisoParameters {
        void *ppt[64];
        int maxfct, mnum, mtype, phase, n, idum, nrhs, msglvl, error;
        int iparm[64];
        double ddum;
    };

    template<typename pt>
    class matrix {
    protected:
        int *ia{}, *ja{};
        double *ent{};
        std::vector<pt> *pts{};
        pardisoParameters pPram{};

    public:
        matrix();

        explicit matrix(std::vector<pt> *pts);

        virtual ~matrix();

        int *getIa() const;

        void setIa(int *pInt);

        int *getJa() const;

        void setJa(int *pInt);

        double *getEnt() const;

        void setEnt(double *pDouble);

        std::vector<pt> *getPts() const;

        void setPts(std::vector<pt> *vector);

        virtual void makeMatrix();

        virtual void factorizeMatrix();

        virtual void calculateMatrix();

        void releaseMatrix();
    };



}


#endif //AGM_MATRIX_H
