//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#include "matrixMulti.h"

template<typename pt>
AGM::matrixMulti<pt>::matrixMulti() = default;

template<typename pt>
AGM::matrixMulti<pt>::matrixMulti(std::vector<pt> *pts, std::vector<pt> *pts0) : matrix<pt>(pts), pts0(pts0) {}

template<typename pt>
void AGM::matrixMulti<pt>::factorizeMatrix() {
    int size = int(matrix<pt>::pts->size());
    matrix<pt>::pPram.n = size * 2;
    matrix<pt>::pPram.nrhs = 2;
    for (auto &i: matrix<pt>::pPram.iparm) i = 0;
    matrix<pt>::pPram.iparm[0] = 1;         /* No solver default */
    matrix<pt>::pPram.iparm[1] = 3;         /* Fill-in reordering from METIS */
    matrix<pt>::pPram.iparm[3] = 0;         /* No iterative-direct algorithm */
    matrix<pt>::pPram.iparm[4] = 0;         /* No user fill-in reducing permutation */
    matrix<pt>::pPram.iparm[5] = 0;         /* Write solution into x */
    matrix<pt>::pPram.iparm[6] = 0;         /* Not in use */
    matrix<pt>::pPram.iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    matrix<pt>::pPram.iparm[8] = 0;         /* Not in use */
    matrix<pt>::pPram.iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    matrix<pt>::pPram.iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    matrix<pt>::pPram.iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    matrix<pt>::pPram.iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    matrix<pt>::pPram.iparm[13] = 0;        /* Output: Number of perturbed pivots */
    matrix<pt>::pPram.iparm[14] = 0;        /* Not in use */
    matrix<pt>::pPram.iparm[15] = 0;        /* Not in use */
    matrix<pt>::pPram.iparm[16] = 0;        /* Not in use */
    matrix<pt>::pPram.iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    matrix<pt>::pPram.iparm[18] = -1;       /* Output: Mflops for LU factorization */
    matrix<pt>::pPram.iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    matrix<pt>::pPram.maxfct = 1;           /* Maximum number of numerical factorizations. */
    matrix<pt>::pPram.mnum = 1;             /* Which factorization to use. */
    matrix<pt>::pPram.msglvl = 1;           /* Print statistical information in file */
    matrix<pt>::pPram.error = 0;            /* Initialize error flag */
    matrix<pt>::pPram.mtype = 11;
    matrix<pt>::pPram.iparm[60] = 1;

    for (auto &i: matrix<pt>::pPram.ppt) i = nullptr;

    matrix<pt>::pPram.phase = 12;
    pardiso(matrix<pt>::pPram.ppt, &matrix<pt>::pPram.maxfct, &matrix<pt>::pPram.mnum, &matrix<pt>::pPram.mtype,
            &matrix<pt>::pPram.phase, &matrix<pt>::pPram.n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja,
            &matrix<pt>::pPram.idum, &matrix<pt>::pPram.nrhs, matrix<pt>::pPram.iparm, &matrix<pt>::pPram.msglvl,
            &matrix<pt>::pPram.ddum, &matrix<pt>::pPram.ddum, &matrix<pt>::pPram.error);
    if (matrix<pt>::pPram.error != 0) {
        printf("\nERROR during solution: %d", matrix<pt>::pPram.error);
        exit(3);
    }
    matrix<pt>::pPram.msglvl = 0;
}

template<typename pt>
void AGM::matrixMulti<pt>::calculateMatrix() {
    int size{int(matrix<pt>::pts->size())};
    auto *rb = new double[2 * matrix<pt>::pPram.n];
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        rb[i] = matrix<pt>::pts->at(i).getRb()[0];
        rb[i + size] = matrix<pt>::pts->at(i).getRb()[1];
        rb[i + 2 * size] = pts0->at(i).getRb()[0];
        rb[i + 3 * size] = pts0->at(i).getRb()[1];
    }
    double x[2 * matrix<pt>::pPram.n];
    for (auto &i: x) i = ZEROVALUE;
    matrix<pt>::pPram.phase = 33;
    pardiso(matrix<pt>::pPram.ppt, &matrix<pt>::pPram.maxfct, &matrix<pt>::pPram.mnum, &matrix<pt>::pPram.mtype,
            &matrix<pt>::pPram.phase, &matrix<pt>::pPram.n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja,
            &matrix<pt>::pPram.idum, &matrix<pt>::pPram.nrhs, matrix<pt>::pPram.iparm, &matrix<pt>::pPram.msglvl, rb, x,
            &matrix<pt>::pPram.error);
    if (matrix<pt>::pPram.error != 0) {
        printf("\nERROR during solution: %d", matrix<pt>::pPram.error);
        exit(3);
    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        matrix<pt>::pts->at(i)["sol"] = x[i];
        matrix<pt>::pts->at(i)["phi"] = x[i + size];
        pts0->at(i)["sol"] = x[i + 2 * size];
        pts0->at(i)["phi"] = x[i + 3 * size];
    }
    delete[] rb;
}

template<typename pt>
AGM::matrixMulti<pt>::~matrixMulti() = default;

template
class AGM::matrixMulti<AGM::point>;

template
class AGM::matrixMulti<AGM::pointHeat>;


