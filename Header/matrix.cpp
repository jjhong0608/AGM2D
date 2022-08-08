//
// Created by 조준홍 on 2022/02/13.
//

#include "matrix.h"

template<typename pt>
AGM::matrix<pt>::matrix() = default;

template<typename pt>
AGM::matrix<pt>::matrix(std::vector<pt> *pts) : pts(pts) {}

template<typename pt>
auto AGM::matrix<pt>::getIa() const -> int * {
    return ia;
}

template<typename pt>
void AGM::matrix<pt>::setIa(int *pInt) {
    matrix::ia = pInt;
}

template<typename pt>
auto AGM::matrix<pt>::getJa() const -> int * {
    return ja;
}

template<typename pt>
void AGM::matrix<pt>::setJa(int *pInt) {
    matrix::ja = pInt;
}

template<typename pt>
auto AGM::matrix<pt>::getEnt() const -> double * {
    return ent;
}

template<typename pt>
void AGM::matrix<pt>::setEnt(double *pDouble) {
    matrix::ent = pDouble;
}

template<typename pt>
auto AGM::matrix<pt>::getPts() const -> std::vector<pt> * {
    return pts;
}

template<typename pt>
void AGM::matrix<pt>::setPts(std::vector<pt> *vector) {
    matrix::pts = vector;
}

template<typename pt>
void AGM::matrix<pt>::makeMatrix() {
    if (!ia) {
        int nEntry{};
        for (const auto &item: *pts) {
            for (const auto &item0: item.getSolMatrixRow()) {
                nEntry += int(item0.size());
            }
        }
        std::cout << "nEntry = " << nEntry << "\n";
        ia = new int[2 * pts->size() + 1];
        ja = new int[nEntry];
        ent = new double[nEntry];
        int iaIdx{}, jaIdx{}, iaValue{1};
        auto assignMatrix = [&](int i) -> void {
            for (const auto &item: *pts) {
                iaValue += int(item.getSolMatrixRow()[i].size());
                ia[iaIdx++] = iaValue;
                for (const auto &ele: item.getSolMatrixRow()[i]) {
                    ja[jaIdx] = ele.idx + 1;
                    ent[jaIdx++] = ele.value;
                }
            }
        };
        ia[iaIdx++] = iaValue;
        for (int i = 0; i < 2; ++i) assignMatrix(i);
    } else {
        printError("AGM::matrix<pt>::makeMatrix()");
    }
}

template<typename pt>
void AGM::matrix<pt>::factorizeMatrix() {
    int size = int(pts->size());
    pPram.n = size * 2;
    pPram.nrhs = 1;
    for (auto &i: pPram.iparm) i = 0;
    pPram.iparm[0] = 1;         /* No solver default */
    pPram.iparm[1] = 3;         /* Fill-in reordering from METIS */
    pPram.iparm[3] = 0;         /* No iterative-direct algorithm */
    pPram.iparm[4] = 0;         /* No user fill-in reducing permutation */
    pPram.iparm[5] = 0;         /* Write solution into x */
    pPram.iparm[6] = 0;         /* Not in use */
    pPram.iparm[7] = 2;         /* Max numbers of iterative refinement steps */
    pPram.iparm[8] = 0;         /* Not in use */
    pPram.iparm[9] = 13;        /* Perturb the pivot elements with 1E-13 */
    pPram.iparm[10] = 1;        /* Use nonsymmetric permutation and scaling MPS */
    pPram.iparm[11] = 0;        /* Conjugate transposed/transpose solve */
    pPram.iparm[12] = 1;        /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
    pPram.iparm[13] = 0;        /* Output: Number of perturbed pivots */
    pPram.iparm[14] = 0;        /* Not in use */
    pPram.iparm[15] = 0;        /* Not in use */
    pPram.iparm[16] = 0;        /* Not in use */
    pPram.iparm[17] = -1;       /* Output: Number of nonzeros in the factor LU */
    pPram.iparm[18] = -1;       /* Output: Mflops for LU factorization */
    pPram.iparm[19] = 0;        /* Output: Numbers of CG Iterations */
    pPram.maxfct = 1;           /* Maximum number of numerical factorizations. */
    pPram.mnum = 1;             /* Which factorization to use. */
    pPram.msglvl = 1;           /* Print statistical information in file */
    pPram.error = 0;            /* Initialize error flag */
    pPram.mtype = 11;
    pPram.iparm[60] = 1;

    for (auto &i: pPram.ppt) i = nullptr;

    pPram.phase = 12;
    pardiso(pPram.ppt, &pPram.maxfct, &pPram.mnum, &pPram.mtype, &pPram.phase, &pPram.n, ent, ia, ja, &pPram.idum,
            &pPram.nrhs, pPram.iparm, &pPram.msglvl, &pPram.ddum, &pPram.ddum, &pPram.error);
    if (pPram.error != 0) {
        printf("\nERROR during solution: %d", pPram.error);
        exit(3);
    }
    pPram.msglvl = 0;
//    pPram.iparm[12] = 2;
}

template<typename pt>
void AGM::matrix<pt>::calculateMatrix() {
    int size{int(pts->size())};
    auto *rb = new double[pPram.n];

    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        rb[i] = pts->at(i).getRb()[0];
        rb[i + size] = pts->at(i).getRb()[1];
    }
    double x[pPram.n];
    for (int i = 0; i < pPram.n; ++i) {
        x[i] = ZEROVALUE;
    }
    pPram.phase = 33;
    pardiso(pPram.ppt, &pPram.maxfct, &pPram.mnum, &pPram.mtype, &pPram.phase, &pPram.n, ent, ia, ja, &pPram.idum,
            &pPram.nrhs, pPram.iparm, &pPram.msglvl, rb, x, &pPram.error);
    if (pPram.error != 0) {
        printf("\nERROR during solution: %d", pPram.error);
        exit(3);
    }
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        pts->at(i)["sol"] = x[i];
        pts->at(i)["phi"] = x[i + size];
    }
//    std::cout << "pressure\n";
//    calculateResidual(x, rb);
    delete[] rb;
}

template<typename pt>
void AGM::matrix<pt>::calculateResidual(const double *x, const double *rb) {
    auto res = new double[2 * pts->size()];
    double maxRes{};
    int maxIdx{};
    #pragma omp parallel for
    for (int i = 0; i < 2 * pts->size(); ++i) {
        res[i] = ZEROVALUE;
        for (int j = ia[i]; j < ia[i + 1]; ++j) {
            res[i] += ent[j - 1] * x[ja[j - 1] - 1];
        }
        res[i] = std::fabs(res[i] - rb[i]);
    }
    for (int i = 0; i < 2 * pts->size(); ++i) {
//        std::cout << "res[" << i << "] = " << res[i] << "\n";
        if (maxRes < res[i]) {
            maxRes = res[i];
            maxIdx = i;
        }
    }

    std::ofstream f("/home/jjhong0608/docker/AGM2D/Matrix/Residual.dat");
    if (!f.is_open()) {
        std::cout << "Residual file is not opened\n";
    }
    for (int i = 0; i < 2 * pts->size(); ++i) {
        f << i << "\t" << pts->at(i % pts->size())[0] << "\t" << pts->at(i % pts->size())[1] << "\t" << res[i]
          << "\n";
    }
    f.close();

    std::cout << "------ Residual of the matrix |Ax - b| -----\n";
    std::cout << "Maximum Residual = " << maxRes << "\n";
    std::cout << maxIdx << "-th point\n";
    std::cout << "(x, y) = (" << pts->at(maxIdx % pts->size())[0] << ", " << pts->at(maxIdx % pts->size())[1] << ")\n";
    std::cout << "boundary condition = " << pts->at(maxIdx % pts->size()).getCondition() << "\n";
    std::cout << "--------------------------------------------\n";

    if (maxRes > 1e2) {
        exit(1);
    }
}

template<typename pt>
void AGM::matrix<pt>::releaseMatrix() {
    pPram.phase = -1;
    pardiso(pPram.ppt, &pPram.maxfct, &pPram.mnum, &pPram.mtype, &pPram.phase, &pPram.n, ent, ia, ja, &pPram.idum,
            &pPram.nrhs, pPram.iparm, &pPram.msglvl, &pPram.ddum, &pPram.ddum, &pPram.error);
}

template<typename pt>
AGM::matrix<pt>::~matrix() {
    delete[] ia;
    delete[] ja;
    delete[] ent;
}

template
class AGM::matrix<AGM::point>;

template
class AGM::matrix<AGM::pointHeat>;

template
class AGM::matrix<AGM::pointAxisymmetric>;