//
// Created by 조준홍 on 1/27/24.
//

#include "matrixStokesNormal.h"

template<typename pt>
AGM::matrixStokesNormal<pt>::matrixStokesNormal() = default;

template<typename pt>
AGM::matrixStokesNormal<pt>::matrixStokesNormal(
    std::vector<pt> *uvel,
    std::vector<pt> *vvel,
    std::vector<pt> *pres,
    int fixedPointIdx) : uvel(uvel), vvel(vvel), matrix<pt>(pres), fixedPointIdx(fixedPointIdx) {
}

template<typename pt>
auto AGM::matrixStokesNormal<pt>::getFixedPointIdx() const -> int {
  return fixedPointIdx;
}

template<typename pt>
void AGM::matrixStokesNormal<pt>::setFixedPointIdx(int fixed_point_idx) {
  matrixStokesNormal::fixedPointIdx = fixed_point_idx;
}

template<typename pt>
void AGM::matrixStokesNormal<pt>::makeMatrix() {
  auto size{int(matrix<pt>::pts->size())};
  auto a{Eigen::SparseMatrix<double, Eigen::RowMajor>(5 * size, 5 * size - 1)};
  for (const auto &item : *uvel) {
    for (const auto &row : item.getSolMatrixRow()[0]) {
      if (!iszero(row.value)) {
        a.insert(item.getIdx(), row.idx) = row.value;
      }
    }
    for (const auto &row : item.getPartMatrixRow()[0]) {
      if (!iszero(row.value)) {
        if (row.idx < fixedPointIdx) {
          a.insert(item.getIdx(), row.idx + 4 * size) = row.value;
        } else if (row.idx > fixedPointIdx) {
          a.insert(item.getIdx(), row.idx + 4 * size - 1) = row.value;
        }
      }
    }
  }
  for (const auto &item : *uvel) {
    for (const auto &row : item.getSolMatrixRow()[1]) {
      if (!iszero(row.value)) {
        a.insert(item.getIdx() + size, row.idx) = row.value;
      }
    }
    for (const auto &row : item.getPartMatrixRow()[0]) {
      if (!iszero(row.value)) {
        if (row.idx < fixedPointIdx) {
          a.insert(item.getIdx() + size, row.idx + 4 * size) = row.value;
        } else if (row.idx > fixedPointIdx) {
          a.insert(item.getIdx() + size, row.idx + 4 * size - 1) = row.value;
        }
      }
    }
  }
  for (const auto &item : *vvel) {
    for (const auto &row : item.getSolMatrixRow()[0]) {
      if (!iszero(row.value)) {
        a.insert(item.getIdx() + 2 * size, row.idx + 2 * size) = row.value;
      }
    }
    for (const auto &row : item.getPartMatrixRow()[1]) {
      if (!iszero(row.value)) {
        if (row.idx - size < fixedPointIdx) {
          a.insert(item.getIdx() + 2 * size, row.idx + 3 * size) = row.value;
        } else if (row.idx - size > fixedPointIdx) {
          a.insert(item.getIdx() + 2 * size, row.idx + 3 * size - 1) = row.value;
        }
      }
    }
  }
  for (const auto &item : *vvel) {
    for (const auto &row : item.getSolMatrixRow()[1]) {
      if (!iszero(row.value)) {
        a.insert(item.getIdx() + 3 * size, row.idx + 2 * size) = row.value;
      }
    }
    for (const auto &row : item.getPartMatrixRow()[1]) {
      if (!iszero(row.value)) {
        if (row.idx - size < fixedPointIdx) {
          a.insert(item.getIdx() + 3 * size, row.idx + 3 * size) = -row.value;
        } else if (row.idx - size > fixedPointIdx) {
          a.insert(item.getIdx() + 3 * size, row.idx + 3 * size - 1) = -row.value;
        }
      }
    }
  }
  auto presMatrix{std::vector<matrixRow>(size)};
  for (int i = 0; i < size; ++i) {
    presMatrix.at(i) = uvel->at(i).getPMatrixRow()[0] + vvel->at(i).getPMatrixRow()[1];
  }

  for (const auto &item : *matrix<pt>::pts) {
    for (const auto &row : presMatrix.at(item.getIdx())) {
      if (!iszero(row.value)) {
        if (row.idx - 4 * size < fixedPointIdx) {
          a.insert(item.getIdx() + 4 * size, row.idx) = row.value;
        } else if (row.idx - 4 * size > fixedPointIdx) {
          a.insert(item.getIdx() + 4 * size, row.idx - 1) = row.value;
        }
      }
    }
  }
  A = a.transpose() * a;
  A.makeCompressed();
  matrix<pt>::ia = new int[A.outerSize() + 1];
  matrix<pt>::ja = new int[A.nonZeros()];
  matrix<pt>::ent = new double[A.nonZeros()];
  auto oPtr{A.outerIndexPtr()};
  for (int j = 0; j < A.outerSize() + 1; ++j) {
    matrix<pt>::ia[j] = *oPtr + 1;
    ++oPtr;
  }
  auto iPtr{A.innerIndexPtr()};
  auto vPtr{A.valuePtr()};
  for (int j = 0; j < A.nonZeros(); ++j) {
    matrix<pt>::ja[j] = *iPtr + 1;
    matrix<pt>::ent[j] = *vPtr;
    ++iPtr;
    ++vPtr;
  }
  A = a.transpose();
  A.makeCompressed();
  iaT = new int[A.outerSize() + 1];
  jaT = new int[A.nonZeros()];
  entT = new double[A.nonZeros()];
  oPtr = A.outerIndexPtr();
  for (int j = 0; j < A.outerSize() + 1; ++j) {
    iaT[j] = *oPtr;
    ++oPtr;
  }
  iPtr = A.innerIndexPtr();
  vPtr = A.valuePtr();
  for (int j = 0; j < A.nonZeros(); ++j) {
    jaT[j] = *iPtr;
    entT[j] = *vPtr;
    ++iPtr;
    ++vPtr;
  }
}

template<typename pt>
void AGM::matrixStokesNormal<pt>::factorizeMatrix() {
  int size = int(matrix<pt>::pts->size());
  matrix<pt>::pPram.n = size * 5 - 1;
  matrix<pt>::pPram.nrhs = 1;
  for (auto &i : matrix<pt>::pPram.iparm) i = 0;
  matrix<pt>::pPram.iparm[0] = 1;   /* No solver default */
  matrix<pt>::pPram.iparm[1] = 3;   /* Fill-in reordering from METIS */
  matrix<pt>::pPram.iparm[3] = 0;   /* No iterative-direct algorithm */
  matrix<pt>::pPram.iparm[4] = 0;   /* No user fill-in reducing permutation */
  matrix<pt>::pPram.iparm[5] = 0;   /* Write solution into x */
  matrix<pt>::pPram.iparm[6] = 0;   /* Not in use */
  matrix<pt>::pPram.iparm[7] = 2;   /* Max numbers of iterative refinement steps */
  matrix<pt>::pPram.iparm[8] = 0;   /* Not in use */
  matrix<pt>::pPram.iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
  matrix<pt>::pPram.iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
  matrix<pt>::pPram.iparm[11] = 0;  /* Conjugate transposed/transpose solve */
  matrix<pt>::pPram.iparm[12] = 1;  /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
  matrix<pt>::pPram.iparm[13] = 0;  /* Output: Number of perturbed pivots */
  matrix<pt>::pPram.iparm[14] = 0;  /* Not in use */
  matrix<pt>::pPram.iparm[15] = 0;  /* Not in use */
  matrix<pt>::pPram.iparm[16] = 0;  /* Not in use */
  matrix<pt>::pPram.iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
  matrix<pt>::pPram.iparm[18] = -1; /* Output: Mflops for LU factorization */
  matrix<pt>::pPram.iparm[19] = 0;  /* Output: Numbers of CG Iterations */
  matrix<pt>::pPram.maxfct = 1;     /* Maximum number of numerical factorizations. */
  matrix<pt>::pPram.mnum = 1;       /* Which factorization to use. */
  matrix<pt>::pPram.msglvl = 1;     /* Print statistical information in file */
  matrix<pt>::pPram.error = 0;      /* Initialize error flag */
  matrix<pt>::pPram.mtype = 11;
  matrix<pt>::pPram.iparm[60] = 1;

  for (auto &i : matrix<pt>::pPram.ppt) i = nullptr;

  matrix<pt>::pPram.phase = 12;
  pardiso(matrix<pt>::pPram.ppt, &matrix<pt>::pPram.maxfct, &matrix<pt>::pPram.mnum, &matrix<pt>::pPram.mtype, &matrix<pt>::pPram.phase, &matrix<pt>::pPram.n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &matrix<pt>::pPram.idum, &matrix<pt>::pPram.nrhs, matrix<pt>::pPram.iparm, &matrix<pt>::pPram.msglvl, &matrix<pt>::pPram.ddum, &matrix<pt>::pPram.ddum, &matrix<pt>::pPram.error);
  if (matrix<pt>::pPram.error != 0) {
    printf("\nERROR during solution: %d", matrix<pt>::pPram.error);
    exit(3);
  }
}

template<typename pt>
void AGM::matrixStokesNormal<pt>::calculateMatrix() {
  int size{int(matrix<pt>::pts->size())};
  auto *rb = new double[matrix<pt>::pPram.n];
  auto *rb0 = new double[matrix<pt>::pPram.n + 1];
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    rb0[i] = uvel->at(i).getRb()[0];
    rb0[i + size] = uvel->at(i).getRb()[1];
    rb0[i + 2 * size] = vvel->at(i).getRb()[0];
    rb0[i + 3 * size] = vvel->at(i).getRb()[1];
    rb0[i + 4 * size] = uvel->at(i).getSrb()[0] + vvel->at(i).getSrb()[1];
  }
#pragma omp parallel for
  for (int j = 0; j < matrix<pt>::pPram.n; ++j) {
    rb[j] = ZEROVALUE;
    for (int k = iaT[j]; k < iaT[j + 1]; ++k) {
      rb[j] += entT[k] * rb0[jaT[k]];
    }
  }

  //    for (int i = 0; i < matrix<pt>::pPram.n + 1; ++i) {
  //        if (std::isnan(rb0[i])) {
  //            printf("Found NaN value at %d, ", i);
  //            printf(" / %d (%d), ", point::getNPts(), int(i / point::getNPts()));
  //            printf("condition = %c\n", uvel->at(i % point::getNPts()).getCondition());
  //        }
  //    }
  //
  //    for (int i = 0; i < matrix<pt>::pPram.n; ++i) {
  //        if (std::isnan(rb[i])) {
  //            printf("Found NaN1 value at %d, ", i);
  //            printf(" / %d (%d), ", point::getNPts(), int(i / point::getNPts()));
  //            printf("condition = %c\n", uvel->at(i % point::getNPts()).getCondition());
  //        }
  //    }
  //
  //    for (int i = 0; i < A.nonZeros(); ++i) {
  //        if (std::isnan(entT[i])) {
  //            printf("Found NaN2 value at %d, ", i);
  //            printf(" / %d (%d), ", point::getNPts(), int(i / point::getNPts()));
  //            printf("condition = %c\n", uvel->at(i % point::getNPts()).getCondition());
  //        }
  //    }
  //    printError("-----");

  double x[matrix<pt>::pPram.n];
  for (int j = 0; j < matrix<pt>::pPram.n; ++j) {
    x[j] = ZEROVALUE;
  }
  matrix<pt>::pPram.phase = 33;
  pardiso(matrix<pt>::pPram.ppt, &matrix<pt>::pPram.maxfct, &matrix<pt>::pPram.mnum, &matrix<pt>::pPram.mtype, &matrix<pt>::pPram.phase, &matrix<pt>::pPram.n, matrix<pt>::ent, matrix<pt>::ia, matrix<pt>::ja, &matrix<pt>::pPram.idum, &matrix<pt>::pPram.nrhs, matrix<pt>::pPram.iparm, &matrix<pt>::pPram.msglvl, rb, x, &matrix<pt>::pPram.error);
  if (matrix<pt>::pPram.error != 0) {
    printf("\nERROR during solution: %d", matrix<pt>::pPram.error);
    exit(3);
  }
#pragma omp parallel for
  for (int i = 0; i < size; ++i) {
    uvel->at(i)["sol"] = x[i];
    uvel->at(i)["phi"] = x[i + size];
    vvel->at(i)["sol"] = x[i + 2 * size];
    vvel->at(i)["phi"] = x[i + 3 * size];
    if (i < fixedPointIdx) {
      matrix<pt>::pts->at(i)["sol"] = x[i + 4 * size];
    } else if (i > fixedPointIdx) {
      matrix<pt>::pts->at(i)["sol"] = x[i + 4 * size - 1];
    } else {
      matrix<pt>::pts->at(i)["sol"] = ZEROVALUE;
    }
  }
  delete[] rb;
  delete[] rb0;
}

template<typename pt>
AGM::matrixStokesNormal<pt>::~matrixStokesNormal() {
  delete[] iaT;
  delete[] jaT;
  delete[] entT;
}

template class AGM::matrixStokesNormal<AGM::pointStokes>;

template class AGM::matrixStokesNormal<AGM::pointAxisymmetricStokes>;