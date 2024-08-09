//
// Created by 조준홍 on 1/27/24.
//

#ifndef AGM_MATRIXSTOKESNORMAL_H
#define AGM_MATRIXSTOKESNORMAL_H

#include "matrixPhi.h"

namespace AGM {
template<typename pt>
class matrixStokesNormal : public matrix<pt> {
 public:
  matrixStokesNormal();

  matrixStokesNormal(std::vector<pt> *uvel, std::vector<pt> *vvel, std::vector<pt> *pres, int fixedPointIdx);

  virtual ~matrixStokesNormal();

  [[nodiscard]] auto getFixedPointIdx() const -> int;

  void setFixedPointIdx(int fixed_point_idx);

  void makeMatrix() override;

  void factorizeMatrix() override;

  void calculateMatrix() override;

 private:
  Eigen::SparseMatrix<double, Eigen::RowMajor> A{};
  std::vector<pt> *uvel{}, *vvel{};
  int fixedPointIdx{-1};
  int *iaT{}, *jaT{};
  double *entT{};
};
}// namespace AGM

#endif//AGM_MATRIXSTOKESNORMAL_H
