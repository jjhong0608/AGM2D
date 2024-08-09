//
// Created by 조준홍 on 2023/08/22.
//

#ifndef AGM_MATRIXPHI_H
#define AGM_MATRIXPHI_H

#include "matrixNormal.h"

namespace AGM {
template<typename pt>
class matrixPhi : public matrix<pt> {
 public:
  matrixPhi();

  matrixPhi(std::vector<pt> *pts);

  virtual ~matrixPhi();

  void makeMatrix() override;

  void factorizeMatrix() override;

  void calculateMatrix() override;
};
}// namespace AGM

#endif//AGM_MATRIXPHI_H
