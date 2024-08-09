//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#ifndef AGM_WRITEFILEMULTIPLE_H
#define AGM_WRITEFILEMULTIPLE_H

#include "Eigen/StdVector"
#include "writeFile.h"

namespace AGM {
template<typename T0, typename T1, typename T2>
class writeFileMultiple {
 public:
  writeFileMultiple();

  virtual ~writeFileMultiple();

  writeFileMultiple(std::vector<T0> *pts0, std::vector<T1> *pts1, std::vector<T2> *pts2);

  auto calculateErrorAtPoint(const std::string &string);

  auto calculateError(const std::string &string) -> double;

  void writeResult(const std::string &string);

  void writeStruct(const std::string &string);

 private:
  T0 *pt0{};
  T1 *pt1{};
  T2 *pt2{};
  std::vector<T0> *pts0{};
  std::vector<T1> *pts1{};
  std::vector<T2> *pts2{};
};

}// namespace AGM

#endif//AGM_WRITEFILEMULTIPLE_H
