//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_VALUE_H
#define AGM_VALUE_H

#include "coordinate.h"

namespace AGM {
class value : public std::array<double, 8> {

 public:
  value();

  virtual ~value();

  auto operator[](const std::string &string) -> double &;

  auto operator[](const std::string &string) const -> const double &;
};
}// namespace AGM

#endif//AGM_VALUE_H
