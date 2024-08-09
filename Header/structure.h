//
// Created by 조준홍 on 2023/01/04.
//

#ifndef AGM_STRUCTURE_H
#define AGM_STRUCTURE_H

#include "pointAxisymmetricStokes.h"

namespace AGM {
auto deltaFtn(double r, double h) -> double;

class structure {
 public:
  structure();

  explicit structure(int idx);

  structure(int idx, const coordinate &xy);

  virtual ~structure();

  [[nodiscard]] auto getIdx() const -> int;

  void setIdx(int i);

  [[nodiscard]] auto getXy() const -> const coordinate &;

  void setXy(const coordinate &coordinate);

  [[nodiscard]] auto getValues() const -> const value &;

  void setValues(const value &value);

  static auto getNStructures() -> int;

  static void setNStructures(int i);

  static auto getHf() -> double;

  static void setHf(double d);

  static auto getHs() -> double;

  static void setHs(double d);

  static auto getStructures() -> std::array<std::vector<structure> *, 2> *;

  static void setStructures(std::array<std::vector<structure> *, 2> *array);

  auto operator-(const point &src) -> double;

  auto operator-(const pointHeat &src) -> double;

  auto operator-(const structure &src) -> double;

  auto operator[](int i) -> double &;

  auto operator[](int i) const -> const double &;

  auto operator[](const std::string &string) -> double &;

  auto operator[](const std::string &string) const -> const double &;

  void copyFluidVelocity(pointHeat *src);

  void structureForceUpdateToFluid(pointHeat *src);

  void updateForce(std::vector<structure> *previousStructures, int component);

 private:
  int idx{};
  coordinate xy{};
  value values{};
  static int nStructures;
  static double hf, hs;
  static std::array<std::vector<structure> *, 2> *structures;
};

}// namespace AGM

#endif//AGM_STRUCTURE_H
