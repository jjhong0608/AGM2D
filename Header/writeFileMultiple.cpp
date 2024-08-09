//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#include "writeFileMultiple.h"

#include "Eigen/StdVector"

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::writeFileMultiple() = default;

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::writeFileMultiple(std::vector<T0> *pts0, std::vector<T1> *pts1,
                                                      std::vector<T2> *pts2) : pts0(pts0), pts1(pts1), pts2(pts2) {}

template<typename T0, typename T1, typename T2>
auto AGM::writeFileMultiple<T0, T1, T2>::calculateErrorAtPoint(const std::string &string) {
  auto printErrorToDouble = [&](const std::string &string1) -> double {
    std::cout << "\nidx = " << pt0->getIdx() << ", (" << (*pt0)[0] << ", " << (*pt0)[1] << ")\n";
    printError("AGM::writeFIleMultiple<T0, T1, T2>::calculateErrorAtPoint()", string1.c_str());
    return ZEROVALUE;
  };
  auto xm{(*pt0)[W] ? (*pt0)[W]->getXy()[0] : (*pt0)[WN] ? (*pt0)[WN]->getXy()[0]
                                                         : printErrorToDouble("xm")};
  auto xp{(*pt0)[E] ? (*pt0)[E]->getXy()[0] : (*pt0)[EN] ? (*pt0)[EN]->getXy()[0]
                                                         : printErrorToDouble("xp")};
  auto ym{(*pt0)[S] ? (*pt0)[S]->getXy()[1] : (*pt0)[SE] ? (*pt0)[SE]->getXy()[1]
                                                         : printErrorToDouble("ym")};
  auto yp{(*pt0)[N] ? (*pt0)[N]->getXy()[1] : (*pt0)[NE] ? (*pt0)[NE]->getXy()[1]
                                                         : printErrorToDouble("yp")};
  auto g{AGM::NavierStokesFunction()};
  auto f = [&](const T0 &point) {
    return std::make_pair(g.u(pointHeat::getTime(), point), g.v(pointHeat::getTime(), point));
  };
  auto velocity0{std::sqrt(std::pow(f(*pt0).first, 2.0) + std::pow(f(*pt1).second, 2.0))};
  auto velocity1{std::sqrt(std::pow((*pt0)["sol"], 2.0) + std::pow((*pt1)["sol"], 2.0))};
  auto numerator{std::pow(velocity1 - velocity0, 2.0) * (xp - xm) * (yp - ym) * 0.25};
  auto denominator{std::pow(velocity0, 2.0) * (xp - xm) * (yp - ym) * 0.25};
  return std::make_pair(numerator, denominator);
}
template<typename T0, typename T1, typename T2>
auto AGM::writeFileMultiple<T0, T1, T2>::calculateError(const std::string &string) -> double {
  auto numerator{ZEROVALUE}, denominator{ZEROVALUE};
  auto error{std::pair<double, double>{}};
  for (int i = 0; i < point::getNPts(); ++i) {
    pt0 = &pts0->at(i);
    pt1 = &pts1->at(i);
    pt2 = &pts2->at(i);
    error = calculateErrorAtPoint(string);
    numerator += error.first;
    denominator += error.second;
  }
  return std::sqrt(numerator / denominator);
}

template<typename T0, typename T1, typename T2>
void AGM::writeFileMultiple<T0, T1, T2>::writeResult(const std::string &string) {
  int bc{};
  std::ofstream f(string);
  if (f.fail()) {
    printError("AGM::writeFileMultiple<T0, T1, T2>::writeResult", "file (%s) does not opened.", string.c_str());
  }
  f.precision(16);
  for (int i = 0; i < point::getNPts(); ++i) {
    bc = pts2->at(i).getCondition() == 'C' ? 0 : pts2->at(i).getCondition() == 'N' ? 1
        : pts2->at(i).getCondition() == 'D'                                        ? 2
        : pts2->at(i).getCondition() == 'I'                                        ? 3
                                                                                   : 4;
    f << std::scientific;
    f << i << "\t";                 // idx
    f << pts0->at(i)[0] << "\t";    // x
    f << pts0->at(i)[1] << "\t";    // y
    f << pts0->at(i)["sol"] << "\t";// u
    f << pts1->at(i)["sol"] << "\t";// v
    f << pts2->at(i)["sol"] << "\t";// p
    f << pts0->at(i)["dx"] << "\t"; // ux
    f << pts0->at(i)["dy"] << "\t"; // uy
    f << pts1->at(i)["dx"] << "\t"; // vx
    f << pts1->at(i)["dy"] << "\t"; // vy
    f << pts2->at(i)["dx"] << "\t"; // px
    f << pts2->at(i)["dy"] << "\t"; // py
    f << pts0->at(i)["phi"] << "\t";// phi
    f << pts1->at(i)["phi"] << "\t";// psi
    f << pts2->at(i)["phi"] << "\t";// phi of p
    f << bc << "\n";
  }
  f.close();
}

template<typename T0, typename T1, typename T2>
void AGM::writeFileMultiple<T0, T1, T2>::writeStruct(const std::string &string) {
  int bc{};
  std::ofstream f(string);
  if (f.fail()) {
    printError("AGM::writeFileMultiple<T0, T1, T2>::writeStruct", "file (%s) does not opened.", string.c_str());
  }
  f.precision(16);
  for (int i = 0; i < pts0->size(); ++i) {
    f << std::scientific;
    f << i << "\t";                 // idx
    f << pts0->at(i)[0] << "\t";    // x
    f << pts0->at(i)[1] << "\t";    // y
    f << pts0->at(i)["sol"] << "\t";// u
    f << pts1->at(i)["sol"] << "\n";// v
  }
  f.close();
}

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::~writeFileMultiple() = default;

template class AGM::writeFileMultiple<AGM::pointStokes, AGM::pointStokes, AGM::pointStokes>;

template class AGM::writeFileMultiple<AGM::pointAxisymmetricStokes, AGM::pointAxisymmetricStokes, AGM::pointAxisymmetricStokes>;

template class AGM::writeFileMultiple<AGM::pointHeat, AGM::pointHeat, AGM::point>;

template class AGM::writeFileMultiple<AGM::point, AGM::point, AGM::point>;

//template class AGM::writeFileMultiple<AGM::structure, AGM::structure, AGM::point>;
