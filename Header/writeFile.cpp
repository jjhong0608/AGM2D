//
// Created by 조준홍 on 2022/02/19.
//

#include "writeFile.h"

template<typename T>
AGM::writeFile<T>::writeFile() = default;

template<typename T>
AGM::writeFile<T>::writeFile(std::vector<T> *pts) : pts(pts) {}

template<typename T>
auto AGM::writeFile<T>::getPt() const -> T * {
  return pt;
}

template<typename T>
void AGM::writeFile<T>::setPt(T *t) {
  writeFile::pt = t;
}

template<typename T>
auto AGM::writeFile<T>::getPts() const -> std::vector<T> * {
  return pts;
}

template<typename T>
void AGM::writeFile<T>::setPts(std::vector<T> *vector) {
  writeFile::pts = vector;
}

template<typename T>
auto AGM::writeFile<T>::calculateErrorAtPoint(const std::string &string) {
  auto printErrorToDouble = [&](const std::string &string1) -> double {
    std::cout << "\nidx = " << pt->getIdx() << ", (" << (*pt)[0] << ", " << (*pt)[1] << ")\n";
    printError("AGM::writeFile<T>::calculateErrorAtPoint()", string1.c_str());
    return ZEROVALUE;
  };
  double xm = (*pt)[W] ? (*pt)[W]->getXy()[0] : (*pt)[WN] ? (*pt)[WN]->getXy()[0]
                                                          : printErrorToDouble("xm");
  double xp = (*pt)[E] ? (*pt)[E]->getXy()[0] : (*pt)[EN] ? (*pt)[EN]->getXy()[0]
                                                          : printErrorToDouble("xp");
  double ym = (*pt)[S] ? (*pt)[S]->getXy()[1] : (*pt)[SE] ? (*pt)[SE]->getXy()[1]
                                                          : printErrorToDouble("ym");
  double yp = (*pt)[N] ? (*pt)[N]->getXy()[1] : (*pt)[NE] ? (*pt)[NE]->getXy()[1]
                                                          : printErrorToDouble("yp");
  auto g{AGM::NavierStokesFunction()};
  //  auto g{AGM::ellipticFunction()};
  //  pointHeat temp{};
  pointAxisymmetricStokes temp{};
  auto f = [&](const point &point) -> double {
    temp.point::operator=(point);
    //    return std::sqrt(g.u(pointHeat::getTime(), temp) + g.v(pointHeat::getTime(), temp));
    //    return g.phi(point);
    return g.v(pointHeat::getTime(), temp);
  };
  auto h{ZEROVALUE};
  //    if (pt->getCondition() == 'D') {
  //        if (isclose(pt->getXy()[1], ZEROVALUE)) {
  //            h = yp - ym;
  //            xm = ZEROVALUE;
  //            xp = h * 2;
  //        } else if (isclose(pt->getXy()[1], UNITVALUE)) {
  //            h = yp - ym;
  //            xm = ZEROVALUE;
  //            xp = h * 2;
  //        } else if (isclose(pt->getXy()[0], ZEROVALUE)) {
  //            h = xp - xm;
  //            ym = ZEROVALUE;
  //            yp = h * 2;
  //        } else if (isclose(pt->getXy()[0], UNITVALUE)) {
  //            h = xp - xm;
  //            ym = ZEROVALUE;
  //            yp = h * 2;
  //        }
  //    }
  double value = string == "grad" ? std::sqrt(std::pow((*pt)["dx"], 2) + std::pow((*pt)["dy"], 2)) : (*pt)[string];
  double numerator{std::pow(value - f(*pt), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
  double denominator{std::pow(f(*pt), 2) * (xp - xm) * (yp - ym) * 2.5E-1};
  return std::make_pair(numerator, denominator);
}

template<typename T>
auto AGM::writeFile<T>::calculateError(const std::string &string) -> double {
  double numerator{}, denominator{};
  auto error{std::pair<double, double>{}};
  for (auto &item : *pts) {
    pt = &item;
    error = calculateErrorAtPoint(string);
    numerator += error.first;
    denominator += error.second;
  }
  return std::sqrt(numerator / denominator);
}

template<typename T>
void AGM::writeFile<T>::writeResult(const std::string &string) {
  int bc{};
  std::ofstream f(string);
  if (!f.is_open()) {
    printError("AGM::writeFile<T>::writeResult", "file (%s) does not opened.", string.c_str());
  }
  f.precision(16);
  for (const auto &item : *pts) {
    bc = item.getCondition() == 'C' ? 0 : item.getCondition() == 'D' ? 1
        : item.getCondition() == 'N'                                 ? 2
        : item.getCondition() == 'I'                                 ? 3
                                                                     : 4;
    f << std::scientific;
    f << item.getIdx() << "\t";
    f << item[0] << "\t";
    f << item[1] << "\t";
    f << item["sol"] << "\t";
    f << item["phi"] << "\t";
    f << item["dx"] << "\t";
    f << item["dy"] << "\t";
    f << bc << "\n";
  }
  f.close();
}

template<typename T>
void AGM::writeFile<T>::writeAxialLines(const std::string &pname, const std::string &xname, const std::string &yname,
                                        std::vector<axialLine> *xline, std::vector<axialLine> *yline) {
  int bc{};
  std::ofstream p(pname);
  if (p.fail()) {
    printError("AGM::writeFile<T>::writeAxialLines", "file \"%s\" is not opened", pname.c_str());
  }
  p.precision(16);
  p << std::scientific;
  for (auto &item : *pts) {
    bc = item.getCondition() == 'C' ? 0 : 1;
    p << item[0] << "\t" << item[1] << "\t" << bc << "\n";
  }
  p.close();
  std::ofstream fx(xname);
  if (fx.fail()) {
    printError("AGM::writeFile<T>::writeAxialLines", "file \"%s\" is not opened", xname.c_str());
  }
  fx.precision(16);
  fx << std::scientific;
  for (auto &item : *xline) {
    fx << item[0] << "\t" << item[1] << "\t" << item[2] << "\t" << item[3] << "\n";
  }
  fx.close();
  std::ofstream fy(yname);
  if (fy.fail()) {
    printError("AGM::writeFile<T>::writeAxialLines", "file \"%s\" is not opened", xname.c_str());
  }
  fy.precision(16);
  fy << std::scientific;
  for (auto &item : *yline) {
    fy << item[0] << "\t" << item[1] << "\t" << item[2] << "\t" << item[3] << "\n";
  }
  fy.close();
}

template class AGM::writeFile<AGM::point>;

template class AGM::writeFile<AGM::pointHeat>;

template class AGM::writeFile<AGM::pointAxisymmetric>;

template class AGM::writeFile<AGM::pointAxisymmetricStokes>;
