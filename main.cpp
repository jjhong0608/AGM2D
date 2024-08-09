#include <sys/resource.h>

#include "Header/solver.h"

auto main() -> int {
  struct rlimit rlim {};
  getrlimit(RLIMIT_STACK, &rlim);
  rlim.rlim_cur = -1;
  rlim.rlim_max = -1;
  setrlimit(RLIMIT_STACK, &rlim);

  kmp_set_warnings_off();
  mkl_set_dynamic(0);
  mkl_set_num_threads(int(mkl_get_max_threads() / 2));

  omp_set_dynamic(0);
  omp_set_num_threads(int(omp_get_max_threads() / 2));

  auto pts{std::vector<AGM::point>{}};
  auto xline{std::vector<AGM::axialLine>{}};
  auto yline{std::vector<AGM::axialLine>{}};
  auto bdline{std::vector<AGM::boundaryLine2D>{}};
  AGM::readFile::loadData("ALG_output", &pts, &xline, &yline);
  AGM::readFile::loadBoundaryData("GeoInfo", &bdline);
  AGM::point::setAxialLines(&xline, 'x');
  AGM::point::setAxialLines(&yline, 'y');
  AGM::point::setBdLine(&bdline);

  for (auto &item : pts) {
    if (item.getCondition() == 'C') {
      auto gap{std::vector<double>{item - *(item[AGM::E]), item - *(item[AGM::W]), item - *(item[AGM::N]), item - *(item[AGM::S])}};
      std::sort(gap.begin(), gap.end());
      if (gap.back() > AGM::point::getAlinMaxGap()) {
        AGM::point::setAlinMaxGap(gap.back());
      }
    }
  }
  std::cout << "-----< information >-----" << "\n";
  std::cout << "# of the points = " << pts.size() << "\n";
  std::cout << "# of the x-axial lines = " << xline.size() << "\n";
  std::cout << "# of the y-axial lines = " << yline.size() << "\n";
  std::cout << "epsilon (1 / Reynolds number) = " << pts.at(0).getMp() << "\n";
  std::cout << "Reynolds number = " << UNITVALUE / pts.at(0).getMp() << "\n";
  std::cout << "-------------------------" << "\n";

  for (auto &item : pts) {
    AGM::point::setPts(&pts);
    item.findStencil();
  }
  for (auto &item : pts) {
    if (item.getCondition() == 'd' || item.getCondition() == 'n') {
      item.findStencil();
      std::cout << "condition = " << item.getCondition() << "\n";
    }
  }

  auto solver{AGM::solver(&pts)};
  //  solver.ellipticSolver();
  //  solver.streamSolver();
  //  solver.FluidStructureInteraction();
  //  solver.StokesSolver();
  //  solver.StokesSolverFull();
  solver.NavierStokesSolver();
  //  solver.axisymmetricStokesSolverFull();
  //  solver.axisymmetricStokesSolver();

  return 0;
}