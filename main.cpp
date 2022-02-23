#include "Header/solver.h"

int main() {
    auto pts{std::vector<AGM::point>{}};
    auto xline{std::vector<AGM::axialLine>{}};
    auto yline{std::vector<AGM::axialLine>{}};
    AGM::readFile::loadData("ALG_output", &pts, &xline, &yline);
    AGM::point::setAxialLines(&xline, 'x');
    AGM::point::setAxialLines(&yline, 'y');

    for (auto &item: pts) {
        item.findStencil();
    }
    auto solver{AGM::solver(&pts)};
    solver.NavierStokesSolver();
//    auto wf{AGM::writeFile<AGM::point>(&pts)};
//    std::cout << "Relative L-2 Error = " << wf.calculateError("sol") << "\n";
//    wf.writeResult("/home/jjhong0608/docker/AGM2D/AGM_Result");

    return 0;
}
