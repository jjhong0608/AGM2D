#include "solver.h"

int main() {
    auto pts{std::vector<AGM::point>{}};
    auto xline{std::vector<AGM::axialLine>{}};
    auto yline{std::vector<AGM::axialLine>{}};
    auto bdline{std::vector<AGM::boundaryLine2D>{}};
    AGM::readFile::loadData("ALG_output", &pts, &xline, &yline);
    AGM::readFile::loadBoundaryData("GeoInfo", &bdline);
    AGM::point::setAxialLines(&xline, 'x');
    AGM::point::setAxialLines(&yline, 'y');
    AGM::point::setBdLine(&bdline);

    for (auto &item: pts) {
        AGM::point::setPts(&pts);
        item.findStencil();
    }
    for (auto &item: pts) {
        if (item.getCondition() == 'd' || item.getCondition() == 'n') {
            item.findStencil();
            std::cout << "condition = " << item.getCondition() << "\n";
        }
    }
    auto solver{AGM::solver(&pts)};
    solver.NavierStokesSolver();
//    auto wf{AGM::writeFile<AGM::point>(&pts)};
//    std::cout << "Relative L-2 Error = " << wf.calculateError("sol") << "\n";
//    wf.writeResult("/home/jjhong0608/docker/AGM2D/air_foil/adaptive/AGM_Result");

//    AGM::point::setNPts(int(pts.size()));
//    auto wf{AGM::writeFileMultiple<AGM::point, AGM::point, AGM::point>(&pts, &pts, &pts)};
//    wf.writeResult("/home/jjhong0608/docker/AGM2D/air_foil/adaptive/AGM_Result");

    return 0;
}
