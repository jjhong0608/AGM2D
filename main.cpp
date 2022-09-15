#include "solver.h"
#include <sys/resource.h>

auto main() -> int {
    struct rlimit rlim{};
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = -1;
    rlim.rlim_max = -1;
    setrlimit(RLIMIT_STACK, &rlim);

    auto pts{std::vector<AGM::point>{}};
    auto xline{std::vector<AGM::axialLine>{}};
    auto yline{std::vector<AGM::axialLine>{}};
    auto bdline{std::vector<AGM::boundaryLine2D>{}};
    AGM::readFile::loadData("ALG_output", &pts, &xline, &yline);
    AGM::readFile::loadBoundaryData("GeoInfo", &bdline);
    AGM::point::setAxialLines(&xline, 'x');
    AGM::point::setAxialLines(&yline, 'y');
    AGM::point::setBdLine(&bdline);

    /* export axial lines */ /*
    auto wf{AGM::writeFile<AGM::point>(&pts)};
    wf.writeAxialLines("/home/jjhong0608/docker/ALG2D/point", "/home/jjhong0608/docker/ALG2D/xaxial",
                       "/home/jjhong0608/docker/ALG2D/yaxial", &xline, &yline);
    return 0;
    */

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

    std::cout << "Total Points number = " << pts.size() << "\n";
    std::cout << "Reynols number = " << UNITVALUE / pts[0].getMp() << "\n";

    auto solver{AGM::solver(&pts)};
    solver.streamSolver();
//    solver.NavierStokesSolver();

    return 0;
}
