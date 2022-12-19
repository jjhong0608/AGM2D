#include "solver.h"
#include <sys/resource.h>

auto main() -> int {
    struct rlimit rlim{};
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = -1;
    rlim.rlim_max = -1;
    setrlimit(RLIMIT_STACK, &rlim);

    mkl_set_dynamic(0);
//    mkl_set_num_threads(2);
    mkl_set_num_threads(mkl_get_max_threads());

    omp_set_dynamic(0);
//    omp_set_num_threads(2);
    omp_set_num_threads(omp_get_max_threads());

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
    wf.writeAxialLines("/home/jjhong0608/docker/Navier-Stokes_Result/2D/3.External_flow_past_circular_cylinder/point",
                       "/home/jjhong0608/docker/Navier-Stokes_Result/2D/3.External_flow_past_circular_cylinder/xaxial",
                       "/home/jjhong0608/docker/Navier-Stokes_Result/2D/3.External_flow_past_circular_cylinder/yaxial",
                       &xline,
                       &yline);
    return 0;
    */

    for (auto &item: pts) {
        AGM::point::setPts(&pts);
        item.findStencil();
    }
//    for (auto &item: pts) {
//        item.setMp(UNITVALUE / 1e1);
//    }
    for (auto &item: pts) {
        if (item.getCondition() == 'd' || item.getCondition() == 'n') {
            item.findStencil();
            std::cout << "condition = " << item.getCondition() << "\n";
        }
    }

    std::cout << "-----< information >-----" << "\n";
    std::cout << "# of the points = " << pts.size() << "\n";
    std::cout << "# of the x-axial lines = " << xline.size() << "\n";
    std::cout << "# of the y-axial lines = " << yline.size() << "\n";
    std::cout << "epsilon (Reynolds number) = " << pts.at(0).getMp() << "\n";
    std::cout << "Reynolds number = " << UNITVALUE / pts.at(0).getMp() << "\n";
    std::cout << "-------------------------" << "\n";

    auto solver{AGM::solver(&pts)};
//    solver.streamSolver();
    solver.NavierStokesSolver();

    return 0;
}
