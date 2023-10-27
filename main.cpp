#include "solver.h"
#include <sys/resource.h>

auto main() -> int {
    struct rlimit rlim{};
    getrlimit(RLIMIT_STACK, &rlim);
    rlim.rlim_cur = -1;
    rlim.rlim_max = -1;
    setrlimit(RLIMIT_STACK, &rlim);

    kmp_set_warnings_off();
    mkl_set_dynamic(0);
//    mkl_set_num_threads(4);
    mkl_set_num_threads(mkl_get_max_threads());

    omp_set_dynamic(0);
//    omp_set_num_threads(4);
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

    for (auto &item: pts) {
        item.setMp(UNITVALUE / 1000.);
//        item.setMp(UNITVALUE / 150.);
//        item.setMp(UNITVALUE / 1000.);
    }

    for (auto &item: pts) {
        if (item.getCondition() == 'C') {
            auto gap{std::vector<double>{
                    item - *(item[AGM::E]),
                    item - *(item[AGM::W]),
                    item - *(item[AGM::N]),
                    item - *(item[AGM::S])
            }};
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
    std::cout << "epsilon (Reynolds number) = " << pts.at(0).getMp() << "\n";
    std::cout << "Reynolds number = " << UNITVALUE / pts.at(0).getMp() << "\n";
    std::cout << "-------------------------" << "\n";
//    exit(1);
    /* Tesla valve */
//    pts.at(65779).setCondition('I');
//    pts.at(66981).setCondition('I');
    /* Tesla valve */

    /* export axial lines */ /*
    auto wf{AGM::writeFile<AGM::point>(&pts)};
    wf.writeAxialLines("/home/jhjo/extHD1/Navier-Stokes_Result/2D/0.Taylor-Green_vortex/point",
                       "/home/jhjo/extHD1/Navier-Stokes_Result/2D/0.Taylor-Green_vortex/xaxial",
                       "/home/jhjo/extHD1/Navier-Stokes_Result/2D/0.Taylor-Green_vortex/yaxial",
                       &xline,
                       &yline);
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

    /* Tesla valve */
    // forward
//    pts.at(65779).setCondition('N');
//    pts.at(66981).setCondition('N');
    // reverse
//    pts.at(65779).setCondition('D');
//    pts.at(66981).setCondition('D');
    /* Tesla valve */

    auto solver{AGM::solver(&pts)};
//    solver.ellipticSolver();
//    solver.streamSolver();
//    solver.FluidStructureInteraction();
    solver.NavierStokesSolver();
//    solver.TwoStepNS();

    return 0;
}