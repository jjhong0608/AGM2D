//
// Created by NIMS-JUNHONG on 2022/02/23.
//

#include "writeFileMultiple.h"
#include "StdVector"

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::writeFileMultiple() = default;

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::writeFileMultiple(std::vector<T0> *pts0, std::vector<T1> *pts1,
                                                      std::vector<T2> *pts2):pts0(pts0), pts1(pts1), pts2(pts2) {}

template<typename T0, typename T1, typename T2>
void AGM::writeFileMultiple<T0, T1, T2>::writeResult(const std::string &string) {
    int bc{};
    std::ofstream f(string);
    if (f.fail()) {
        printError("AGM::writeFileMultiple<T0, T1, T2>::writeResult", "file (%s) does not opened.", string.c_str());
    }
    f.precision(16);
    for (int i = 0; i < point::getNPts(); ++i) {
        bc = pts2->at(i).getCondition() == 'C' ? 0 : pts2->at(i).getCondition() == 'N' ? 1 :
                                                     pts2->at(i).getCondition() == 'D' ? 2 :
                                                     pts2->at(i).getCondition() == 'I' ? 3 : 4;
        f << std::scientific;
        f << i << "\t";                    // idx
        f << pts0->at(i)[0] << "\t";       // x
        f << pts0->at(i)[1] << "\t";       // y
        f << pts0->at(i)["sol"] << "\t";   // u
        f << pts1->at(i)["sol"] << "\t";   // v
        f << pts2->at(i)["sol"] << "\t";   // p
        f << pts0->at(i)["dx"] << "\t";    // ux
        f << pts0->at(i)["dy"] << "\t";    // uy
        f << pts1->at(i)["dx"] << "\t";    // vx
        f << pts1->at(i)["dy"] << "\t";    // vy
        f << pts2->at(i)["dx"] << "\t";   // px
        f << pts2->at(i)["dy"] << "\t";   // py
        f << pts0->at(i)["phi"] << "\t";   // phi
        f << pts1->at(i)["phi"] << "\t";   // psi
        f << pts2->at(i)["phi"] << "\t";   // phi of p
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
        f << i << "\t";                    // idx
        f << pts0->at(i)[0] << "\t";       // x
        f << pts0->at(i)[1] << "\t";       // y
        f << pts0->at(i)["sol"] << "\t";   // u
        f << pts1->at(i)["sol"] << "\n";   // v
    }
    f.close();
}

template<typename T0, typename T1, typename T2>
AGM::writeFileMultiple<T0, T1, T2>::~writeFileMultiple() = default;

template
class AGM::writeFileMultiple<AGM::pointHeat, AGM::pointHeat, AGM::point>;

template
class AGM::writeFileMultiple<AGM::point, AGM::point, AGM::point>;

template
class AGM::writeFileMultiple<AGM::structure, AGM::structure, AGM::point>;
