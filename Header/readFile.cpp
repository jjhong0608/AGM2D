//
// Created by 조준홍 on 2022/02/18.
//

#include "readFile.h"

void AGM::readFile::loadData(const std::string &filename, std::vector<point> *pts, std::vector<axialLine> *alineX,
                             std::vector<axialLine> *alineY) {
    std::ifstream f(filename);
    std::string line{}, tempString{};
    int nRegion{}, idx{}, index{}, tempInteger{}, nAxial{};
    int nCross{}, nBound{}, nXAxial{}, nYAxial{};
    coordinate normal{};
    double mp{}, bv{};
    char bc{};
    point pt{};
    std::vector<std::vector<int>> xline{}, yline{};
    std::vector<int> aline{};

    if (!f.is_open()) {
        printError("AGM::readFile::loadData", "No Axial data file: \"%s\"\nPlease check file name", filename.c_str());
    }

    std::cout << "Axial file: \"" << filename << "\" open\n";

    while (!f.eof()) {
        std::getline(f, line);
        if (line.empty()) std::getline(f, line);
        if (line.find("ENDREGION") != std::string::npos) {
            idx = 0;
            std::getline(f, line);
        }
        if (line.find("REGION") != std::string::npos) {
            ++nRegion;
            std::cout << "REGION " << nRegion << " reading...\n";
            std::getline(f, line);
            f >> tempString >> tempString >> tempString >> mp;
        }
        if (line.find('#') != std::string::npos) {
            ++idx;
            switch (idx) {
                case 1:
                    nCross = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nCross; ++i) {
                        f >> tempInteger >> pt[0] >> pt[1];
                        pt.setMp(mp);
                        pt.setCondition('C');
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 2:
                    nBound = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    for (int i = 0; i < nBound; ++i) {
                        f >> tempInteger >> pt[0] >> pt[1] >> bc >> bv >> normal[0] >> normal[1];
                        pt.setMp(mp);
                        if (bc != 'D' && bc != 'N' && bc != 'I' && bc != 'F') {
                            printError("AGM::ReadFile::loadAxialData", "boundary condition (which is %s) is wrong", bc);
                        }
                        pt.setCondition(bc);
                        pt["bdv"] = bv;
                        pt.setNormal(normal);
                        pt.setIdx(index++);
                        pts->emplace_back(pt);
                    }
                    break;
                case 3:
                    nXAxial = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nAxial = 0;
                    while (nAxial != nXAxial) {
                        f >> tempString;
                        while (tempString != "/") {
                            aline.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        xline.emplace_back(aline);
                        aline.clear();
                        ++nAxial;
                    }
                    break;
                case 4:
                    nYAxial = std::stoi(line.substr(line.find('=') + 2, line.size()));
                    nAxial = 0;
                    while (nAxial != nYAxial) {
                        f >> tempString;
                        while (tempString != "/") {
                            aline.emplace_back(std::stoi(tempString));
                            f >> tempString;
                        }
                        yline.emplace_back(aline);
                        aline.clear();
                        ++nAxial;
                    }
                    break;
                default:
                    printError("AGM::ReadFile::loadAxialData", "idx (which is %s) is wrong", idx);
            }
        }
    }
    f.close();

    auto xAxialLine = axialLine('x');
    auto yAxialLine = axialLine('y');

    for (const auto &item: xline) {
        xAxialLine[0] = pts->at(item.front())[0];
        xAxialLine[1] = pts->at(item.back())[0];
        xAxialLine[2] = pts->at(item.front())[1];
        xAxialLine[3] = pts->at(item.back())[1];
        for (const auto &item0: item) {
            xAxialLine.emplace_back(&(pts->at(item0)));
        }
        alineX->emplace_back(xAxialLine);
        xAxialLine.clear();
    }
    for (const auto &item: yline) {
        yAxialLine[0] = pts->at(item.front())[0];
        yAxialLine[1] = pts->at(item.back())[0];
        yAxialLine[2] = pts->at(item.front())[1];
        yAxialLine[3] = pts->at(item.back())[1];
        for (const auto &item0: item) {
            yAxialLine.emplace_back(&(pts->at(item0)));
        }
        alineY->emplace_back(yAxialLine);
        yAxialLine.clear();
    }
    for (auto &item: *alineX) {
        for (const auto &item0: item) {
            item0->setAxialLine(&item, 'x');
        }
    }
    for (auto &item: *alineY) {
        for (const auto &item0: item) {
            item0->setAxialLine(&item, 'y');
        }
    }

    for (auto &item: *alineX) {
        for (auto &item0: item) {
            if (item0 == item.front()) {
                if (item0->getCondition() != 'I') {
                    (*item0)[W] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[N] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[S] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[N] = (*item0)[S] = item0;
                }
            } else if (item0 == item.back()) {
                (**std::prev(&item0))[E] = item0;
                if (item0->getCondition() != 'I') {
                    (*item0)[E] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[N] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[S] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[N] = (*item0)[S] = item0;
                }
                (*item0)[W] = *std::prev(&item0);
            } else {
                (**std::prev(&item0))[E] = item0;
                (*item0)[W] = *std::prev(&item0);
            }
        }
    }
    for (auto &item: *alineY) {
        for (auto &item0: item) {
            if (item0 == item.front()) {
                if (item0->getCondition() != 'I') {
                    (*item0)[S] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[E] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[W] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[E] = (*item0)[W] = item0;
                }
            } else if (item0 == item.back()) {
                (**std::prev(&item0))[N] = item0;
                if (item0->getCondition() != 'I') {
                    (*item0)[N] = item0;
                    if (ispositive(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[E] = item0;
                    if (isnegative(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[W] = item0;
                    if (iszero(item0->getNormal()[1]) && !item0->getAxialLine('y')) (*item0)[E] = (*item0)[W] = item0;
                }
                (*item0)[S] = *std::prev(&item0);
            } else {
                (**std::prev(&item0))[N] = item0;
                (*item0)[S] = *std::prev(&item0);
            }
        }
    }
}

void AGM::readFile::loadBoundaryData(const std::string &filename, std::vector<AGM::boundaryLine2D> *bdLine) {
    std::ifstream f(filename);
    std::string word{}, name{};
    auto readSectionLine = [&f, &bdLine]() -> void {
        std::string word{}, name{}, sx{}, sy{}, ex{}, ey{}, bc{}, bv{};
        auto start{vector{}}, end{vector{}};
        auto line{boundaryLine2D{}};

        start.resize(2);
        end.resize(2);

        f >> name;
        while (word != "ENDSECTION") {
            if (word == "LINE") {
                f >> sx >> sy >> ex >> ey >> bc;
                if (bc != "I") f >> bv;
                start[0] = std::stod(sx);
                start[1] = std::stod(sy);
                end[0] = std::stod(ex);
                end[1] = std::stod(ey);
                line.setStart(start);
                line.setAnEnd(end);
                line.setCondition(bc.c_str()[0]);
                if (line.getCondition() != 'I') line.setBoundaryValue(std::stod(bv));
                line.calcProperties();
                bdLine->emplace_back(line);
            }
            f >> word;
        }
    };
    if (!f.is_open()) {
        printError("AGM::readFile::loadBoundaryData", "No Boundary Data file: \"%s\"\n Please check file name",
                   filename.c_str());
    }
    std::cout << "Boundary Data file: \"" << filename << "\" open" << "\n";
    while (!f.eof()) {
        f >> word;
        if (word == "SECTION") {
            readSectionLine();
        }
        word = "";
    }
    f.close();
}
