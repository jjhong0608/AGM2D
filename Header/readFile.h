//
// Created by 조준홍 on 2022/02/18.
//

#ifndef AGM_READFILE_H
#define AGM_READFILE_H

#include "matrixNormal.h"

namespace AGM {
    class readFile {
    public:
        static void loadData(const std::string &filename, std::vector<point> *pts, std::vector<axialLine> *alineX,
                             std::vector<axialLine> *alineY);

        static void loadBoundaryData(const std::string &filename, std::vector<boundaryLine2D> *bdLine);
    };

}


#endif //AGM_READFILE_H
