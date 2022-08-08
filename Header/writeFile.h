//
// Created by 조준홍 on 2022/02/19.
//

#ifndef AGM_WRITEFILE_H
#define AGM_WRITEFILE_H

#include "readFile.h"

namespace AGM {
    template<typename T>
    class writeFile {
    protected:
        T *pt{};
        std::vector<T> *pts{};

    public:
        writeFile();

        explicit writeFile(std::vector<T> *pts);

        auto getPt() const -> T *;

        void setPt(T *t);

        auto getPts() const -> std::vector<T> *;

        void setPts(std::vector<T> *vector);

        auto calculateErrorAtPoint(const std::string &string);

        auto calculateError(const std::string &string) -> double;

        void writeResult(const std::string &string);

        void writeAxialLines(const std::string &pname, const std::string &xname, const std::string &yname,
                             std::vector<axialLine> *xline, std::vector<axialLine> *yline);
    };

}


#endif //AGM_WRITEFILE_H
