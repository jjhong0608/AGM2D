//
// Created by NIMS-JUNHONG on 2022/01/21.
//

#ifndef AGM_AXIALLINE_H
#define AGM_AXIALLINE_H

#include "GreenfunctionReactionDiffusion.h"

namespace AGM {
    class point;

    class axialLine : public std::vector<point *> {
    private:
        char mark{};
        std::array<double, 4> coordinate{};

    public:
        axialLine();

        explicit axialLine(char mark);

        virtual ~axialLine();

        [[nodiscard]] auto getMark() const -> char;

        void setMark(char i);

        auto operator[](int i) -> double &;

        auto operator-(axialLine &line) -> double;
    };
}


#endif //AGM_AXIALLINE_H
