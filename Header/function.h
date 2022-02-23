//
// Created by 조준홍 on 2022/02/13.
//

#ifndef AGM_FUNCTION_H
#define AGM_FUNCTION_H

#include "pointHeat.h"

namespace AGM {
    class ellipticFunction {
    public:
        ellipticFunction();

        virtual ~ellipticFunction();

        double u(const point &pt);

        double phi(const point &pt);

        double f(const point &pt);

        double ux(const point &pt);

        double uy(const point &pt);

        void assignBoundaryValue(point &pt);
    };

    class heatFunction {
    public:
        heatFunction();

        virtual ~heatFunction();

        static double initialTime();

        static double terminalTime();

        static double deltaTime();

        double u(double t, const point &pt);

        double phi(double t, const point &pt);

        double f(double t, const point &pt);

        double ux(double t, const point &pt);

        double uy(double t, const point &pt);

        void assignPreviousValue(AGM::value &value, point &pt);

        void assignBoundaryValue(point &pt);
    };

    class NavierStokesFunction {
    public:
        NavierStokesFunction();

        virtual ~NavierStokesFunction();

        static double initialTime();

        static double terminalTime();

        static double deltaTime();

        double u(double t, const point &pt);

        double v(double t, const point &pt);

        double p(double t, const point &pt);

        double phi(double t, const point &pt);

        double psi(double t, const point &pt);

        double ux(double t, const point &pt);

        double uy(double t, const point &pt);

        double vx(double t, const point &pt);

        double vy(double t, const point &pt);

        double px(double t, const point &pt);

        double py(double t, const point &pt);

        double f1(double t, const point &pt);

        double f2(double t, const point &pt);

        void assignPreviousValue(AGM::value &pu, AGM::value &pv, AGM::value &pp, point &uvel, point &vvel, point &pres);

        void assignBoundaryValue(point &uvel, point &vvel);
    };
}


#endif //AGM_FUNCTION_H
