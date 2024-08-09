//
// Created by 조준홍 on 12/3/23.
//

#include "GreenfunctionAxisymmetricStokes.h"

AGM::GreenfunctionAxisymmetricStokes::GreenfunctionAxisymmetricStokes(double tm, double tau, double tp, double mpl, double mpr) : Greenfunction(tm, tau, tp, mpl, mpr) {
  a1 = tau * (tp * tp - tau * tau) / ((tau * tau - tm * tm) * (tp * tp + tau * tau) + (tau * tau + tm * tm) * (tp * tp - tau * tau));
  a2 = tau * (tau * tau - tm * tm) / ((tau * tau - tm * tm) * (tp * tp + tau * tau) + (tau * tau + tm * tm) * (tp * tp - tau * tau));
  b1 = tm * tm;
  b2 = tp * tp;
  c1 = -(tp * tp + tau * tau) / (2. * tau * tau * (tp * tp - tm * tm));
  c2 = (tau * tau + tm * tm) / (2. * tau * tau * (tp * tp - tm * tm));
  d1 = a1 * b1;
  d2 = a2 * b2;
  e1 = c1 * b1;
  e2 = c2 * b2;
}

auto AGM::GreenfunctionAxisymmetricStokes::green_function(double t) const -> double {
  if (t < tau || isclose(tm, t)) {
    return a1 * (t - b1 / t);
  } else {
    return a2 * (b2 / t - t);
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::green_function_t(double t) const -> double {
  if (t < tau || isclose(tm, t)) {
    if (iszero(tm)) {
      return 2. * a1;
    }
    return a1 * (UNITVALUE + b1 / std::pow(t, 2.));
  } else {
    return -a2 * (UNITVALUE + b2 / std::pow(t, 2.));
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::green_function_tau(double t) const -> double {
  if (t < tau || isclose(tm, t)) {
    return c1 * (t - b1 / t);
  } else {
    return c2 * (b2 / t - t);
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::green_function_ttau(double t) const -> double {
  if (t < tau || isclose(tm, t)) {
    if (iszero(tm)) {
      return 1. * c1;
    }
    return c1 * (UNITVALUE + b1 / std::pow(t, 2.));
  } else {
    return -c2 * (UNITVALUE + b2 / std::pow(t, 2.));
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_square(char i) const -> double {
  switch (i) {
    case 'l': return a1 / 4. * (std::pow(tau, 4.) - std::pow(tm, 4.)) - d1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.));
    case 'r': return a2 / 4. * (std::pow(tau, 4.) - std::pow(tp, 4.)) - d2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_linear(char i) const -> double {
  switch (i) {
    case 'l': return a1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) - d1 * (tau - tm);
    case 'r': return a2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) - d2 * (tau - tp);
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linear", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_const(char i) const -> double {
  switch (i) {
    case 'l':
      if (iszero(tm)) {
        return a1 / 2. * std::pow(tau, 2.) - d1 * std::log(tau);
      } else {
        return a1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) - d1 * (std::log(tau) - std::log(tm));
      }
    case 'r': return a2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) - d2 * (std::log(tau) - std::log(tp));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_square_t(char i) const -> double {
  switch (i) {
    case 'l': return a1 / 4. * (std::pow(tau, 4.) - std::pow(tm, 4.)) + d1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.));
    case 'r': return a2 / 4. * (std::pow(tau, 4.) - std::pow(tp, 4.)) + d2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square_t", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  switch (i) {
    case 'l': return a1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) + d1 * (tau - tm);
    case 'r': return a2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) + d2 * (tau - tp);
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square_t", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_linear_t(char i) const -> double {
  switch (i) {
    case 'l': return a1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) + d1 * (tau - tm);
    case 'r': return a2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) + d2 * (tau - tp);
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linear_t", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  //    switch (i) {
  //        case 'l':
  //            if (iszero(tm)) {
  //                return a1 / 2. * std::pow(tau, 2.) + d1 * std::log(tau);
  //            } else {
  //                return a1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) + d1 * (std::log(tau) - std::log(tm));
  //            }
  //        case 'r':
  //            return a2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) + d2 * (std::log(tau) - std::log(tp));
  //        default:
  //            printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linear_t", "char i (which is %c) wrong", i);
  //            return ZEROVALUE;
  //    }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_const_t(char i) const -> double {
  switch (i) {
    case 'l':
      if (iszero(tm)) {
        return a1 / 2. * std::pow(tau, 2.) + d1 * std::log(tau);
      } else {
        return a1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) + d1 * (std::log(tau) - std::log(tm));
      }
    case 'r': return a2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) + d2 * (std::log(tau) - std::log(tp));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const_t", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  //    switch (i) {
  //        case 'l':
  //            if (iszero(tm)) {
  //                return a1 * tau - d1 / tau;
  //            } else {
  //                return a1 * (tau - tm) - d1 * (UNITVALUE / tau - UNITVALUE / tm);
  //            }
  //        case 'r':
  //            return a2 * (tau - tp) - d2 * (UNITVALUE / tau - UNITVALUE / tp);
  //        default:
  //            printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const_t", "char i (which is %c) wrong", i);
  //            return ZEROVALUE;
  //    }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_square_tau(char i) const -> double {
  switch (i) {
    case 'l': return c1 / 4. * (std::pow(tau, 4.) - std::pow(tm, 4.)) - e1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.));
    case 'r': return c2 / 4. * (std::pow(tau, 4.) - std::pow(tp, 4.)) - e2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square_tau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_linear_tau(char i) const -> double {
  switch (i) {
    case 'l': return c1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) - e1 * (tau - tm);
    case 'r': return c2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) - e2 * (tau - tp);
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linea_tau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_const_tau(char i) const -> double {
  switch (i) {
    case 'l':
      if (iszero(tm)) {
        return c1 / 2. * std::pow(tau, 2.) - e1 * std::log(tau);
      } else {
        return c1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) - e1 * (std::log(tau) - std::log(tm));
      }
    case 'r': return c2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) - e2 * (std::log(tau) - std::log(tp));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const_tau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_square_ttau(char i) const -> double {
  switch (i) {
    case 'l': return c1 / 4. * (std::pow(tau, 4.) - std::pow(tm, 4.)) + e1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.));
    case 'r': return c2 / 4. * (std::pow(tau, 4.) - std::pow(tp, 4.)) + e2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square_ttau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  //    switch (i) {
  //        case 'l':
  //            return c1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) + e1 * (tau - tm);
  //        case 'r':
  //            return c2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) + e2 * (tau - tp);
  //        default:
  //            printError("AGM::GreenfunctionAxisymmetricStokes::integrate_square_ttau", "char i (which is %c) wrong", i);
  //            return ZEROVALUE;
  //    }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_linear_ttau(char i) const -> double {
  switch (i) {
    case 'l': return c1 / 3. * (std::pow(tau, 3.) - std::pow(tm, 3.)) + e1 * (tau - tm);
    case 'r': return c2 / 3. * (std::pow(tau, 3.) - std::pow(tp, 3.)) + e2 * (tau - tp);
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linear_ttau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  //    switch (i) {
  //        case 'l':
  //            if (iszero(tm)) {
  //                return c1 / 2. * std::pow(tau, 2.) + e1 * std::log(tau);
  //            } else {
  //                return c1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) + e1 * (std::log(tau) - std::log(tm));
  //            }
  //        case 'r':
  //            return c2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) + e2 * (std::log(tau) - std::log(tp));
  //        default:
  //            printError("AGM::GreenfunctionAxisymmetricStokes::integrate_linear_ttau", "char i (which is %c) wrong", i);
  //            return ZEROVALUE;
  //    }
}

auto AGM::GreenfunctionAxisymmetricStokes::integrate_const_ttau(char i) const -> double {
  switch (i) {
    case 'l':
      if (iszero(tm)) {
        return c1 / 2. * std::pow(tau, 2.) + e1 * std::log(tau);
      } else {
        return c1 / 2. * (std::pow(tau, 2.) - std::pow(tm, 2.)) + e1 * (std::log(tau) - std::log(tm));
      }
    case 'r': return c2 / 2. * (std::pow(tau, 2.) - std::pow(tp, 2.)) + e2 * (std::log(tau) - std::log(tp));
    default:
      printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const_ttau", "char i (which is %c) wrong", i);
      return ZEROVALUE;
  }

  //    switch (i) {
  //        case 'l':
  //            if (iszero(tm)) {
  //                return c1 * tau - e1 / tau;
  //            } else {
  //                return c1 * (tau - tm) - e1 * (UNITVALUE / tau - UNITVALUE / tm);
  //            }
  //        case 'r':
  //            return c2 * (tau - tp) - e2 * (UNITVALUE / tau - UNITVALUE / tp);
  //        default:
  //            printError("AGM::GreenfunctionAxisymmetricStokes::integrate_const_ttau", "char i (which is %c) wrong", i);
  //            return ZEROVALUE;
  //    }
}
