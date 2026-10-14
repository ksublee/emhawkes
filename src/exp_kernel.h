#ifndef EMHAWKES_EXP_KERNEL_H
#define EMHAWKES_EXP_KERNEL_H

#include <cmath>

namespace emhawkes {
inline long double integrated_component_beta(double component, double beta, double dt,
                                             double decay, double integral) {
    if (component == 0.0) return 0.0L;
    const long double time = dt, rate = beta, x = rate * time;
    if (x < 0.01L) {
        // Avoid cancellation near zero, including beta * dt underflow.
        return (component * time) * (time * (-0.5L + x * (1.0L / 3.0L + x * (-1.0L / 8.0L +
            x * (1.0L / 30.0L + x * (-1.0L / 144.0L + x / 840.0L))))));
    }
    // Multiply the component before dividing twice by beta. This also works
    // where long double has the same exponent range as double.
    return (time * (component * static_cast<long double>(decay)) - integral) / rate;
}

inline double integrated_component(double component, double beta, double dt,
                                   double scaled_time) {
    // Match the R small-decay limit and avoid zero times an infinite weight.
    if (component == 0.0) return 0.0;
    const double weight = scaled_time == 0.0 ? dt : -std::expm1(-scaled_time) / beta;
    return component * weight;
}
}

#endif
