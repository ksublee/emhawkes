#ifndef EMHAWKES_EXP_KERNEL_H
#define EMHAWKES_EXP_KERNEL_H

#include <cmath>

namespace emhawkes {
inline double integrated_component(double component, double beta, double dt,
                                   double scaled_time) {
    // Match the R small-decay limit and avoid zero times an infinite weight.
    if (component == 0.0) return 0.0;
    const double weight = scaled_time == 0.0 ? dt : -std::expm1(-scaled_time) / beta;
    return component * weight;
}
}

#endif
