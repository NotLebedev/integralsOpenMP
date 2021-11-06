#include <cmath>

#include "arcsin.h"

extern bool invalidate;

/**
 * Arc sine via series expansion
 * Integral of arcsin(x) is sqrt(1 - x^2) + x * arcsin(x)
 * Integral from 0 to 1 is (pi - 2) / 2 = 0.570796
 *
 * Defines:
 * SI_SERIES_PARALLEL -- define to enable parallel calculation of series
 * SI_SERIES_SIZE -- defines number of terms to use in series expansion, default 100
 * SI_PREPARE_COEFFICIENTS -- define to prepare array with precalculated coefficients
 **/
data_t arcsin(data_t x) {
#ifndef SI_SERIES_SIZE
#define SI_SERIES_SIZE 80
#endif

#ifdef SI_PREPARE_COEFFICIENTS
    static data_t coeffs[SI_SERIES_SIZE] = {0};
    static bool coeffs_ready = false;
#pragma omp threadprivate(coeffs)
#pragma omp threadprivate(coeffs_ready)
    if (invalidate)
        coeffs_ready = false;
    if (!coeffs_ready) {
        for (int n = 0; n < SI_SERIES_SIZE; n++) {
            coeffs[n] = tgamma(2.0 * n + 1) / (powl(4.0, n) * powl(tgamma(n + 1), 2.0) * (2.0 * n + 1));
        }
        coeffs_ready = true;
    }
#endif

    data_t result = 0.0;
#ifdef SI_SERIES_PARALLEL
#pragma omp parallel for shared(x), reduction(+:result), default(none)
#endif
    for (int n = 0; n < SI_SERIES_SIZE; n++) {
#ifdef SI_PREPARE_COEFFICIENTS
        result += powl(x, 2 * n + 1) * coeffs[n];
#else
        result += powl(x, 2 * n + 1) * tgamma(2.0 * n + 1) / (powl(4.0, n) * powl(tgamma(n + 1), 2.0) * (2.0 * n + 1));
#endif
    }

    return result;
}