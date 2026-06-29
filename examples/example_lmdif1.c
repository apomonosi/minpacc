/*
 * example_lmdif1 - fit exponential model using finite-difference Jacobian.
 *
 * Same problem as example_lmder1 but uses lmdif1 (no analytical Jacobian needed).
 */

#include <stdio.h>
#include <math.h>
#include "minpack.h"

static const double y[15] = {
    1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1, 3.9e-1,
    3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34e0, 2.1e0,  4.39e0
};

static void fcn(int m, int n, const double *x, double *fvec, int *iflag, void *udata)
{
    (void)n; (void)udata;
    if (*iflag == 0) return;
    for (int i = 0; i < m; i++) {
        double t1 = i + 1;
        double t2 = 16 - i - 1;
        double t3 = i >= 8 ? t2 : t1;
        fvec[i] = y[i] - (x[0] + t1/(x[1]*t2 + x[2]*t3));
    }
}

int main(void)
{
    const int m = 15, n = 3;
    double x[3] = {1.0, 1.0, 1.0};
    double fvec[15];
    int ipvt[3];
    int lwa = m*n + 5*n + m;
    double wa[lwa];
    int info = 0;
    double tol = sqrt(minpack_dpmpar(1));

    minpack_lmdif1(fcn, m, n, x, fvec, tol, &info, ipvt, wa, lwa, NULL);

    double fnorm = 0.0;
    for (int i = 0; i < m; i++) fnorm += fvec[i]*fvec[i];
    fnorm = sqrt(fnorm);

    printf("Final L2 norm of residuals: %15.7e\n", fnorm);
    printf("Exit parameter:             %d\n", info);
    printf("Final approximate solution:\n");
    printf("  x[0] = %15.7e  (expected ~0.08241058)\n", x[0]);
    printf("  x[1] = %15.7e  (expected ~1.133037)\n",   x[1]);
    printf("  x[2] = %15.7e  (expected ~2.343695)\n",   x[2]);
    return 0;
}
