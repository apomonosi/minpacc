/*
 * example_lmder1 - fit exponential model to 15 data points.
 *
 * Minimize sum of squares of  f_i(x) = y_i - (x[0] + t_i/(x[1]*s_i + x[2]*t_i))
 * where y_i are observations, t_i = i+1, s_i = 16-i-1.
 * User supplies analytical Jacobian via lmder1.
 */

#include <stdio.h>
#include <math.h>
#include "minpack.h"

static const double y[15] = {
    1.4e-1, 1.8e-1, 2.2e-1, 2.5e-1, 2.9e-1, 3.2e-1, 3.5e-1, 3.9e-1,
    3.7e-1, 5.8e-1, 7.3e-1, 9.6e-1, 1.34e0, 2.1e0,  4.39e0
};

/* Column-major Jacobian: element (row i, col j) = fjac[j*ldfjac + i] */
static void fcn(int m, int n, const double *x, double *fvec, double *fjac,
                int ldfjac, int *iflag, void *udata)
{
    (void)n; (void)udata;
    if (*iflag == 1) {
        for (int i = 0; i < m; i++) {
            double t1 = i + 1;
            double t2 = 16 - i - 1;
            double t3 = i >= 8 ? t2 : t1;
            fvec[i] = y[i] - (x[0] + t1/(x[1]*t2 + x[2]*t3));
        }
    } else if (*iflag == 2) {
        for (int i = 0; i < m; i++) {
            double t1 = i + 1;
            double t2 = 16 - i - 1;
            double t3 = i >= 8 ? t2 : t1;
            double t4 = (x[1]*t2 + x[2]*t3) * (x[1]*t2 + x[2]*t3);
            fjac[i]            = -1.0;
            fjac[i + ldfjac]   =  t1*t2/t4;
            fjac[i + 2*ldfjac] =  t1*t3/t4;
        }
    }
}

int main(void)
{
    const int m = 15, n = 3;
    double x[3] = {1.0, 1.0, 1.0};
    double xp[3], fvec[15], fvecp[15], err[15];
    double fjac[15*3];
    int ipvt[3];
    double wa[5*3 + 15];
    int info = 0;
    double tol = sqrt(minpack_dpmpar(1));

    /* Check derivatives */
    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 1, err);
    info = 1; fcn(m, n, x, fvec, fjac, m, &info, NULL);
    info = 2; fcn(m, n, x, fvec, fjac, m, &info, NULL);
    info = 1; fcn(m, n, xp, fvecp, fjac, m, &info, NULL);
    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 2, err);
    printf("Derivative check (1.0=correct, 0.0=incorrect):\n");
    for (int i = 0; i < m; i++) printf("  err[%2d] = %.4f\n", i, err[i]);

    /* Solve */
    x[0] = x[1] = x[2] = 1.0;
    minpack_lmder1(fcn, m, n, x, fvec, fjac, m, tol, &info, ipvt, wa, 30, NULL);

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
