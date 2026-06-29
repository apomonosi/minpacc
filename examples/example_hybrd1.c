/*
 * example_hybrd1 - solve tridiagonal system of 9 nonlinear equations.
 *
 * Find x[0..8] satisfying:
 *   (3 - 2*x[0])*x[0]              - 2*x[1]              = -1
 *   -x[k-1] + (3 - 2*x[k])*x[k]   - 2*x[k+1]            = -1  (k=1..7)
 *             -x[7]  + (3 - 2*x[8])*x[8]                 = -1
 */

#include <stdio.h>
#include <math.h>
#include "minpack.h"

static void fcn(int n, const double *x, double *fvec, int *iflag, void *udata)
{
    (void)udata;
    if (*iflag == 0) return;
    for (int k = 0; k < n; k++) {
        double temp  = (3.0 - 2.0*x[k])*x[k];
        double temp1 = k > 0   ? x[k-1] : 0.0;
        double temp2 = k < n-1 ? x[k+1] : 0.0;
        fvec[k] = temp - temp1 - 2.0*temp2 + 1.0;
    }
}

int main(void)
{
    const int n   = 9;
    const int lwa = (n*(3*n+13))/2;
    double x[9]   = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double fvec[9], wa[lwa];
    int info = 0;

    double tol = sqrt(minpack_dpmpar(1));

    minpack_hybrd1(fcn, n, x, fvec, tol, &info, wa, lwa, NULL);

    double fnorm = 0.0;
    for (int i = 0; i < n; i++) fnorm += fvec[i]*fvec[i];
    fnorm = sqrt(fnorm);

    printf("Final L2 norm of residuals: %15.7e\n", fnorm);
    printf("Exit parameter:             %d\n", info);
    printf("Final approximate solution:\n");
    for (int i = 0; i < n; i++) printf("  x[%d] = %15.7e\n", i, x[i]);

    /*
     * Expected output:
     *   Final L2 norm of residuals:  1.1926364e-08
     *   Exit parameter:             1
     *   x = -0.5706545  -0.6816283  -0.7017325
     *       -0.7042129  -0.7013690  -0.6918656
     *       -0.6657920  -0.5960342  -0.4164121
     */
    return 0;
}
