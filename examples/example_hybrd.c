/*
 * example_hybrd - solve tridiagonal system using full hybrd interface.
 *
 * Find x[0..8] satisfying the same tridiagonal equations as example_hybrd1,
 * using the expert interface with banded Jacobian (ml=mu=1).
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
    const int n     = 9;
    const int ldfjac = n;
    const int lr    = (n*(n+1))/2;
    double x[9]     = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double diag[9]  = {1,1,1,1,1,1,1,1,1};
    double fvec[9], fjac[9*9], r[45], qtf[9];
    double wa1[9], wa2[9], wa3[9], wa4[9];
    int info = 0, nfev = 0;

    double xtol   = sqrt(minpack_dpmpar(1));
    int maxfev    = 2000;
    int ml = 1, mu = 1;
    double epsfcn = 0.0;
    int mode      = 2;
    double factor = 100.0;
    int nprint    = 0;

    minpack_hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn,
                  diag, mode, factor, nprint,
                  &info, &nfev, fjac, ldfjac, r, lr,
                  qtf, wa1, wa2, wa3, wa4, NULL);

    double fnorm = 0.0;
    for (int i = 0; i < n; i++) fnorm += fvec[i]*fvec[i];
    fnorm = sqrt(fnorm);

    printf("Final L2 norm of residuals: %15.7e\n", fnorm);
    printf("Number of function evaluations: %d\n", nfev);
    printf("Exit parameter:             %d\n", info);
    printf("Final approximate solution:\n");
    for (int i = 0; i < n; i++) printf("  x[%d] = %15.7e\n", i, x[i]);
    return 0;
}
