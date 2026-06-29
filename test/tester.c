#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <assert.h>

#include "testsuite.h"
#include "minpack.h"

static double
enorm(const int n, const double* x)
{
    double norm = 0.0;
    for (int i = 0; i < n; i++)
        norm += x[i]*x[i];
    return sqrt(norm);
}

static void
trial_hybrd_fcn(int n, const double* x, double* fvec, int* iflag, void* udata)
{
    (void)udata;
    if (*iflag == 0) return;
    for (int k = 0; k < n; k++) {
        double temp  = (3.0 - 2.0*x[k])*x[k];
        double temp1 = k != 0   ? x[k - 1] : 0.0;
        double temp2 = k != n-1 ? x[k + 1] : 0.0;
        fvec[k] = temp - temp1 - 2.0*temp2 + 1.0;
    }
}

/* Column-major Jacobian: element (row i, col j) = fjac[j*ldfjac + i] */
static void
trial_hybrj_fcn(int n, const double* x, double* fvec, double* fjac, int ldfjac,
                int* iflag, void* udata)
{
    if (*iflag == 1) {
        trial_hybrd_fcn(n, x, fvec, iflag, udata);
    } else {
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                fjac[j*ldfjac + i] = 0.0;
        for (int k = 0; k < n; k++) {
            fjac[k*ldfjac + k] = 3.0 - 4.0*x[k];      /* (row=k, col=k) */
            if (k != 0)   fjac[(k-1)*ldfjac + k] = -1.0; /* (row=k, col=k-1) */
            if (k != n-1) fjac[(k+1)*ldfjac + k] = -2.0; /* (row=k, col=k+1) */
        }
    }
}

static void
trial_lmder_fcn(int m, int n, const double* x, double* fvec, double* fjac,
                int ldfjac, int* iflag, void* data)
{
    assert(data);
    double* y = (double*)data;
    assert(m == 15);
    assert(n == 3);

    if (*iflag == 1) {
        for (int i = 0; i < m; i++) {
            double tmp1 = i + 1;
            double tmp2 = 16 - i - 1;
            double tmp3 = i >= 8 ? tmp2 : tmp1;
            fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
        }
    } else if (*iflag == 2) {
        /* column-major: col j has leading dimension ldfjac */
        assert(ldfjac == m);
        for (int i = 0; i < m; i++) {
            double tmp1 = i + 1;
            double tmp2 = 16 - i - 1;
            double tmp3 = i >= 8 ? tmp2 : tmp1;
            double tmp4 = (x[1]*tmp2 + x[2]*tmp3) * (x[1]*tmp2 + x[2]*tmp3);
            fjac[i]            = -1.0;
            fjac[i + ldfjac]   = tmp1*tmp2/tmp4;
            fjac[i + 2*ldfjac] = tmp1*tmp3/tmp4;
        }
    }
}

static void
trial_lmdif_fcn(int m, int n, const double* x, double* fvec, int* iflag, void* data)
{
    assert(data);
    double* y = (double*)data;
    assert(m == 15);
    assert(n == 3);
    if (*iflag == 0) return;
    for (int i = 0; i < m; i++) {
        double tmp1 = i + 1;
        double tmp2 = 16 - i - 1;
        double tmp3 = i >= 8 ? tmp2 : tmp1;
        fvec[i] = y[i] - (x[0] + tmp1/(x[1]*tmp2 + x[2]*tmp3));
    }
}

static void
trial_lmstr_fcn(int m, int n, const double* x, double* fvec, double* fjrow,
                int* iflag, void* udata)
{
    (void)udata; (void)m; (void)n;
    if (*iflag == 1) {
        fvec[0] = 10.0 * (x[1] - x[0]*x[0]);
        fvec[1] = 1.0 - x[0];
    } else if (*iflag == 2) {
        fjrow[0] = -20.0 * x[0];
        fjrow[1] = 10.0;
    } else {
        fjrow[0] = -1.0;
        fjrow[1] = 0.0;
    }
}

static int
test_hybrd1(void)
{
    int n = 9;
    int info = 0;
    int lwa = 180;
    double x[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double fvec[9], wa[180];
    double tol = sqrt(minpack_dpmpar(1));
    double ref[9] = {-0.5706545,-0.6816283,-0.7017325,
                     -0.7042129,-0.7013690,-0.6918656,
                     -0.6657920,-0.5960342,-0.4164121};

    minpack_hybrd1(trial_hybrd_fcn, n, x, fvec, tol, &info, wa, lwa, NULL);
    if (!check(info, 1, "hybrd1: info")) return 1;
    if (!check(enorm(n, fvec), 0.0, tol, "hybrd1: residual")) return 1;
    for (int i = 0; i < 9; i++)
        if (!check(x[i], ref[i], 10*tol, "hybrd1: solution")) return 1;
    return 0;
}

static int
test_hybrd(void)
{
    int n = 9;
    int info = 0, nfev = 0;
    double x[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double diag[9] = {1,1,1,1,1,1,1,1,1};
    double fvec[9], fjac[9*9];
    double r[45], qtf[9], wa1[9], wa2[9], wa3[9], wa4[9];
    double tol = sqrt(minpack_dpmpar(1));
    double ref[9] = {-0.5706545,-0.6816283,-0.7017325,
                     -0.7042129,-0.7013690,-0.6918656,
                     -0.6657920,-0.5960342,-0.4164121};

    minpack_hybrd(trial_hybrd_fcn, n, x, fvec, tol, 2000, 1, 1, 0.0, diag, 2, 100.0, 0,
                  &info, &nfev, fjac, 9, r, 45, qtf, wa1, wa2, wa3, wa4, NULL);
    if (!check(info, 1, "hybrd: info")) return 1;
    if (!check(nfev, 14, "hybrd: nfev")) return 1;
    if (!check(enorm(n, fvec), 0.0, tol, "hybrd: residual")) return 1;
    for (int i = 0; i < 9; i++)
        if (!check(x[i], ref[i], 10*tol, "hybrd: solution")) return 1;
    return 0;
}

static int
test_hybrj1(void)
{
    int n = 9;
    int info = 0, lwa = 180;
    double x[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double fvec[9], fjac[9*9], wa[180];
    double tol = sqrt(minpack_dpmpar(1));
    double ref[9] = {-0.5706545,-0.6816283,-0.7017325,
                     -0.7042129,-0.7013690,-0.6918656,
                     -0.6657920,-0.5960342,-0.4164121};

    minpack_hybrj1(trial_hybrj_fcn, n, x, fvec, fjac, n, tol, &info, wa, lwa, NULL);
    if (!check(info, 1, "hybrj1: info")) return 1;
    if (!check(enorm(n, fvec), 0.0, tol, "hybrj1: residual")) return 1;
    for (int i = 0; i < 9; i++)
        if (!check(x[i], ref[i], 10*tol, "hybrj1: solution")) return 1;
    return 0;
}

static int
test_hybrj(void)
{
    int n = 9;
    int info = 0, nfev = 0, njev = 0;
    double x[9] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
    double diag[9] = {1,1,1,1,1,1,1,1,1};
    double fvec[9], fjac[9*9];
    double r[45], qtf[9], wa1[9], wa2[9], wa3[9], wa4[9];
    double tol = sqrt(minpack_dpmpar(1));
    double ref[9] = {-0.5706545,-0.6816283,-0.7017325,
                     -0.7042129,-0.7013690,-0.6918656,
                     -0.6657920,-0.5960342,-0.4164121};

    minpack_hybrj(trial_hybrj_fcn, n, x, fvec, fjac, n, tol, 2000, diag, 2, 100.0, 0,
                  &info, &nfev, &njev, r, 45, qtf, wa1, wa2, wa3, wa4, NULL);
    if (!check(info, 1, "hybrj: info")) return 1;
    if (!check(nfev, 11, "hybrj: nfev")) return 1;
    if (!check(njev, 1, "hybrj: njev")) return 1;
    if (!check(enorm(n, fvec), 0.0, tol, "hybrj: residual")) return 1;
    for (int i = 0; i < 9; i++)
        if (!check(x[i], ref[i], 10*tol, "hybrj: solution")) return 1;
    return 0;
}

static int
test_lmder1(void)
{
    const double y[15] = {1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
                          3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0};
    const int m = 15, n = 3;
    int info = 0;
    double x[3] = {1,1,1}, xp[3];
    double fvec[15], fvecp[15], err[15];
    double fjac[15*3];
    int ipvt[3];
    double wa[5*3+15];
    double tol = sqrt(minpack_dpmpar(1));

    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 1, err);
    info = 1; trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, (void*)y);
    info = 2; trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, (void*)y);
    info = 1; trial_lmder_fcn(m, n, xp, fvecp, fjac, m, &info, (void*)y);
    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 2, err);
    for (int i = 0; i < 15; i++)
        if (!check(err[i], 1.0, tol, "lmder1: chkder")) return 1;

    minpack_lmder1(trial_lmder_fcn, m, n, x, fvec, fjac, m, tol, &info, ipvt, wa, 30, (void*)y);
    if (!check(info, 1, "lmder1: info")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "lmder1: x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "lmder1: x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "lmder1: x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "lmder1: residual")) return 1;
    return 0;
}

static int
test_lmder(void)
{
    const double y[15] = {1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
                          3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0};
    const int m = 15, n = 3;
    int info = 0, nfev = 0, njev = 0;
    double x[3] = {1,1,1}, xp[3];
    double fvec[15], fvecp[15], err[15];
    double fjac[15*3];
    int ipvt[3];
    double diag[3], qtf[3], wa1[3], wa2[3], wa3[3], wa4[15];
    double tol = sqrt(minpack_dpmpar(1));

    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 1, err);
    info = 1; trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, (void*)y);
    info = 2; trial_lmder_fcn(m, n, x, fvec, fjac, m, &info, (void*)y);
    info = 1; trial_lmder_fcn(m, n, xp, fvecp, fjac, m, &info, (void*)y);
    minpack_chkder(m, n, x, fvec, fjac, m, xp, fvecp, 2, err);
    for (int i = 0; i < 15; i++)
        if (!check(err[i], 1.0, tol, "lmder: chkder")) return 1;

    minpack_lmder(trial_lmder_fcn, m, n, x, fvec, fjac, m, tol, tol, 0.0, 2000, diag, 1,
                  100.0, 0, &info, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4, (void*)y);
    if (!check(info, 1, "lmder: info")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "lmder: x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "lmder: x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "lmder: x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "lmder: residual")) return 1;
    return 0;
}

static int
test_lmdif1(void)
{
    double y[15] = {1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
                    3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0};
    const int m = 15, n = 3;
    double x[3] = {1,1,1}, fvec[15];
    int info = 0;
    double tol = sqrt(minpack_dpmpar(1));
    int ipvt[3];
    int lwa = m*n + 5*n + m;
    double wa[lwa];

    minpack_lmdif1(trial_lmdif_fcn, 15, 3, x, fvec, tol, &info, ipvt, wa, lwa, (void*)y);
    if (!check(info, 1, "lmdif1: info")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "lmdif1: x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "lmdif1: x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "lmdif1: x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "lmdif1: residual")) return 1;
    return 0;
}

static int
test_lmdif(void)
{
    double y[15] = {1.4e-1,1.8e-1,2.2e-1,2.5e-1,2.9e-1,3.2e-1,3.5e-1,3.9e-1,
                    3.7e-1,5.8e-1,7.3e-1,9.6e-1,1.34e0,2.1e0,4.39e0};
    const int m = 15, n = 3;
    double x[3] = {1,1,1}, fvec[15];
    int info = 0, nfev = 0;
    double tol = sqrt(minpack_dpmpar(1));
    int ipvt[3];
    double fjac[m*n], diag[3], qtf[3], wa1[3], wa2[3], wa3[3], wa4[15];

    minpack_lmdif(trial_lmdif_fcn, 15, 3, x, fvec, tol, tol, 0.0, 2000, 0.0, diag, 1,
                  100.0, 0, &info, &nfev, fjac, 15, ipvt, qtf, wa1, wa2, wa3, wa4, (void*)y);
    if (!check(info, 1, "lmdif: info")) return 1;
    if (!check(x[0], 0.8241058e-1, 100*tol, "lmdif: x[0]")) return 1;
    if (!check(x[1], 0.1133037e+1, 100*tol, "lmdif: x[1]")) return 1;
    if (!check(x[2], 0.2343695e+1, 100*tol, "lmdif: x[2]")) return 1;
    if (!check(enorm(m, fvec), 0.9063596e-1, tol, "lmdif: residual")) return 1;
    return 0;
}

static int
test_lmstr1(void)
{
    const int m = 2, n = 2;
    double x[2] = {-1.2, 1.0}, fvec[2], fjac[4];
    int info = 0;
    double tol = sqrt(minpack_dpmpar(1));
    int ipvt[2];
    int lwa = m*n + 5*n + m;
    double wa[lwa];

    minpack_lmstr1(trial_lmstr_fcn, 2, 2, x, fvec, fjac, 2, tol, &info, ipvt, wa, lwa, NULL);
    if (!check(info, 4, "lmstr1: info")) return 1;
    if (!check(x[0], 1.0, 100*tol, "lmstr1: x[0]")) return 1;
    if (!check(x[1], 1.0, 100*tol, "lmstr1: x[1]")) return 1;
    if (!check(enorm(m, fvec), 0.0, tol, "lmstr1: residual")) return 1;
    return 0;
}

static int
test_lmstr(void)
{
    const int m = 2;
    double x[2] = {-1.2, 1.0}, fvec[2], fjac[4], diag[2];
    int info = 0, nfev = 0, njev = 0;
    double tol = sqrt(minpack_dpmpar(1));
    int ipvt[2];
    double qtf[2], wa1[2], wa2[2], wa3[2], wa4[2];

    minpack_lmstr(trial_lmstr_fcn, 2, 2, x, fvec, fjac, 2, tol, tol, 0.0, 2000, diag, 1,
                  100.0, 0, &info, &nfev, &njev, ipvt, qtf, wa1, wa2, wa3, wa4, NULL);
    if (!check(info, 4, "lmstr: info")) return 1;
    if (!check(nfev, 21, "lmstr: nfev")) return 1;
    if (!check(njev, 16, "lmstr: njev")) return 1;
    if (!check(x[0], 1.0, 100*tol, "lmstr: x[0]")) return 1;
    if (!check(x[1], 1.0, 100*tol, "lmstr: x[1]")) return 1;
    if (!check(enorm(m, fvec), 0.0, tol, "lmstr: residual")) return 1;
    return 0;
}

int
main(void)
{
    int stat = 0;
    stat += run("hybrd1", test_hybrd1);
    stat += run("hybrd ", test_hybrd);
    stat += run("hybrj1", test_hybrj1);
    stat += run("hybrj ", test_hybrj);
    stat += run("lmder1", test_lmder1);
    stat += run("lmder ", test_lmder);
    stat += run("lmdif1", test_lmdif1);
    stat += run("lmdif ", test_lmdif);
    stat += run("lmstr1", test_lmstr1);
    stat += run("lmstr ", test_lmstr);

    if (stat > 0)
        fprintf(stderr, "[FAIL] %d test(s) failed\n", stat);
    else
        fprintf(stderr, "[PASS] all tests passed\n");
    return stat == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
