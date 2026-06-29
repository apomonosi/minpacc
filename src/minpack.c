/* minpack.c - C translation of Minpack nonlinear solver library
 *
 * Original Minpack: Argonne National Laboratory, March 1980.
 *   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.
 * Fortran modernization: Jacob Williams, Sept 2021.
 * C translation: direct port (not f2c), 2024.
 *
 * All 2-D arrays use column-major (Fortran) storage order.
 * Element (row i, col j), 0-based, with leading dim ld: a[j*ld + i]
 */

#include "minpack.h"

#include <math.h>
#include <float.h>

/* -------------------------------------------------------------------------
 * Column-major 2D access macro (1-based indices, matching Fortran source).
 * Define locally in each function to bind the specific array + leading dim.
 * Example:  #define A(i,j)  a[((j)-1)*(lda) + ((i)-1)]
 * ------------------------------------------------------------------------- */

static const double epsmch = DBL_EPSILON;

/* =========================================================================
 * MACHINE CONSTANTS
 * ========================================================================= */

double minpack_dpmpar(int i)
{
    switch (i) {
    case 1: return DBL_EPSILON;
    case 2: return DBL_MIN;
    case 3: return DBL_MAX;
    default: return 0.0;
    }
}

/* =========================================================================
 * INTERNAL UTILITIES (static — not exposed in header)
 * ========================================================================= */

/*
 * enorm - Euclidean norm with three-level accumulation to prevent
 * overflow/underflow.  x is 0-based.
 */
static double enorm(int n, const double *x)
{
    static const double rdwarf = 3.834e-20;
    static const double rgiant = 1.304e19;

    double s1 = 0.0, s2 = 0.0, s3 = 0.0;
    double x1max = 0.0, x3max = 0.0;
    double agiant = rgiant / (double)n;

    for (int i = 0; i < n; i++) {
        double xabs = fabs(x[i]);
        if (xabs > rdwarf && xabs < agiant) {
            s2 += xabs * xabs;
        } else if (xabs <= rdwarf) {
            if (xabs <= x3max) {
                if (xabs != 0.0)
                    s3 += (xabs / x3max) * (xabs / x3max);
            } else {
                s3 = 1.0 + s3 * (x3max / xabs) * (x3max / xabs);
                x3max = xabs;
            }
        } else if (xabs <= x1max) {
            s1 += (xabs / x1max) * (xabs / x1max);
        } else {
            s1 = 1.0 + s1 * (x1max / xabs) * (x1max / xabs);
            x1max = xabs;
        }
    }

    if (s1 != 0.0) {
        return x1max * sqrt(s1 + (s2 / x1max) / x1max);
    } else if (s2 == 0.0) {
        return x3max * sqrt(s3);
    } else {
        if (s2 >= x3max)
            return sqrt(s2 * (1.0 + (x3max / s2) * (x3max * s3)));
        return sqrt(x3max * ((s2 / x3max) + (x3max * s3)));
    }
}

/*
 * dogleg - compute convex combination of Gauss-Newton and scaled gradient
 * directions that minimises (a*x - b) subject to ||d*x|| <= delta.
 *
 * r(lr)     : upper triangular matrix R stored row-wise (compact form)
 * diag(n)   : diagonal scaling matrix D
 * qtb(n)    : first n components of Q^T * b
 * wa1, wa2  : work arrays length n
 * x(n)      : output step direction
 */
static void dogleg(int n, const double *r, int lr,
                   const double *diag, const double *qtb, double delta,
                   double *x, double *wa1, double *wa2)
{
    (void)lr; /* used only for interface consistency */

    /* compute the Gauss-Newton direction by back-substitution */
    int jj = (n * (n + 1)) / 2 + 1; /* matches Fortran: jj = n*(n+1)/2 + 1 */
    for (int k = 1; k <= n; k++) {
        int j   = n - k + 1;  /* 1-based */
        int jp1 = j + 1;
        jj -= k;               /* 1-based diagonal index of column j */
        int l = jj + 1;        /* 1-based index of element above diagonal */
        double sum = 0.0;
        if (n >= jp1) {
            for (int i = jp1; i <= n; i++) {
                sum += r[l - 1] * x[i - 1];
                l++;
            }
        }
        double temp = r[jj - 1];          /* r(jj) diagonal */
        if (temp == 0.0) {
            l = j;
            for (int i = 1; i <= j; i++) {
                if (fabs(r[l - 1]) > fabs(temp)) temp = fabs(r[l - 1]);
                l += n - i;
            }
            temp = epsmch * temp;
            if (temp == 0.0) temp = epsmch;
        }
        x[j - 1] = (qtb[j - 1] - sum) / temp;
    }

    /* test whether the Gauss-Newton direction is acceptable */
    for (int j = 0; j < n; j++) {
        wa1[j] = 0.0;
        wa2[j] = diag[j] * x[j];
    }
    double qnorm = enorm(n, wa2);

    if (qnorm <= delta) return; /* Gauss-Newton step is within trust region */

    /* compute the scaled gradient direction */
    int l = 0;
    for (int j = 1; j <= n; j++) {
        double temp = qtb[j - 1];
        for (int i = j; i <= n; i++) {
            wa1[i - 1] += r[l] * temp;
            l++;
        }
        wa1[j - 1] /= diag[j - 1];
    }

    double gnorm = enorm(n, wa1);
    double sgnorm = 0.0;
    double alpha  = delta / qnorm;

    if (gnorm != 0.0) {
        for (int j = 0; j < n; j++)
            wa1[j] = (wa1[j] / gnorm) / diag[j];

        l = 0;
        for (int j = 1; j <= n; j++) {
            double sum = 0.0;
            for (int i = j; i <= n; i++) {
                sum += r[l] * wa1[i - 1];
                l++;
            }
            wa2[j - 1] = sum;
        }
        double temp = enorm(n, wa2);
        sgnorm = (gnorm / temp) / temp;
        alpha  = 0.0;

        if (sgnorm < delta) {
            double bnorm = enorm(n, qtb);
            temp = (bnorm / gnorm) * (bnorm / qnorm) * (sgnorm / delta);
            temp = temp - (delta / qnorm) * (sgnorm / delta) * (sgnorm / delta)
                 + sqrt((temp - (delta / qnorm)) * (temp - (delta / qnorm))
                        + (1.0 - (delta / qnorm) * (delta / qnorm))
                          * (1.0 - (sgnorm / delta) * (sgnorm / delta)));
            alpha = ((delta / qnorm) * (1.0 - (sgnorm / delta) * (sgnorm / delta))) / temp;
        }
    }

    double temp2 = (1.0 - alpha) * (sgnorm < delta ? sgnorm : delta);
    for (int j = 0; j < n; j++)
        x[j] = temp2 * wa1[j] + alpha * x[j];
}

/*
 * fdjac1 - forward-difference Jacobian approximation for n×n systems.
 * Supports banded Jacobians via ml (subdiags) and mu (superdiags).
 * fcn callback uses Fortran-style interface adapted for C (no udata here;
 * called via wrapper that carries udata).
 */
typedef void (*fcn_n_t)(int n, const double *x, double *fvec,
                        int *iflag, void *udata);

static void fdjac1(fcn_n_t fcn, void *udata,
                   int n, double *x, const double *fvec,
                   double *fjac, int ldfjac, int *iflag,
                   int ml, int mu, double epsfcn,
                   double *wa1, double *wa2)
{
#define FJ1(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    double eps  = sqrt(fabs(epsfcn) > epsmch ? fabs(epsfcn) : epsmch);
    int msum = ml + mu + 1;

    if (msum < n) {
        /* banded Jacobian */
        for (int k = 1; k <= msum; k++) {
            for (int j = k; j <= n; j += msum) {
                wa2[j - 1] = x[j - 1];
                double h = eps * fabs(wa2[j - 1]);
                if (h == 0.0) h = eps;
                x[j - 1] = wa2[j - 1] + h;
            }
            fcn(n, x, wa1, iflag, udata);
            if (*iflag < 0) { *iflag = -1; goto done1; }
            for (int j = k; j <= n; j += msum) {
                x[j - 1] = wa2[j - 1];
                double h = eps * fabs(wa2[j - 1]);
                if (h == 0.0) h = eps;
                for (int i = 1; i <= n; i++) {
                    FJ1(i, j) = 0.0;
                    if (i >= j - mu && i <= j + ml)
                        FJ1(i, j) = (wa1[i - 1] - fvec[i - 1]) / h;
                }
            }
        }
    } else {
        /* dense Jacobian */
        for (int j = 1; j <= n; j++) {
            double temp = x[j - 1];
            double h    = eps * fabs(temp);
            if (h == 0.0) h = eps;
            x[j - 1] = temp + h;
            fcn(n, x, wa1, iflag, udata);
            if (*iflag < 0) { *iflag = -1; goto done1; }
            x[j - 1] = temp;
            for (int i = 1; i <= n; i++)
                FJ1(i, j) = (wa1[i - 1] - fvec[i - 1]) / h;
        }
    }
done1:;
#undef FJ1
}

/*
 * fdjac2 - forward-difference Jacobian approximation for m×n systems.
 */
typedef void (*fcn_mn_t)(int m, int n, const double *x, double *fvec,
                         int *iflag, void *udata);

static void fdjac2(fcn_mn_t fcn, void *udata,
                   int m, int n, double *x, const double *fvec,
                   double *fjac, int ldfjac, int *iflag,
                   double epsfcn, double *wa)
{
#define FJ2(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    double eps = sqrt(fabs(epsfcn) > epsmch ? fabs(epsfcn) : epsmch);
    for (int j = 1; j <= n; j++) {
        double temp = x[j - 1];
        double h    = eps * fabs(temp);
        if (h == 0.0) h = eps;
        x[j - 1] = temp + h;
        fcn(m, n, x, wa, iflag, udata);
        if (*iflag < 0) goto done2;
        x[j - 1] = temp;
        for (int i = 1; i <= m; i++)
            FJ2(i, j) = (wa[i - 1] - fvec[i - 1]) / h;
    }
done2:;
#undef FJ2
}

/*
 * qrfac - QR factorisation with optional column pivoting.
 * a[lda*n] column-major; on output: strict upper triangle = R,
 * lower triangle = Householder vectors for Q.
 */
static void qrfac(int m, int n, double *a, int lda,
                  int pivot, int *ipvt, int lipvt,
                  double *rdiag, double *acnorm, double *wa)
{
#define AA(i,j)  a[((j)-1)*(lda)+((i)-1)]
    (void)lipvt;
    static const double p05 = 5.0e-2;

    for (int j = 1; j <= n; j++) {
        acnorm[j - 1] = enorm(m, &a[(j - 1) * lda]);
        rdiag[j - 1]  = acnorm[j - 1];
        wa[j - 1]     = rdiag[j - 1];
        if (pivot) ipvt[j - 1] = j;
    }

    int minmn = m < n ? m : n;
    for (int j = 1; j <= minmn; j++) {
        if (pivot) {
            int kmax = j;
            for (int k = j; k <= n; k++)
                if (rdiag[k - 1] > rdiag[kmax - 1]) kmax = k;
            if (kmax != j) {
                for (int i = 1; i <= m; i++) {
                    double temp = AA(i, j);
                    AA(i, j)    = AA(i, kmax);
                    AA(i, kmax) = temp;
                }
                rdiag[kmax - 1] = rdiag[j - 1];
                wa[kmax - 1]    = wa[j - 1];
                int k           = ipvt[j - 1];
                ipvt[j - 1]     = ipvt[kmax - 1];
                ipvt[kmax - 1]  = k;
            }
        }

        /* Householder transformation for column j */
        double ajnorm = enorm(m - j + 1, &a[(j - 1) * lda + (j - 1)]);
        if (ajnorm != 0.0) {
            if (AA(j, j) < 0.0) ajnorm = -ajnorm;
            for (int i = j; i <= m; i++) AA(i, j) /= ajnorm;
            AA(j, j) += 1.0;

            int jp1 = j + 1;
            if (n >= jp1) {
                for (int k = jp1; k <= n; k++) {
                    double sum = 0.0;
                    for (int i = j; i <= m; i++)
                        sum += AA(i, j) * AA(i, k);
                    double temp = sum / AA(j, j);
                    for (int i = j; i <= m; i++)
                        AA(i, k) -= temp * AA(i, j);

                    if (pivot && rdiag[k - 1] != 0.0) {
                        temp = AA(j, k) / rdiag[k - 1];
                        rdiag[k - 1] *= sqrt(1.0 - temp * temp > 0.0
                                             ? 1.0 - temp * temp : 0.0);
                        if (p05 * (rdiag[k - 1] / wa[k - 1])
                                * (rdiag[k - 1] / wa[k - 1]) <= epsmch) {
                            rdiag[k - 1] = enorm(m - j, &a[(k - 1) * lda + j]);
                            wa[k - 1]    = rdiag[k - 1];
                        }
                    }
                }
            }
        }
        rdiag[j - 1] = -ajnorm;
    }
#undef AA
}

/*
 * qform - accumulate orthogonal Q from its factored Householder form.
 * q[ldq*m] square column-major, overwritten on output.
 */
static void qform(int m, int n, double *q, int ldq, double *wa)
{
#define QQ(i,j)  q[((j)-1)*(ldq)+((i)-1)]
    int minmn = m < n ? m : n;

    /* zero strict upper triangle in first min(m,n) columns */
    if (minmn >= 2) {
        for (int j = 2; j <= minmn; j++)
            for (int i = 1; i < j; i++)
                QQ(i, j) = 0.0;
    }

    /* initialise remaining columns to identity */
    int np1 = n + 1;
    if (m >= np1) {
        for (int j = np1; j <= m; j++) {
            for (int i = 1; i <= m; i++) QQ(i, j) = 0.0;
            QQ(j, j) = 1.0;
        }
    }

    /* accumulate Q from factored form */
    for (int l = 1; l <= minmn; l++) {
        int k = minmn - l + 1;
        for (int i = k; i <= m; i++) {
            wa[i - 1] = QQ(i, k);
            QQ(i, k)  = 0.0;
        }
        QQ(k, k) = 1.0;
        if (wa[k - 1] != 0.0) {
            for (int j = k; j <= m; j++) {
                double sum = 0.0;
                for (int i = k; i <= m; i++)
                    sum += QQ(i, j) * wa[i - 1];
                double temp = sum / wa[k - 1];
                for (int i = k; i <= m; i++)
                    QQ(i, j) -= temp * wa[i - 1];
            }
        }
    }
#undef QQ
}

/*
 * qrsolv - solve damped QR least-squares system.
 * r[ldr*n] column-major, modified in lower triangle to hold S.
 */
static void qrsolv(int n, double *r, int ldr, const int *ipvt,
                   const double *diag, const double *qtb,
                   double *x, double *sdiag, double *wa)
{
#define RR(i,j)  r[((j)-1)*(ldr)+((i)-1)]
    /* copy R and Q^T*b; save diagonal of R in x */
    for (int j = 1; j <= n; j++) {
        for (int i = j; i <= n; i++) RR(i, j) = RR(j, i);
        x[j - 1]  = RR(j, j);
        wa[j - 1] = qtb[j - 1];
    }

    /* eliminate diagonal matrix D using Givens rotations */
    for (int j = 1; j <= n; j++) {
        int l = ipvt[j - 1];
        if (diag[l - 1] != 0.0) {
            for (int k = j; k <= n; k++) sdiag[k - 1] = 0.0;
            sdiag[j - 1] = diag[l - 1];
            double qtbpj = 0.0;
            for (int k = j; k <= n; k++) {
                if (sdiag[k - 1] != 0.0) {
                    double cs, sn;
                    if (fabs(RR(k, k)) >= fabs(sdiag[k - 1])) {
                        double tan = sdiag[k - 1] / RR(k, k);
                        cs = 0.5 / sqrt(0.25 + 0.25 * tan * tan);
                        sn = cs * tan;
                    } else {
                        double cotan = RR(k, k) / sdiag[k - 1];
                        sn = 0.5 / sqrt(0.25 + 0.25 * cotan * cotan);
                        cs = sn * cotan;
                    }
                    RR(k, k) = cs * RR(k, k) + sn * sdiag[k - 1];
                    double temp = cs * wa[k - 1] + sn * qtbpj;
                    qtbpj      = -sn * wa[k - 1] + cs * qtbpj;
                    wa[k - 1]  = temp;
                    int kp1 = k + 1;
                    if (n >= kp1) {
                        for (int i = kp1; i <= n; i++) {
                            temp        = cs * RR(i, k) + sn * sdiag[i - 1];
                            sdiag[i-1]  = -sn * RR(i, k) + cs * sdiag[i - 1];
                            RR(i, k)    = temp;
                        }
                    }
                }
            }
        }
        sdiag[j - 1] = RR(j, j);
        RR(j, j)     = x[j - 1];
    }

    /* solve triangular system for z, least-squares if singular */
    int nsing = n;
    for (int j = 1; j <= n; j++) {
        if (sdiag[j - 1] == 0.0 && nsing == n) nsing = j - 1;
        if (nsing < n) wa[j - 1] = 0.0;
    }
    if (nsing >= 1) {
        for (int k = 1; k <= nsing; k++) {
            int j   = nsing - k + 1;
            double sum = 0.0;
            int jp1 = j + 1;
            if (nsing >= jp1) {
                for (int i = jp1; i <= nsing; i++)
                    sum += RR(i, j) * wa[i - 1];
            }
            wa[j - 1] = (wa[j - 1] - sum) / sdiag[j - 1];
        }
    }

    for (int j = 1; j <= n; j++)
        x[ipvt[j - 1] - 1] = wa[j - 1];
#undef RR
}

/*
 * lmpar - determine Levenberg-Marquardt parameter.
 * r[ldr*n] column-major, upper triangle = R from QR; output: lower triangle = S.
 */
static void lmpar(int n, double *r, int ldr, const int *ipvt,
                  const double *diag, const double *qtb,
                  double delta, double *par,
                  double *x, double *sdiag, double *wa1, double *wa2)
{
#define RP(i,j)  r[((j)-1)*(ldr)+((i)-1)]
    static const double p1   = 1.0e-1;
    static const double p001 = 1.0e-3;
    double dwarf = DBL_MIN;

    /* compute Gauss-Newton direction (handle rank deficiency) */
    int nsing = n;
    for (int j = 1; j <= n; j++) {
        wa1[j - 1] = qtb[j - 1];
        if (RP(j, j) == 0.0 && nsing == n) nsing = j - 1;
        if (nsing < n) wa1[j - 1] = 0.0;
    }
    if (nsing >= 1) {
        for (int k = 1; k <= nsing; k++) {
            int j  = nsing - k + 1;
            wa1[j - 1] /= RP(j, j);
            double temp  = wa1[j - 1];
            int jm1 = j - 1;
            if (jm1 >= 1) {
                for (int i = 1; i <= jm1; i++)
                    wa1[i - 1] -= RP(i, j) * temp;
            }
        }
    }
    for (int j = 1; j <= n; j++)
        x[ipvt[j - 1] - 1] = wa1[j - 1];

    /* evaluate function at origin; test Gauss-Newton direction */
    int iter = 0;
    for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
    double dxnorm = enorm(n, wa2);
    double fp     = dxnorm - delta;

    if (fp <= p1 * delta) {
        if (iter == 0) *par = 0.0;
        goto done_lmpar;
    }

    double parl = 0.0;
    if (nsing >= n) {
        for (int j = 1; j <= n; j++) {
            int l    = ipvt[j - 1];
            wa1[j-1] = diag[l - 1] * (wa2[l - 1] / dxnorm);
        }
        for (int j = 1; j <= n; j++) {
            double sum = 0.0;
            int jm1 = j - 1;
            if (jm1 >= 1) {
                for (int i = 1; i <= jm1; i++)
                    sum += RP(i, j) * wa1[i - 1];
            }
            wa1[j - 1] = (wa1[j - 1] - sum) / RP(j, j);
        }
        double temp = enorm(n, wa1);
        parl = ((fp / delta) / temp) / temp;
    }

    /* upper bound */
    for (int j = 1; j <= n; j++) {
        double sum = 0.0;
        for (int i = 1; i <= j; i++)
            sum += RP(i, j) * qtb[i - 1];
        int l    = ipvt[j - 1];
        wa1[j-1] = sum / diag[l - 1];
    }
    double gnorm = enorm(n, wa1);
    double paru  = gnorm / delta;
    if (paru == 0.0) paru = dwarf / (delta < p1 ? delta : p1);

    *par = *par > parl ? *par : parl;
    *par = *par < paru ? *par : paru;
    if (*par == 0.0) *par = gnorm / dxnorm;

    /* iterate to find par */
    for (;;) {
        iter++;
        if (*par == 0.0) *par = dwarf > p001 * paru ? dwarf : p001 * paru;
        double temp = sqrt(*par);
        for (int j = 0; j < n; j++) wa1[j] = temp * diag[j];
        qrsolv(n, r, ldr, ipvt, wa1, qtb, x, sdiag, wa2);
        for (int j = 0; j < n; j++) wa2[j] = diag[j] * x[j];
        dxnorm = enorm(n, wa2);
        double fp_old = fp;
        fp = dxnorm - delta;

        if (fabs(fp) <= p1 * delta
            || (parl == 0.0 && fp <= fp_old && fp_old < 0.0)
            || iter == 10) {
            if (iter == 0) *par = 0.0;
            break;
        }

        /* Newton correction */
        for (int j = 1; j <= n; j++) {
            int l    = ipvt[j - 1];
            wa1[j-1] = diag[l - 1] * (wa2[l - 1] / dxnorm);
        }
        for (int j = 1; j <= n; j++) {
            wa1[j - 1] /= sdiag[j - 1];
            double temp  = wa1[j - 1];
            int jp1 = j + 1;
            if (n >= jp1) {
                for (int i = jp1; i <= n; i++)
                    wa1[i - 1] -= RP(i, j) * temp;
            }
        }
        double enorm_wa1 = enorm(n, wa1);
        double parc = ((fp / delta) / enorm_wa1) / enorm_wa1;

        if (fp > 0.0) parl = *par > parl ? *par : parl;
        if (fp < 0.0) paru = *par < paru ? *par : paru;
        *par = *par + parc;
        if (*par < parl) *par = parl;
    }
done_lmpar:;
#undef RP
}

/*
 * r1mpyq - postmultiply m×n matrix a by product of Givens rotations
 * encoded in v and w.
 */
static void r1mpyq(int m, int n, double *a, int lda,
                   const double *v, const double *w)
{
#define AM(i,j)  a[((j)-1)*(lda)+((i)-1)]
    int nm1 = n - 1;
    if (nm1 < 1) return;

    /* first set of rotations */
    for (int nmj = 1; nmj <= nm1; nmj++) {
        int j = n - nmj;
        double cs, sn;
        if (fabs(v[j - 1]) > 1.0) {
            cs = 1.0 / v[j - 1];
            sn = sqrt(1.0 - cs * cs);
        } else {
            sn = v[j - 1];
            cs = sqrt(1.0 - sn * sn);
        }
        for (int i = 1; i <= m; i++) {
            double temp = cs * AM(i, j) - sn * AM(i, n);
            AM(i, n)    = sn * AM(i, j) + cs * AM(i, n);
            AM(i, j)    = temp;
        }
    }

    /* second set of rotations */
    for (int j = 1; j <= nm1; j++) {
        double cs, sn;
        if (fabs(w[j - 1]) > 1.0) {
            cs = 1.0 / w[j - 1];
            sn = sqrt(1.0 - cs * cs);
        } else {
            sn = w[j - 1];
            cs = sqrt(1.0 - sn * sn);
        }
        for (int i = 1; i <= m; i++) {
            double temp = cs * AM(i, j) + sn * AM(i, n);
            AM(i, n)    = -sn * AM(i, j) + cs * AM(i, n);
            AM(i, j)    = temp;
        }
    }
#undef AM
}

/*
 * r1updt - rank-1 update of lower-trapezoidal matrix stored by columns.
 * s(ls) compact storage; u(m), v(n) input; w(m) output.
 * sing: set true if any diagonal of output S is zero.
 */
static void r1updt(int m, int n, double *s, int ls,
                   const double *u, double *v, double *w, int *sing)
{
    (void)ls;
    static const double p5  = 0.5;
    static const double p25 = 0.25;
    double giant = DBL_MAX;

    /* index of last diagonal element (1-based) */
    int jj = (n * (2 * m - n + 1)) / 2 - (m - n);

    /* copy last column of S into w */
    int l = jj;
    for (int i = n; i <= m; i++) {
        w[i - 1] = s[l - 1];
        l++;
    }

    int nm1 = n - 1;

    /* rotate v into a multiple of e_n, introducing spike into w */
    if (nm1 >= 1) {
        for (int nmj = 1; nmj <= nm1; nmj++) {
            int j = n - nmj;
            jj   -= (m - j + 1);
            w[j - 1] = 0.0;
            if (v[j - 1] != 0.0) {
                double cs, sn, tau;
                if (fabs(v[n - 1]) >= fabs(v[j - 1])) {
                    double tan = v[j - 1] / v[n - 1];
                    cs  = p5 / sqrt(p25 + p25 * tan * tan);
                    sn  = cs * tan;
                    tau = sn;
                } else {
                    double cotan = v[n - 1] / v[j - 1];
                    sn  = p5 / sqrt(p25 + p25 * cotan * cotan);
                    cs  = sn * cotan;
                    tau = 1.0;
                    if (fabs(cs) * giant > 1.0) tau = 1.0 / cs;
                }
                v[n - 1] = sn * v[j - 1] + cs * v[n - 1];
                v[j - 1] = tau;
                l = jj;
                for (int i = j; i <= m; i++) {
                    double temp = cs * s[l - 1] - sn * w[i - 1];
                    w[i - 1]   = sn * s[l - 1] + cs * w[i - 1];
                    s[l - 1]   = temp;
                    l++;
                }
            }
        }
    }

    /* add spike from rank-1 update to w */
    for (int i = 0; i < m; i++) w[i] += v[n - 1] * u[i];

    /* eliminate the spike */
    *sing = 0;
    if (nm1 >= 1) {
        for (int j = 1; j <= nm1; j++) {
            if (w[j - 1] != 0.0) {
                double cs, sn, tau;
                if (fabs(s[jj - 1]) >= fabs(w[j - 1])) {
                    double tan = w[j - 1] / s[jj - 1];
                    cs  = p5 / sqrt(p25 + p25 * tan * tan);
                    sn  = cs * tan;
                    tau = sn;
                } else {
                    double cotan = s[jj - 1] / w[j - 1];
                    sn  = p5 / sqrt(p25 + p25 * cotan * cotan);
                    cs  = sn * cotan;
                    tau = 1.0;
                    if (fabs(cs) * giant > 1.0) tau = 1.0 / cs;
                }
                l = jj;
                for (int i = j; i <= m; i++) {
                    double temp = cs * s[l - 1] + sn * w[i - 1];
                    w[i - 1]   = -sn * s[l - 1] + cs * w[i - 1];
                    s[l - 1]   = temp;
                    l++;
                }
                w[j - 1] = tau;
            }
            if (s[jj - 1] == 0.0) *sing = 1;
            jj += (m - j + 1);
        }
    }

    /* move w back into last column of S */
    l = jj;
    for (int i = n; i <= m; i++) {
        s[l - 1] = w[i - 1];
        l++;
    }
    if (s[jj - 1] == 0.0) *sing = 1;
}

/*
 * rwupdt - update n×n upper-triangular R and RHS b when a new row w is added.
 * Returns cos/sin arrays for the n Givens rotations applied.
 */
static void rwupdt(int n, double *r, int ldr,
                   const double *w, double *b, double *alpha,
                   double *gcos, double *gsin)
{
#define RW(i,j)  r[((j)-1)*(ldr)+((i)-1)]
    static const double p5  = 0.5;
    static const double p25 = 0.25;

    for (int j = 1; j <= n; j++) {
        double rowj = w[j - 1];
        int jm1 = j - 1;
        if (jm1 >= 1) {
            for (int i = 1; i <= jm1; i++) {
                double temp = gcos[i - 1] * RW(i, j) + gsin[i - 1] * rowj;
                rowj        = -gsin[i - 1] * RW(i, j) + gcos[i - 1] * rowj;
                RW(i, j)    = temp;
            }
        }
        gcos[j - 1] = 1.0;
        gsin[j - 1] = 0.0;
        if (rowj != 0.0) {
            if (fabs(RW(j, j)) >= fabs(rowj)) {
                double tan  = rowj / RW(j, j);
                gcos[j - 1] = p5 / sqrt(p25 + p25 * tan * tan);
                gsin[j - 1] = gcos[j - 1] * tan;
            } else {
                double cotan = RW(j, j) / rowj;
                gsin[j - 1]  = p5 / sqrt(p25 + p25 * cotan * cotan);
                gcos[j - 1]  = gsin[j - 1] * cotan;
            }
            RW(j, j)    = gcos[j - 1] * RW(j, j) + gsin[j - 1] * rowj;
            double temp = gcos[j - 1] * b[j - 1] + gsin[j - 1] * (*alpha);
            *alpha      = -gsin[j - 1] * b[j - 1] + gcos[j - 1] * (*alpha);
            b[j - 1]    = temp;
        }
    }
#undef RW
}

/* =========================================================================
 * PUBLIC SOLVERS
 * ========================================================================= */

/* -------------------------------------------------------------------------
 * chkder - verify user-supplied Jacobian
 * ------------------------------------------------------------------------- */
void minpack_chkder(int m, int n, const double *x, double *fvec,
                    double *fjac, int ldfjac, double *xp, double *fvecp,
                    int mode, double *err)
{
#define FJC(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    double eps    = sqrt(epsmch);
    double factor = 100.0;
    double epsf   = factor * epsmch;
    double epslog = log10(eps);

    if (mode == 2) {
        for (int i = 0; i < m; i++) err[i] = 0.0;
        for (int j = 1; j <= n; j++) {
            double temp = fabs(x[j - 1]);
            if (temp == 0.0) temp = 1.0;
            for (int i = 1; i <= m; i++)
                err[i - 1] += temp * FJC(i, j);
        }
        for (int i = 1; i <= m; i++) {
            double temp = 1.0;
            if (fvec[i - 1] != 0.0 && fvecp[i - 1] != 0.0
                && fabs(fvecp[i - 1] - fvec[i - 1]) >= epsf * fabs(fvec[i - 1])) {
                temp = eps * fabs((fvecp[i - 1] - fvec[i - 1]) / eps - err[i - 1])
                     / (fabs(fvec[i - 1]) + fabs(fvecp[i - 1]));
            }
            err[i - 1] = 1.0;
            if (temp > epsmch && temp < eps)
                err[i - 1] = (log10(temp) - epslog) / epslog;
            if (temp >= eps)
                err[i - 1] = 0.0;
        }
    } else { /* mode == 1 */
        for (int j = 0; j < n; j++) {
            double temp = eps * fabs(x[j]);
            if (temp == 0.0) temp = eps;
            xp[j] = x[j] + temp;
        }
    }
#undef FJC
}

/* -------------------------------------------------------------------------
 * hybrd - solve n×n nonlinear system, numerical Jacobian
 * ------------------------------------------------------------------------- */
void minpack_hybrd(minpack_func fcn, int n, double *x, double *fvec,
                   double xtol, int maxfev, int ml, int mu, double epsfcn,
                   double *diag, int mode, double factor, int nprint,
                   int *info, int *nfev, double *fjac, int ldfjac,
                   double *r, int lr, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata)
{
#define FH(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    static const double p1    = 1.0e-1;
    static const double p5    = 5.0e-1;
    static const double p001  = 1.0e-3;
    static const double p0001 = 1.0e-4;

    *info = 0;
    int iflag = 0;
    *nfev = 0;

    /* validate inputs */
    if (n <= 0 || xtol < 0.0 || maxfev <= 0 || ml < 0 || mu < 0
        || factor <= 0.0 || ldfjac < n || lr < (n * (n + 1)) / 2)
        goto done_hybrd;
    if (mode == 2) {
        for (int j = 0; j < n; j++)
            if (diag[j] <= 0.0) goto done_hybrd;
    }

    iflag = 1;
    fcn(n, x, fvec, &iflag, udata);
    *nfev = 1;
    if (iflag < 0) goto done_hybrd;

    double fnorm = enorm(n, fvec);
    int msum = ml + mu + 1 < n ? ml + mu + 1 : n;
    int iter = 1, ncsuc = 0, ncfail = 0, nslow1 = 0, nslow2 = 0;
    double delta = 0.0, xnorm = 0.0;

    /* outer loop */
    for (;;) {
        int jeval = 1;
        iflag = 2;
        fdjac1(fcn, udata, n, x, fvec, fjac, ldfjac, &iflag,
               ml, mu, epsfcn, wa1, wa2);
        *nfev += msum;
        if (iflag < 0) goto done_hybrd;

        int iwa[1];
        qrfac(n, n, fjac, ldfjac, 0, iwa, 1, wa1, wa2, wa3);

        if (iter == 1) {
            if (mode != 2) {
                for (int j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) diag[j] = 1.0;
                }
            }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if (delta == 0.0) delta = factor;
        }

        /* Q^T * fvec → qtf */
        for (int i = 0; i < n; i++) qtf[i] = fvec[i];
        for (int j = 1; j <= n; j++) {
            if (FH(j, j) != 0.0) {
                double sum = 0.0;
                for (int i = j; i <= n; i++) sum += FH(i, j) * qtf[i - 1];
                double temp = -sum / FH(j, j);
                for (int i = j; i <= n; i++) qtf[i - 1] += FH(i, j) * temp;
            }
        }

        /* copy upper triangle of QR into r (row-wise compact) */
        int sing = 0;
        for (int j = 1; j <= n; j++) {
            int l   = j;
            int jm1 = j - 1;
            if (jm1 >= 1) {
                for (int i = 1; i <= jm1; i++) {
                    r[l - 1] = FH(i, j);
                    l += n - i;
                }
            }
            r[l - 1] = wa1[j - 1];
            if (wa1[j - 1] == 0.0) sing = 1;
        }

        qform(n, n, fjac, ldfjac, wa1);

        if (mode != 2) {
            for (int j = 0; j < n; j++)
                if (wa2[j] > diag[j]) diag[j] = wa2[j];
        }

        /* inner loop */
        for (;;) {
            if (nprint > 0 && (iter - 1) % nprint == 0) {
                iflag = 0;
                fcn(n, x, fvec, &iflag, udata);
                if (iflag < 0) goto done_hybrd;
            }

            dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);

            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm(n, wa3);
            if (iter == 1) delta = delta < pnorm ? delta : pnorm;

            iflag = 1;
            fcn(n, wa2, wa4, &iflag, udata);
            (*nfev)++;
            if (iflag < 0) goto done_hybrd;

            double fnorm1 = enorm(n, wa4);
            double actred = -1.0;
            if (fnorm1 < fnorm) actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);

            int l = 1;
            for (int i = 1; i <= n; i++) {
                double sum = 0.0;
                for (int j = i; j <= n; j++) {
                    sum += r[l - 1] * wa1[j - 1];
                    l++;
                }
                wa3[i - 1] = qtf[i - 1] + sum;
            }
            double temp_e = enorm(n, wa3);
            double prered = 0.0;
            if (temp_e < fnorm) prered = 1.0 - (temp_e / fnorm) * (temp_e / fnorm);

            double ratio = 0.0;
            if (prered > 0.0) ratio = actred / prered;

            if (ratio >= p1) {
                ncfail = 0;
                ncsuc++;
                if (ratio >= p5 || ncsuc > 1)
                    delta = delta > pnorm / p5 ? delta : pnorm / p5;
                if (fabs(ratio - 1.0) <= p1) delta = pnorm / p5;
            } else {
                ncsuc = 0;
                ncfail++;
                delta *= p5;
            }

            if (ratio >= p0001) {
                for (int j = 0; j < n; j++) {
                    x[j]    = wa2[j];
                    wa2[j]  = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                xnorm  = enorm(n, wa2);
                fnorm  = fnorm1;
                iter++;
            }

            nslow1++;
            if (actred >= p001) nslow1 = 0;
            if (jeval) nslow2++;
            if (actred >= p1) nslow2 = 0;

            if (delta <= xtol * xnorm || fnorm == 0.0) *info = 1;
            if (*info != 0) goto done_hybrd;

            if (*nfev >= maxfev)                               *info = 2;
            if (p1 * (p1 * delta > pnorm ? p1 * delta : pnorm)
                <= epsmch * xnorm)                             *info = 3;
            if (nslow2 == 5)                                   *info = 4;
            if (nslow1 == 10)                                  *info = 5;
            if (*info != 0) goto done_hybrd;

            if (ncfail == 2) break; /* cycle outer: recompute Jacobian */

            for (int j = 1; j <= n; j++) {
                double sum = 0.0;
                for (int i = 1; i <= n; i++) sum += FH(i, j) * wa4[i - 1];
                wa2[j - 1] = (sum - wa3[j - 1]) / pnorm;
                wa1[j - 1] = diag[j - 1] * ((diag[j - 1] * wa1[j - 1]) / pnorm);
                if (ratio >= p0001) qtf[j - 1] = sum;
            }

            r1updt(n, n, r, lr, wa1, wa2, wa3, &sing);
            r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
            r1mpyq(1, n, qtf, 1, wa2, wa3);
            jeval = 0;
        } /* inner */
    } /* outer */

done_hybrd:
    if (iflag < 0) *info = iflag;
    iflag = 0;
    if (nprint > 0) fcn(n, x, fvec, &iflag, udata);
#undef FH
}

/* -------------------------------------------------------------------------
 * hybrd1 - simplified driver for hybrd
 * ------------------------------------------------------------------------- */
void minpack_hybrd1(minpack_func fcn, int n, double *x, double *fvec,
                    double tol, int *info, double *wa, int lwa, void *udata)
{
    *info = 0;
    if (n <= 0 || tol < 0.0 || lwa < (n * (3 * n + 13)) / 2) return;

    int maxfev = 200 * (n + 1);
    double xtol = tol;
    int ml = n - 1, mu = n - 1;
    double epsfcn = 0.0;
    for (int j = 0; j < n; j++) wa[j] = 1.0;
    int nprint = 0;
    int lr = (n * (n + 1)) / 2;
    int index = 6 * n + lr;
    int nfev;

    minpack_hybrd(fcn, n, x, fvec, xtol, maxfev, ml, mu, epsfcn,
                  wa, 2, 100.0, nprint, info, &nfev,
                  wa + index, n, wa + 6 * n, lr, wa + n,
                  wa + 2 * n, wa + 3 * n, wa + 4 * n, wa + 5 * n,
                  udata);
    if (*info == 5) *info = 4;
}

/* -------------------------------------------------------------------------
 * hybrj - solve n×n nonlinear system, user-supplied Jacobian
 * ------------------------------------------------------------------------- */
void minpack_hybrj(minpack_fcn_hybrj fcn, int n, double *x, double *fvec,
                   double *fjac, int ldfjac, double xtol, int maxfev,
                   double *diag, int mode, double factor, int nprint,
                   int *info, int *nfev, int *njev,
                   double *r, int lr, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata)
{
#define FRJ(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    static const double p1    = 1.0e-1;
    static const double p5    = 5.0e-1;
    static const double p001  = 1.0e-3;
    static const double p0001 = 1.0e-4;

    *info = 0;
    int iflag = 0;
    *nfev = 0;
    *njev = 0;

    if (n <= 0 || ldfjac < n || xtol < 0.0 || maxfev <= 0
        || factor <= 0.0 || lr < (n * (n + 1)) / 2)
        goto done_hybrj;
    if (mode == 2) {
        for (int j = 0; j < n; j++)
            if (diag[j] <= 0.0) goto done_hybrj;
    }

    iflag = 1;
    fcn(n, x, fvec, fjac, ldfjac, &iflag, udata);
    *nfev = 1;
    if (iflag < 0) goto done_hybrj;
    double fnorm = enorm(n, fvec);

    int iter = 1, ncsuc = 0, ncfail = 0, nslow1 = 0, nslow2 = 0;
    double delta = 0.0, xnorm = 0.0;

    for (;;) { /* outer */
        int jeval = 1;
        iflag = 2;
        fcn(n, x, fvec, fjac, ldfjac, &iflag, udata);
        (*njev)++;
        if (iflag < 0) goto done_hybrj;

        int iwa[1];
        qrfac(n, n, fjac, ldfjac, 0, iwa, 1, wa1, wa2, wa3);

        if (iter == 1) {
            if (mode != 2) {
                for (int j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) diag[j] = 1.0;
                }
            }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if (delta == 0.0) delta = factor;
        }

        for (int i = 0; i < n; i++) qtf[i] = fvec[i];
        for (int j = 1; j <= n; j++) {
            if (FRJ(j, j) != 0.0) {
                double sum = 0.0;
                for (int i = j; i <= n; i++) sum += FRJ(i, j) * qtf[i - 1];
                double temp = -sum / FRJ(j, j);
                for (int i = j; i <= n; i++) qtf[i - 1] += FRJ(i, j) * temp;
            }
        }

        int sing = 0;
        for (int j = 1; j <= n; j++) {
            int l   = j;
            int jm1 = j - 1;
            if (jm1 >= 1) {
                for (int i = 1; i <= jm1; i++) {
                    r[l - 1] = FRJ(i, j);
                    l += n - i;
                }
            }
            r[l - 1] = wa1[j - 1];
            if (wa1[j - 1] == 0.0) sing = 1;
        }

        qform(n, n, fjac, ldfjac, wa1);

        if (mode != 2) {
            for (int j = 0; j < n; j++)
                if (wa2[j] > diag[j]) diag[j] = wa2[j];
        }

        for (;;) { /* inner */
            if (nprint > 0 && (iter - 1) % nprint == 0) {
                iflag = 0;
                fcn(n, x, fvec, fjac, ldfjac, &iflag, udata);
                if (iflag < 0) goto done_hybrj;
            }

            dogleg(n, r, lr, diag, qtf, delta, wa1, wa2, wa3);
            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm(n, wa3);
            if (iter == 1) delta = delta < pnorm ? delta : pnorm;

            iflag = 1;
            fcn(n, wa2, wa4, fjac, ldfjac, &iflag, udata);
            (*nfev)++;
            if (iflag < 0) goto done_hybrj;

            double fnorm1 = enorm(n, wa4);
            double actred = -1.0;
            if (fnorm1 < fnorm) actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);

            int l = 1;
            for (int i = 1; i <= n; i++) {
                double sum = 0.0;
                for (int j = i; j <= n; j++) {
                    sum += r[l - 1] * wa1[j - 1];
                    l++;
                }
                wa3[i - 1] = qtf[i - 1] + sum;
            }
            double tmp_e = enorm(n, wa3);
            double prered = 0.0;
            if (tmp_e < fnorm) prered = 1.0 - (tmp_e / fnorm) * (tmp_e / fnorm);

            double ratio = 0.0;
            if (prered > 0.0) ratio = actred / prered;

            if (ratio >= p1) {
                ncfail = 0;
                ncsuc++;
                if (ratio >= p5 || ncsuc > 1)
                    delta = delta > pnorm / p5 ? delta : pnorm / p5;
                if (fabs(ratio - 1.0) <= p1) delta = pnorm / p5;
            } else {
                ncsuc = 0;
                ncfail++;
                delta *= p5;
            }

            if (ratio >= p0001) {
                for (int j = 0; j < n; j++) {
                    x[j]    = wa2[j];
                    wa2[j]  = diag[j] * x[j];
                    fvec[j] = wa4[j];
                }
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }

            nslow1++;
            if (actred >= p001) nslow1 = 0;
            if (jeval) nslow2++;
            if (actred >= p1) nslow2 = 0;

            if (delta <= xtol * xnorm || fnorm == 0.0) *info = 1;
            if (*info != 0) goto done_hybrj;

            if (*nfev >= maxfev)                                   *info = 2;
            if (p1 * (p1 * delta > pnorm ? p1 * delta : pnorm)
                <= epsmch * xnorm)                                 *info = 3;
            if (nslow2 == 5)                                       *info = 4;
            if (nslow1 == 10)                                      *info = 5;
            if (*info != 0) goto done_hybrj;

            if (ncfail == 2) break; /* cycle outer */

            for (int j = 1; j <= n; j++) {
                double sum = 0.0;
                for (int i = 1; i <= n; i++) sum += FRJ(i, j) * wa4[i - 1];
                wa2[j - 1] = (sum - wa3[j - 1]) / pnorm;
                wa1[j - 1] = diag[j - 1] * ((diag[j - 1] * wa1[j - 1]) / pnorm);
                if (ratio >= p0001) qtf[j - 1] = sum;
            }

            r1updt(n, n, r, lr, wa1, wa2, wa3, &sing);
            r1mpyq(n, n, fjac, ldfjac, wa2, wa3);
            r1mpyq(1, n, qtf, 1, wa2, wa3);
            jeval = 0;
        } /* inner */
    } /* outer */

done_hybrj:
    if (iflag < 0) *info = iflag;
    iflag = 0;
    if (nprint > 0) fcn(n, x, fvec, fjac, ldfjac, &iflag, udata);
#undef FRJ
}

/* -------------------------------------------------------------------------
 * hybrj1 - simplified driver for hybrj
 * ------------------------------------------------------------------------- */
void minpack_hybrj1(minpack_fcn_hybrj fcn, int n, double *x, double *fvec,
                    double *fjac, int ldfjac, double tol, int *info,
                    double *wa, int lwa, void *udata)
{
    *info = 0;
    if (n <= 0 || ldfjac < n || tol < 0.0 || lwa < (n * (n + 13)) / 2)
        return;

    int maxfev = 100 * (n + 1);
    double xtol = tol;
    int lr = (n * (n + 1)) / 2;
    for (int j = 0; j < n; j++) wa[j] = 1.0;
    int nfev, njev;

    minpack_hybrj(fcn, n, x, fvec, fjac, ldfjac, xtol, maxfev,
                  wa, 2, 100.0, 0, info, &nfev, &njev,
                  wa + 6 * n, lr, wa + n,
                  wa + 2 * n, wa + 3 * n, wa + 4 * n, wa + 5 * n,
                  udata);
    if (*info == 5) *info = 4;
}

/* -------------------------------------------------------------------------
 * lmder - Levenberg-Marquardt, user-supplied Jacobian
 * ------------------------------------------------------------------------- */
void minpack_lmder(minpack_fcn_lmder fcn, int m, int n, double *x,
                   double *fvec, double *fjac, int ldfjac,
                   double ftol, double xtol, double gtol,
                   int maxfev, double *diag, int mode, double factor,
                   int nprint, int *info, int *nfev, int *njev,
                   int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata)
{
#define FLD(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    static const double p1    = 1.0e-1;
    static const double p5    = 5.0e-1;
    static const double p25   = 2.5e-1;
    static const double p75   = 7.5e-1;
    static const double p0001 = 1.0e-4;

    *info = 0;
    int iflag = 0;
    *nfev = 0;
    *njev = 0;

    if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0
        || gtol < 0.0 || maxfev <= 0 || factor <= 0.0)
        goto done_lmder;
    if (mode == 2) {
        for (int j = 0; j < n; j++)
            if (diag[j] <= 0.0) goto done_lmder;
    }

    iflag = 1;
    fcn(m, n, x, fvec, fjac, ldfjac, &iflag, udata);
    *nfev = 1;
    if (iflag < 0) goto done_lmder;
    double fnorm = enorm(m, fvec);
    double par   = 0.0;
    int    iter  = 1;
    double delta = 0.0, xnorm = 0.0;

    for (;;) { /* outer */
        iflag = 2;
        fcn(m, n, x, fvec, fjac, ldfjac, &iflag, udata);
        (*njev)++;
        if (iflag < 0) goto done_lmder;

        if (nprint > 0 && (iter - 1) % nprint == 0) {
            iflag = 0;
            fcn(m, n, x, fvec, fjac, ldfjac, &iflag, udata);
            if (iflag < 0) goto done_lmder;
        }

        qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

        if (iter == 1) {
            if (mode != 2) {
                for (int j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) diag[j] = 1.0;
                }
            }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if (delta == 0.0) delta = factor;
        }

        /* Q^T * fvec, first n components in qtf */
        for (int i = 0; i < m; i++) wa4[i] = fvec[i];
        for (int j = 1; j <= n; j++) {
            if (FLD(j, j) != 0.0) {
                double sum = 0.0;
                for (int i = j; i <= m; i++) sum += FLD(i, j) * wa4[i - 1];
                double temp = -sum / FLD(j, j);
                for (int i = j; i <= m; i++) wa4[i - 1] += FLD(i, j) * temp;
            }
            FLD(j, j)    = wa1[j - 1];
            qtf[j - 1]   = wa4[j - 1];
        }

        /* scaled gradient norm */
        double gnorm = 0.0;
        if (fnorm != 0.0) {
            for (int j = 1; j <= n; j++) {
                int l = ipvt[j - 1];
                if (wa2[l - 1] != 0.0) {
                    double sum = 0.0;
                    for (int i = 1; i <= j; i++)
                        sum += FLD(i, j) * (qtf[i - 1] / fnorm);
                    double v = fabs(sum / wa2[l - 1]);
                    if (v > gnorm) gnorm = v;
                }
            }
        }
        if (gnorm <= gtol) { *info = 4; goto done_lmder; }

        if (mode != 2) {
            for (int j = 0; j < n; j++)
                if (wa2[j] > diag[j]) diag[j] = wa2[j];
        }

        for (;;) { /* inner */
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par,
                  wa1, wa2, wa3, wa4);

            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm(n, wa3);
            if (iter == 1) delta = delta < pnorm ? delta : pnorm;

            iflag = 1;
            fcn(m, n, wa2, wa4, fjac, ldfjac, &iflag, udata);
            (*nfev)++;
            if (iflag < 0) goto done_lmder;
            double fnorm1 = enorm(m, wa4);

            double actred = -1.0;
            if (p1 * fnorm1 < fnorm)
                actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);

            for (int j = 1; j <= n; j++) {
                wa3[j - 1] = 0.0;
                int l      = ipvt[j - 1];
                double tmp = wa1[l - 1];
                for (int i = 1; i <= j; i++)
                    wa3[i - 1] += FLD(i, j) * tmp;
            }
            double temp1  = enorm(n, wa3) / fnorm;
            double temp2  = (sqrt(par) * pnorm) / fnorm;
            double prered = temp1 * temp1 + temp2 * temp2 / p5;
            double dirder = -(temp1 * temp1 + temp2 * temp2);

            double ratio = 0.0;
            if (prered != 0.0) ratio = actred / prered;

            if (ratio <= p25) {
                double tmp;
                if (actred >= 0.0) tmp = p5;
                else tmp = p5 * dirder / (dirder + p5 * actred);
                if (p1 * fnorm1 >= fnorm || tmp < p1) tmp = p1;
                delta = tmp * (delta < pnorm / p1 ? delta : pnorm / p1);
                par  /= tmp;
            } else if (par == 0.0 || ratio >= p75) {
                delta = pnorm / p5;
                par  *= p5;
            }

            if (ratio >= p0001) {
                for (int j = 0; j < n; j++) {
                    x[j]   = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (int i = 0; i < m; i++) fvec[i] = wa4[i];
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0)
                *info = 1;
            if (delta <= xtol * xnorm) {
                if (*info == 1) *info = 3;
                else            *info = 2;
            }
            if (*info != 0) goto done_lmder;

            if (*nfev >= maxfev)                                          *info = 5;
            if (fabs(actred) <= epsmch && prered <= epsmch
                && p5 * ratio <= 1.0)                                     *info = 6;
            if (delta <= epsmch * xnorm)                                  *info = 7;
            if (gnorm <= epsmch)                                          *info = 8;
            if (*info != 0) goto done_lmder;

            if (ratio >= p0001) break; /* exit inner */
        } /* inner */
    } /* outer */

done_lmder:
    if (iflag < 0) *info = iflag;
    iflag = 0;
    if (nprint > 0) fcn(m, n, x, fvec, fjac, ldfjac, &iflag, udata);
#undef FLD
}

/* -------------------------------------------------------------------------
 * lmder1 - simplified driver for lmder
 * ------------------------------------------------------------------------- */
void minpack_lmder1(minpack_fcn_lmder fcn, int m, int n, double *x,
                    double *fvec, double *fjac, int ldfjac, double tol,
                    int *info, int *ipvt, double *wa, int lwa, void *udata)
{
    *info = 0;
    if (n <= 0 || m < n || ldfjac < m || tol < 0.0 || lwa < 5 * n + m)
        return;

    int maxfev = 100 * (n + 1);
    int nfev, njev;

    minpack_lmder(fcn, m, n, x, fvec, fjac, ldfjac,
                  tol, tol, 0.0, maxfev,
                  wa, 1, 100.0, 0,
                  info, &nfev, &njev, ipvt, wa + n,
                  wa + 2 * n, wa + 3 * n, wa + 4 * n, wa + 5 * n,
                  udata);
    if (*info == 8) *info = 4;
}

/* -------------------------------------------------------------------------
 * lmdif - Levenberg-Marquardt, numerical Jacobian
 * ------------------------------------------------------------------------- */
void minpack_lmdif(minpack_func2 fcn, int m, int n, double *x, double *fvec,
                   double ftol, double xtol, double gtol,
                   int maxfev, double epsfcn, double *diag, int mode,
                   double factor, int nprint, int *info, int *nfev,
                   double *fjac, int ldfjac, int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata)
{
#define FDF(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    static const double p1    = 1.0e-1;
    static const double p5    = 5.0e-1;
    static const double p25   = 2.5e-1;
    static const double p75   = 7.5e-1;
    static const double p0001 = 1.0e-4;

    *info = 0;
    int iflag = 0;
    *nfev = 0;

    if (n <= 0 || m < n || ldfjac < m || ftol < 0.0 || xtol < 0.0
        || gtol < 0.0 || maxfev <= 0 || factor <= 0.0)
        goto done_lmdif;
    if (mode == 2) {
        for (int j = 0; j < n; j++)
            if (diag[j] <= 0.0) goto done_lmdif;
    }

    iflag = 1;
    fcn(m, n, x, fvec, &iflag, udata);
    *nfev = 1;
    if (iflag < 0) goto done_lmdif;

    double fnorm = enorm(m, fvec);
    double par   = 0.0;
    int    iter  = 1;
    double delta = 0.0, xnorm = 0.0;

    for (;;) { /* outer */
        iflag = 2;
        fdjac2(fcn, udata, m, n, x, fvec, fjac, ldfjac, &iflag, epsfcn, wa4);
        *nfev += n;
        if (iflag < 0) goto done_lmdif;

        if (nprint > 0 && (iter - 1) % nprint == 0) {
            iflag = 0;
            fcn(m, n, x, fvec, &iflag, udata);
            if (iflag < 0) goto done_lmdif;
        }

        qrfac(m, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);

        if (iter == 1) {
            if (mode != 2) {
                for (int j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) diag[j] = 1.0;
                }
            }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if (delta == 0.0) delta = factor;
        }

        for (int i = 0; i < m; i++) wa4[i] = fvec[i];
        for (int j = 1; j <= n; j++) {
            if (FDF(j, j) != 0.0) {
                double sum = 0.0;
                for (int i = j; i <= m; i++) sum += FDF(i, j) * wa4[i - 1];
                double temp = -sum / FDF(j, j);
                for (int i = j; i <= m; i++) wa4[i - 1] += FDF(i, j) * temp;
            }
            FDF(j, j)   = wa1[j - 1];
            qtf[j - 1]  = wa4[j - 1];
        }

        double gnorm = 0.0;
        if (fnorm != 0.0) {
            for (int j = 1; j <= n; j++) {
                int l = ipvt[j - 1];
                if (wa2[l - 1] != 0.0) {
                    double sum = 0.0;
                    for (int i = 1; i <= j; i++)
                        sum += FDF(i, j) * (qtf[i - 1] / fnorm);
                    double v = fabs(sum / wa2[l - 1]);
                    if (v > gnorm) gnorm = v;
                }
            }
        }
        if (gnorm <= gtol) { *info = 4; goto done_lmdif; }

        if (mode != 2) {
            for (int j = 0; j < n; j++)
                if (wa2[j] > diag[j]) diag[j] = wa2[j];
        }

        for (;;) { /* inner */
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par,
                  wa1, wa2, wa3, wa4);

            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm(n, wa3);
            if (iter == 1) delta = delta < pnorm ? delta : pnorm;

            iflag = 1;
            fcn(m, n, wa2, wa4, &iflag, udata);
            (*nfev)++;
            if (iflag < 0) goto done_lmdif;
            double fnorm1 = enorm(m, wa4);

            double actred = -1.0;
            if (p1 * fnorm1 < fnorm)
                actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);

            for (int j = 1; j <= n; j++) {
                wa3[j - 1] = 0.0;
                int l      = ipvt[j - 1];
                double tmp = wa1[l - 1];
                for (int i = 1; i <= j; i++)
                    wa3[i - 1] += FDF(i, j) * tmp;
            }
            double temp1  = enorm(n, wa3) / fnorm;
            double temp2  = (sqrt(par) * pnorm) / fnorm;
            double prered = temp1 * temp1 + temp2 * temp2 / p5;
            double dirder = -(temp1 * temp1 + temp2 * temp2);

            double ratio = 0.0;
            if (prered != 0.0) ratio = actred / prered;

            if (ratio <= p25) {
                double tmp;
                if (actred >= 0.0) tmp = p5;
                else tmp = p5 * dirder / (dirder + p5 * actred);
                if (p1 * fnorm1 >= fnorm || tmp < p1) tmp = p1;
                delta = tmp * (delta < pnorm / p1 ? delta : pnorm / p1);
                par  /= tmp;
            } else if (par == 0.0 || ratio >= p75) {
                delta = pnorm / p5;
                par  *= p5;
            }

            if (ratio >= p0001) {
                for (int j = 0; j < n; j++) {
                    x[j]   = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (int i = 0; i < m; i++) fvec[i] = wa4[i];
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0)
                *info = 1;
            if (delta <= xtol * xnorm) {
                if (*info == 1) *info = 3;
                else            *info = 2;
            }
            if (*info != 0) goto done_lmdif;

            if (*nfev >= maxfev)                                          *info = 5;
            if (fabs(actred) <= epsmch && prered <= epsmch
                && p5 * ratio <= 1.0)                                     *info = 6;
            if (delta <= epsmch * xnorm)                                  *info = 7;
            if (gnorm <= epsmch)                                          *info = 8;
            if (*info != 0) goto done_lmdif;

            if (ratio >= p0001) break;
        } /* inner */
    } /* outer */

done_lmdif:
    if (iflag < 0) *info = iflag;
    iflag = 0;
    if (nprint > 0) fcn(m, n, x, fvec, &iflag, udata);
#undef FDF
}

/* -------------------------------------------------------------------------
 * lmdif1 - simplified driver for lmdif
 * ------------------------------------------------------------------------- */
void minpack_lmdif1(minpack_func2 fcn, int m, int n, double *x, double *fvec,
                    double tol, int *info, int *iwa, double *wa, int lwa,
                    void *udata)
{
    *info = 0;
    if (n <= 0 || m < n || tol < 0.0 || lwa < m * n + 5 * n + m) return;

    int maxfev = 200 * (n + 1);
    int mp5n   = m + 5 * n;
    int nfev;

    minpack_lmdif(fcn, m, n, x, fvec,
                  tol, tol, 0.0, maxfev, 0.0,
                  wa, 1, 100.0, 0,
                  info, &nfev,
                  wa + mp5n, m, iwa, wa + n,
                  wa + 2 * n, wa + 3 * n, wa + 4 * n, wa + 5 * n,
                  udata);
    if (*info == 8) *info = 4;
}

/* -------------------------------------------------------------------------
 * lmstr - Levenberg-Marquardt, user provides one Jacobian row at a time
 * ------------------------------------------------------------------------- */
void minpack_lmstr(minpack_fcn_lmstr fcn, int m, int n, double *x,
                   double *fvec, double *fjac, int ldfjac,
                   double ftol, double xtol, double gtol,
                   int maxfev, double *diag, int mode, double factor,
                   int nprint, int *info, int *nfev, int *njev,
                   int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata)
{
#define FLS(i,j)  fjac[((j)-1)*(ldfjac)+((i)-1)]
    static const double p1    = 1.0e-1;
    static const double p5    = 5.0e-1;
    static const double p25   = 2.5e-1;
    static const double p75   = 7.5e-1;
    static const double p0001 = 1.0e-4;

    *info = 0;
    int iflag = 0;
    *nfev = 0;
    *njev = 0;

    if (n <= 0 || m < n || ldfjac < n || ftol < 0.0 || xtol < 0.0
        || gtol < 0.0 || maxfev <= 0 || factor <= 0.0)
        goto done_lmstr;
    if (mode == 2) {
        for (int j = 0; j < n; j++)
            if (diag[j] <= 0.0) goto done_lmstr;
    }

    iflag = 1;
    fcn(m, n, x, fvec, wa3, &iflag, udata);
    *nfev = 1;
    if (iflag < 0) goto done_lmstr;
    double fnorm = enorm(m, fvec);
    double par   = 0.0;
    int    iter  = 1;
    double delta = 0.0, xnorm = 0.0;

    for (;;) { /* outer */
        if (nprint > 0 && (iter - 1) % nprint == 0) {
            iflag = 0;
            fcn(m, n, x, fvec, wa3, &iflag, udata);
            if (iflag < 0) goto done_lmstr;
        }

        /* build QR row-by-row */
        for (int j = 1; j <= n; j++) {
            qtf[j - 1] = 0.0;
            for (int i = 1; i <= n; i++) FLS(i, j) = 0.0;
        }
        iflag = 2;
        for (int i = 1; i <= m; i++) {
            fcn(m, n, x, fvec, wa3, &iflag, udata);
            if (iflag < 0) goto done_lmstr;
            double temp = fvec[i - 1];
            rwupdt(n, fjac, ldfjac, wa3, qtf, &temp, wa1, wa2);
            iflag++;
        }
        (*njev)++;

        /* handle rank-deficient case */
        int sing = 0;
        for (int j = 1; j <= n; j++) {
            if (FLS(j, j) == 0.0) sing = 1;
            ipvt[j - 1] = j;
            wa2[j - 1]  = enorm(j, &fjac[(j - 1) * ldfjac]);
        }
        if (sing) {
            qrfac(n, n, fjac, ldfjac, 1, ipvt, n, wa1, wa2, wa3);
            for (int j = 1; j <= n; j++) {
                if (FLS(j, j) != 0.0) {
                    double sum = 0.0;
                    for (int i = j; i <= n; i++) sum += FLS(i, j) * qtf[i - 1];
                    double temp = -sum / FLS(j, j);
                    for (int i = j; i <= n; i++) qtf[i - 1] += FLS(i, j) * temp;
                }
                FLS(j, j) = wa1[j - 1];
            }
        }

        if (iter == 1) {
            if (mode != 2) {
                for (int j = 0; j < n; j++) {
                    diag[j] = wa2[j];
                    if (wa2[j] == 0.0) diag[j] = 1.0;
                }
            }
            for (int j = 0; j < n; j++) wa3[j] = diag[j] * x[j];
            xnorm = enorm(n, wa3);
            delta = factor * xnorm;
            if (delta == 0.0) delta = factor;
        }

        double gnorm = 0.0;
        if (fnorm != 0.0) {
            for (int j = 1; j <= n; j++) {
                int l = ipvt[j - 1];
                if (wa2[l - 1] != 0.0) {
                    double sum = 0.0;
                    for (int i = 1; i <= j; i++)
                        sum += FLS(i, j) * (qtf[i - 1] / fnorm);
                    double v = fabs(sum / wa2[l - 1]);
                    if (v > gnorm) gnorm = v;
                }
            }
        }
        if (gnorm <= gtol) { *info = 4; goto done_lmstr; }

        if (mode != 2) {
            for (int j = 0; j < n; j++)
                if (wa2[j] > diag[j]) diag[j] = wa2[j];
        }

        for (;;) { /* inner */
            lmpar(n, fjac, ldfjac, ipvt, diag, qtf, delta, &par,
                  wa1, wa2, wa3, wa4);

            for (int j = 0; j < n; j++) {
                wa1[j] = -wa1[j];
                wa2[j] = x[j] + wa1[j];
                wa3[j] = diag[j] * wa1[j];
            }
            double pnorm = enorm(n, wa3);
            if (iter == 1) delta = delta < pnorm ? delta : pnorm;

            iflag = 1;
            fcn(m, n, wa2, wa4, wa3, &iflag, udata);
            (*nfev)++;
            if (iflag < 0) goto done_lmstr;
            double fnorm1 = enorm(m, wa4);

            double actred = -1.0;
            if (p1 * fnorm1 < fnorm)
                actred = 1.0 - (fnorm1 / fnorm) * (fnorm1 / fnorm);

            for (int j = 1; j <= n; j++) {
                wa3[j - 1] = 0.0;
                int l      = ipvt[j - 1];
                double tmp = wa1[l - 1];
                for (int i = 1; i <= j; i++)
                    wa3[i - 1] += FLS(i, j) * tmp;
            }
            double temp1  = enorm(n, wa3) / fnorm;
            double temp2  = (sqrt(par) * pnorm) / fnorm;
            double prered = temp1 * temp1 + temp2 * temp2 / p5;
            double dirder = -(temp1 * temp1 + temp2 * temp2);

            double ratio = 0.0;
            if (prered != 0.0) ratio = actred / prered;

            if (ratio <= p25) {
                double tmp;
                if (actred >= 0.0) tmp = p5;
                else tmp = p5 * dirder / (dirder + p5 * actred);
                if (p1 * fnorm1 >= fnorm || tmp < p1) tmp = p1;
                delta = tmp * (delta < pnorm / p1 ? delta : pnorm / p1);
                par  /= tmp;
            } else if (par == 0.0 || ratio >= p75) {
                delta = pnorm / p5;
                par  *= p5;
            }

            if (ratio >= p0001) {
                for (int j = 0; j < n; j++) {
                    x[j]   = wa2[j];
                    wa2[j] = diag[j] * x[j];
                }
                for (int i = 0; i < m; i++) fvec[i] = wa4[i];
                xnorm = enorm(n, wa2);
                fnorm = fnorm1;
                iter++;
            }

            if (fabs(actred) <= ftol && prered <= ftol && p5 * ratio <= 1.0)
                *info = 1;
            if (delta <= xtol * xnorm) {
                if (*info == 1) *info = 3;
                else            *info = 2;
            }
            if (*info != 0) goto done_lmstr;

            if (*nfev >= maxfev)                                          *info = 5;
            if (fabs(actred) <= epsmch && prered <= epsmch
                && p5 * ratio <= 1.0)                                     *info = 6;
            if (delta <= epsmch * xnorm)                                  *info = 7;
            if (gnorm <= epsmch)                                          *info = 8;
            if (*info != 0) goto done_lmstr;

            if (ratio >= p0001) break;
        } /* inner */
    } /* outer */

done_lmstr:
    if (iflag < 0) *info = iflag;
    iflag = 0;
    if (nprint > 0) fcn(m, n, x, fvec, wa3, &iflag, udata);
#undef FLS
}

/* -------------------------------------------------------------------------
 * lmstr1 - simplified driver for lmstr
 * ------------------------------------------------------------------------- */
void minpack_lmstr1(minpack_fcn_lmstr fcn, int m, int n, double *x,
                    double *fvec, double *fjac, int ldfjac, double tol,
                    int *info, int *ipvt, double *wa, int lwa, void *udata)
{
    *info = 0;
    if (n <= 0 || m < n || ldfjac < n || tol < 0.0 || lwa < 5 * n + m)
        return;

    int maxfev = 100 * (n + 1);
    int nfev, njev;

    minpack_lmstr(fcn, m, n, x, fvec, fjac, ldfjac,
                  tol, tol, 0.0, maxfev,
                  wa, 1, 100.0, 0,
                  info, &nfev, &njev, ipvt, wa + n,
                  wa + 2 * n, wa + 3 * n, wa + 4 * n, wa + 5 * n,
                  udata);
    if (*info == 8) *info = 4;
}
