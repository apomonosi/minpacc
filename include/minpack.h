/* minpack.h - Public C API for Minpack nonlinear solver library
 *
 * Original Minpack: Argonne National Laboratory, 1980.
 *   Burton S. Garbow, Kenneth E. Hillstrom, Jorge J. More.
 * Fortran modernization: Jacob Williams, 2021.
 * C translation: 2024.
 *
 * All 2-D array arguments (fjac, r) use COLUMN-MAJOR storage
 * with leading dimension as specified (matching Fortran convention).
 * Element (row i, col j) of array a with leading dim lda:
 *   a[j * lda + i]   (both i and j are 0-based)
 */

#ifndef MINPACK_H
#define MINPACK_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Callback typedefs
 * All arrays passed to callbacks are 0-based C arrays.
 * For 2-D arrays (fjac), column-major storage is used.
 *
 * iflag semantics:
 *   iflag=0  : printing iteration (do not alter fvec/fjac)
 *   iflag=1  : evaluate fvec at x
 *   iflag=2  : evaluate fjac at x (or Jacobian row for lmstr)
 *   iflag=i  : (lmstr only) evaluate row i-1 of Jacobian
 *   Set *iflag < 0 to request early termination.
 */

/* For hybrd / hybrd1 / fdjac1 */
typedef void (*minpack_func)(int n, const double *x, double *fvec,
                             int *iflag, void *udata);

/* For lmdif / lmdif1 / fdjac2 */
typedef void (*minpack_func2)(int m, int n, const double *x, double *fvec,
                              int *iflag, void *udata);

/* For hybrj / hybrj1 (provides both fvec and column-major fjac) */
typedef void (*minpack_fcn_hybrj)(int n, const double *x, double *fvec,
                                  double *fjac, int ldfjac, int *iflag,
                                  void *udata);

/* For lmder / lmder1 (provides both fvec and column-major fjac) */
typedef void (*minpack_fcn_lmder)(int m, int n, const double *x, double *fvec,
                                  double *fjac, int ldfjac, int *iflag,
                                  void *udata);

/* For lmstr / lmstr1 (provides fvec and one row of Jacobian at a time) */
typedef void (*minpack_fcn_lmstr)(int m, int n, const double *x, double *fvec,
                                  double *fjrow, int *iflag, void *udata);

/*
 * minpack_dpmpar - machine constants
 *   i=1: machine epsilon (unit roundoff)
 *   i=2: smallest positive magnitude (DBL_MIN)
 *   i=3: largest magnitude (DBL_MAX)
 */
double minpack_dpmpar(int i);

/*
 * minpack_chkder - check user-supplied Jacobian for consistency
 *
 * Call twice: first with mode=1 to compute xp from x,
 * then with mode=2 to check derivatives using fvec, fjac, fvecp.
 * err[i] near 1.0 means derivative i is correct.
 */
void minpack_chkder(int m, int n, const double *x, double *fvec,
                    double *fjac, int ldfjac, double *xp, double *fvecp,
                    int mode, double *err);

/*
 * minpack_hybrd - solve n nonlinear equations in n variables
 * Uses forward-difference Jacobian approximation.
 * Workspace: r[lr], qtf[n], fjac[ldfjac*n], wa1..wa4[n]
 */
void minpack_hybrd(minpack_func fcn, int n, double *x, double *fvec,
                   double xtol, int maxfev, int ml, int mu, double epsfcn,
                   double *diag, int mode, double factor, int nprint,
                   int *info, int *nfev, double *fjac, int ldfjac,
                   double *r, int lr, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata);

/*
 * minpack_hybrd1 - simplified driver for hybrd
 * Workspace wa[lwa] where lwa >= (n*(3*n+13))/2
 */
void minpack_hybrd1(minpack_func fcn, int n, double *x, double *fvec,
                    double tol, int *info, double *wa, int lwa, void *udata);

/*
 * minpack_hybrj - solve n nonlinear equations in n variables
 * User provides both fvec and Jacobian.
 */
void minpack_hybrj(minpack_fcn_hybrj fcn, int n, double *x, double *fvec,
                   double *fjac, int ldfjac, double xtol, int maxfev,
                   double *diag, int mode, double factor, int nprint,
                   int *info, int *nfev, int *njev,
                   double *r, int lr, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata);

/*
 * minpack_hybrj1 - simplified driver for hybrj
 * Workspace wa[lwa] where lwa >= (n*(n+13))/2
 */
void minpack_hybrj1(minpack_fcn_hybrj fcn, int n, double *x, double *fvec,
                    double *fjac, int ldfjac, double tol, int *info,
                    double *wa, int lwa, void *udata);

/*
 * minpack_lmder - minimize sum of squares of m functions in n variables
 * User provides both fvec and Jacobian.
 */
void minpack_lmder(minpack_fcn_lmder fcn, int m, int n, double *x,
                   double *fvec, double *fjac, int ldfjac,
                   double ftol, double xtol, double gtol,
                   int maxfev, double *diag, int mode, double factor,
                   int nprint, int *info, int *nfev, int *njev,
                   int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata);

/*
 * minpack_lmder1 - simplified driver for lmder
 * Workspace wa[lwa] where lwa >= 5*n + m
 */
void minpack_lmder1(minpack_fcn_lmder fcn, int m, int n, double *x,
                    double *fvec, double *fjac, int ldfjac, double tol,
                    int *info, int *ipvt, double *wa, int lwa, void *udata);

/*
 * minpack_lmdif - minimize sum of squares of m functions in n variables
 * Uses forward-difference Jacobian approximation.
 */
void minpack_lmdif(minpack_func2 fcn, int m, int n, double *x, double *fvec,
                   double ftol, double xtol, double gtol,
                   int maxfev, double epsfcn, double *diag, int mode,
                   double factor, int nprint, int *info, int *nfev,
                   double *fjac, int ldfjac, int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata);

/*
 * minpack_lmdif1 - simplified driver for lmdif
 * Workspace: iwa[n] integer, wa[lwa] where lwa >= m*n + 5*n + m
 */
void minpack_lmdif1(minpack_func2 fcn, int m, int n, double *x, double *fvec,
                    double tol, int *info, int *iwa, double *wa, int lwa,
                    void *udata);

/*
 * minpack_lmstr - minimize sum of squares, user provides one Jacobian row at a time
 * (minimal storage variant)
 */
void minpack_lmstr(minpack_fcn_lmstr fcn, int m, int n, double *x,
                   double *fvec, double *fjac, int ldfjac,
                   double ftol, double xtol, double gtol,
                   int maxfev, double *diag, int mode, double factor,
                   int nprint, int *info, int *nfev, int *njev,
                   int *ipvt, double *qtf,
                   double *wa1, double *wa2, double *wa3, double *wa4,
                   void *udata);

/*
 * minpack_lmstr1 - simplified driver for lmstr
 * Workspace wa[lwa] where lwa >= 5*n + m
 */
void minpack_lmstr1(minpack_fcn_lmstr fcn, int m, int n, double *x,
                    double *fvec, double *fjac, int ldfjac, double tol,
                    int *info, int *ipvt, double *wa, int lwa, void *udata);

#ifdef __cplusplus
}
#endif

#endif /* MINPACK_H */
