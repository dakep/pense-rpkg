//
//  olsreg.c
//  pense
//
//  Created by David Kepplinger on 2016-02-13.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include "config.h"
#include "olsreg.h"
#include "BLAS.h"

static BLAS_INT BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static BLAS_CHAR BLAS_UPLO_UPPER = "U";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_DIAG_NO = "N";

int computeOLSCoefs(const double *restrict Xtr, const double *restrict y,
                    const int nobs, const int nvar,
					double *restrict coefs, double *restrict Xsqrt)
{
    int lapackInfo;
    /* We won't handle rank-deficient cases here!! */

    /* Xsqrt = t(X) %*% X */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
        BLAS_1F, Xtr, nvar, Xtr, nvar,
        BLAS_0F, Xsqrt, nvar);

    /* Xsqrt = chol(Xsqrt) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, Xsqrt, nvar, lapackInfo);

    if (lapackInfo != 0) {
        return lapackInfo;
    }

    /* coefEst = t(X) %*% y */
    BLAS_DGEMV(BLAS_TRANS_NO, nvar, nobs, BLAS_1F, Xtr, nvar, y, BLAS_1L, BLAS_0F, coefs,
               BLAS_1L);

    /* coefEst = inv(t(chol(Xsqrt))) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO, nvar, Xsqrt, nvar, coefs,
               BLAS_1L);

    /* coefEst = inv(chol(Xsqrt)) %*% coefEst */
    BLAS_DTRSV(BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO, nvar, Xsqrt, nvar, coefs,
               BLAS_1L);

    return 0;
}


void computeResiduals(const double *restrict Xtr, const double *restrict y,
					  const int nobs, const int nvar,
					  const double *restrict coefs, double *restrict residuals)
{
    memcpy(residuals, y, nobs * sizeof(double));

    BLAS_DGEMV(BLAS_TRANS_TRANS, nvar, nobs,
               BLAS_M1F, Xtr, nvar,
               coefs, BLAS_1L,
               BLAS_1F, residuals, BLAS_1L);
}

