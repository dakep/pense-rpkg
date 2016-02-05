/*
 * BLAS.h
 * penseinit
 *
 * Created by David Kepplinger on 2016-22-01.
 * Copyright (c) 2016 David Kepplinger. All rights reserved.
 *
 * This file is part of the R package penseinit.
 *
 * penseinit is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * penseinit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with R. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef penseinit_BLAS_h
#define penseinit_BLAS_h

#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>

#define BLAS_DGEMV(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)							\
	F77_NAME(dgemv)(trans, &m, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DSYMV(uplo, n, alpha, a, lda, x, incx, beta, y, incy)								\
	F77_NAME(dsymv)(uplo, &n, &alpha, a, &lda, x, &incx, &beta, y, &incy)

#define BLAS_DTRSV(uplo, trans, diag, n, a, lda, x, incx)                                       \
    F77_NAME(dtrsv)(uplo, trans, diag, &n, a, &lda, x, &incx)

#define BLAS_DGEMM(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)				\
	F77_NAME(dgemm)(transa, transb, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, c, &ldc)

#define BLAS_DDOT(n, x, incx, y, incy)															\
	F77_NAME(ddot)(&n, x, &incx, y, &incy)

#define BLAS_DSCAL(n, alpha, x, incx)															\
	F77_NAME(dscal)(&n, &alpha, x, &incx)

#define BLAS_DTRSM(side, uplo, transa, diag, m, n, alpha, a, lda, b, ldb)                       \
    F77_NAME(dtrsm)(side, uplo, transa, diag, &m, &n, &alpha, a, &lda, b, &ldb)


#define LAPACK_DPOTRF(uplo, n, a, lda, info)                                                    \
    F77_NAME(dpotrf)(uplo, &n, a, &lda, &info)

#define LAPACK_DSYEVR_ALL_VECTORS(uplo, n, a, lda, abstol, nevalues, evalues, evectors,         \
                                  ldevectors, isuppevectors, work, lwork, iwork, liwork, info)  \
    F77_NAME(dsyevr)("V", "A", uplo, &n, a, &lda, NULL, NULL, NULL, NULL, &abstol,              \
                     &nevalues, evalues, evectors, &ldevectors, isuppevectors, work, &lwork,    \
                     iwork, &liwork, &info)


#define LAPACK_DSYEVR_RANGE(uplo, n, a, lda, rangeLower, rangeUpper, abstol, nevalues,          \
                            evalues, evectors, ldevectors, isuppevectors, work, lwork, iwork,   \
                            liwork, info)  \
    F77_NAME(dsyevr)("V", "V", uplo, &n, a, &lda, &rangeLower, &rangeUpper, NULL, NULL,         \
                     &abstol, &nevalues, evalues, evectors, &ldevectors, isuppevectors, work,   \
                     &lwork, iwork, &liwork, &info)


#define LAPACK_DGELS(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)                      \
    F77_NAME(dgels)(trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info)


#endif
