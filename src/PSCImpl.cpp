//
//  PSC_OLS.cpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-24.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <cstring>
#include <cfloat>

#include <Rcpp.h>

#include "config.h"
#include "PSC.hpp"
#include "BLAS.h"
#include "LapackException.hpp"

static const int BLAS_M1L = -1;
static const int BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;
static const double BLAS_M1F = -1.0;

static const double LAPACK_EV_ABSTOL = 1e-18;
static const double LAPACK_EV_RANGE_LOWER = 1e-12;
static const double LAPACK_EV_RANGE_UPPER = DBL_MAX;

static const char * const BLAS_SIDE_LEFT = "L";
static const char * const BLAS_UPLO_UPPER = "U";
static const char * const BLAS_TRANS_NO = "N";
static const char * const BLAS_TRANS_TRANS = "T";
static const char * const BLAS_DIAG_NO = "N";


#ifdef DEBUG
static void print_matf(int dr, int dc, const double *A, const char *header) {
    int r, c, i = 0;
    Rcpp::Rcout << "\n" << header << ":" << std::endl;
    for (r = 0; r < dr; ++r) {
        for (c = 0, i = r; c < dc; ++c, i += dr) {
            Rprintf("%7.4f ", A[i]);
        }
        Rprintf("\n");
    }
}
#endif

PSC::PSC()
{
    this->EVDworkMemSize = 1;
    this->EVIworkMemSize = 1;

    this->EVDworkMem = new double[this->EVDworkMemSize];
    this->EVIworkMem = new int[this->EVIworkMemSize];
    this->evalues = new double[1];
    this->evectors = new double[1];
    this->evectorsSupport = new int[1];
}

PSC::~PSC()
{
    delete[] this->EVDworkMem;
    delete[] this->EVIworkMem;
    delete[] evalues;
    delete[] evectors;
    delete[] evectorsSupport;
}

void PSC::setData(const Data &data)
{
    int nvar = data.numVar();
    if (nvar > this->data.numVar()) {
        delete[] evalues;
        delete[] evectors;
        delete[] evectorsSupport;

        this->evalues = new double[nvar];
        this->evectors = new double[nvar * nvar];
        this->evectorsSupport = new int[2 * nvar];
    }

    this->data = data;
}

int PSC::doEigenDecomposition(const char* const uplo, double *RESTRICT matrix, const int n)
{
    int nevalues, lapackInfo;
    /* Compute needed size of working space */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, this->evalues,
                        this->evectors, n, this->evectorsSupport,
                        this->EVDworkMem, BLAS_M1L,
                        this->EVIworkMem, BLAS_M1L, lapackInfo);

    if (lapackInfo != 0) {
        throw LapackException("Eigenvalue decomposition workspace query failed", lapackInfo);
    }

    if (this->EVDworkMem[0] > this->EVDworkMemSize) {
        this->EVDworkMemSize = (int) this->EVDworkMem[0];
        delete[] this->EVDworkMem;
        this->EVDworkMem = new double[this->EVDworkMemSize];
    }

    if (this->EVIworkMem[0] > this->EVIworkMemSize) {
        this->EVIworkMemSize = (int) this->EVIworkMem[0];
        delete[] this->EVIworkMem;
        this->EVIworkMem = new int[this->EVIworkMemSize];
    }

    /* Actually perform eigenvalue decomposition */
    LAPACK_DSYEVR_RANGE(uplo, n, matrix, n,
                        LAPACK_EV_RANGE_LOWER, LAPACK_EV_RANGE_UPPER, LAPACK_EV_ABSTOL,
                        nevalues, this->evalues,
                        this->evectors, n, this->evectorsSupport,
                        this->EVDworkMem, this->EVDworkMemSize,
                        this->EVIworkMem, this->EVIworkMemSize, lapackInfo);


    if (lapackInfo != 0) {
        throw LapackException("Eigenvalue decomposition failed", lapackInfo);
    }

    return nevalues;
}



PSC_OLS::PSC_OLS() : gotResiduals(FALSE), initialized(FALSE)
{}

PSC_OLS::~PSC_OLS()
{
    if (this->initialized) {
        delete[] this->Xsqrt;
        delete[] this->XsqrtInvX;
        delete[] this->Z;
        delete[] this->H;
        delete[] this->residuals;
    }
}


void PSC_OLS::setData(const Data &data) {
    if (data.numVar() > this->data.numVar()) {
        if (this->initialized) {
            delete[] this->Xsqrt;
        }
        this->Xsqrt = new double[data.numVar() * data.numVar()];
    }

    if (data.numObs() > this->data.numObs()) {
        if (this->initialized) {
            delete[] this->H;
            delete[] this->residuals;
        }
        this->H = new double[data.numObs() * data.numObs()];
        this->residuals = new double[data.numObs()];
    }

    if (data.numObs() * data.numVar() > this->data.numObs() * this->data.numVar()) {
        if (this->initialized) {
            delete[] this->Z;
            delete[] this->XsqrtInvX;
        }
        this->Z = new double[data.numObs() * data.numVar()];
        this->XsqrtInvX = new double[data.numObs() * data.numVar()];
    }

    PSC::setData(data);
    this->initialized = TRUE;
}

void PSC_OLS::setResiduals(const double *RESTRICT residuals)
{
    memcpy(this->residuals, residuals, this->data.numObs() * sizeof(double));
    this->gotResiduals = true;
}

int PSC_OLS::computePSC() {
    int lapackInfo = 0;
    int i;
    int nvar = this->data.numVar();
    int nobs = this->data.numObs();
    int nevalues = 0;
    /*
     * G can point to the same address as XtXinvX, because we don't need
     * XtXinvX anymore when G is computed.
     * The same is true for Q and H, as well as Z and XtXinvX
     */
    double *RESTRICT XtXinvX = this->Z;
    double *RESTRICT G = this->Z;
    double *RESTRICT Q = this->H;
    double *iterG, *iterH;
    double Wii;

    memcpy(XsqrtInvX, this->data.getXtr(), nobs * nvar * sizeof(double));

    /* Xsqrt = X %*% t(X) */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
        BLAS_1F, this->data.getXtr(), nvar, this->data.getXtr(), nvar,
        BLAS_0F, this->Xsqrt, nvar);

    /* Xsqrt = chol(Xsqrt) */
    LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, Xsqrt, nvar, lapackInfo);

    if (lapackInfo != 0) {
        throw LapackException("Could not compute Cholesky decomposition.", lapackInfo);
    }

    /* t(XsqrtInvX) = t(inv(Xsqrt)) %*% t(X) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, this->Xsqrt, nvar, this->XsqrtInvX, nvar);

    memcpy(XtXinvX, this->XsqrtInvX, nobs * nvar * sizeof(double));

    /* XtXinvX = inv(Xsqrt) %*% t(XsqrtInvX) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, this->Xsqrt, nvar, XtXinvX, nvar);

    /* H = t(X) %*% XtXinvX */
    BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nobs, nvar,
               BLAS_1F, this->data.getXtr(), nvar, XtXinvX, nvar,
               BLAS_0F, this->H, nobs);

    /* Did we already get the residuals? */
    if (!this->gotResiduals) {
        memcpy(residuals, this->data.getY(), nobs * sizeof(double));
        /* residuals = -H %*% y + y --> H is symmetric! */
        BLAS_DSYMV(BLAS_UPLO_UPPER, nobs,
                   BLAS_M1F, this->H, nobs, this->data.getY(), BLAS_1L,
                   BLAS_1F, this->residuals, BLAS_1L);

    }

    memcpy(G, XsqrtInvX, nobs * nvar * sizeof(double));

    /* t(G) = t(XsqrtInvX) %*% diag(W[i, i]) */
    for (i = 0, iterH = H, iterG = G; i < nobs; ++i, iterH += nobs + 1, iterG += nvar) {
        /* W[i, i] = */
        Wii = residuals[i] / (1 - *iterH);

        /* G[ , i] = W[i, i] * XsqrtInvX[, i] */
        BLAS_DSCAL(nvar, Wii, iterG, BLAS_1L); // Do loop by hand!
    }

    /* Q = t(G) %*% G */
    BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
               BLAS_1F, G, nvar, G, nvar,
               BLAS_0F, Q, nvar);

    nevalues = this->doEigenDecomposition(BLAS_UPLO_UPPER, Q, nvar);

    if (nevalues > 0) {
        BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nevalues, nvar,
                   BLAS_1F, this->XsqrtInvX, nvar, this->getEigenvectors(), nvar,
                   BLAS_0F, this->Z, nobs);
    }

    /* Invalidate the current residuals */
    this->gotResiduals = false;

    return nevalues;
}



PSC_EN::PSC_EN() : gotResiduals(FALSE), initialized(FALSE)
{}

PSC_EN::~PSC_EN()
{
    if (this->initialized) {
        delete[] this->Z;
        delete[] this->residMat;
    }
}


void PSC_EN::setData(const Data &data) {
//    if (data.numVar() > this->data.numVar()) {
//        if (this->initialized) {
//        }
//    }

    if (data.numObs() > this->data.numObs()) {
        if (this->initialized) {
            delete[] this->residMat;
        }
        this->residMat = new double[data.numObs() * data.numObs()];
    }

    if (data.numObs() * data.numVar() > this->data.numObs() * this->data.numVar()) {
        if (this->initialized) {
            delete[] this->Z;
        }
        this->Z = new double[data.numObs() * data.numVar()];
    }

    /*
     * We will let the last column of residMat be the one where the true residuals
     * be stored
     */
    this->residuals = this->residMat + (this->data.numObs() * (this->data.numObs() - 1));

    PSC::setData(data);
    this->initialized = TRUE;
}

void PSC_EN::setResiduals(const double *RESTRICT residuals)
{
    memcpy(this->residuals, residuals, this->data.numObs() * sizeof(double));
    this->gotResiduals = true;
}

int PSC_EN::computePSC() {
    int i;
    int nvar = this->data.numVar();
    int nobs = this->data.numObs();
    int nevalues = 0;

//    nevalues = this->doEigenDecomposition(BLAS_UPLO_UPPER, Q, nvar);

    if (nevalues > 0) {
        BLAS_DGEMM(BLAS_TRANS_TRANS, BLAS_TRANS_NO, nobs, nevalues, nvar,
                   BLAS_1F, this->residMat, nvar, this->getEigenvectors(), nvar,
                   BLAS_0F, this->Z, nobs);
    }

    /* Invalidate the current residuals */
    this->gotResiduals = false;

    return nevalues;
}

