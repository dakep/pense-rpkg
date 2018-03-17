//
//  PSC_OLS.cpp
//  pense
//
//  Created by David Kepplinger on 2016-01-24.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include "config.h"

#include <cstring>
#include <cfloat>

#include "BLAS.h"
#include <RcppArmadillo.h>

#include "PSCxx.hpp"
#include "LapackException.hpp"

static BLAS_INT BLAS_M1L = -1;
static BLAS_INT BLAS_1L = 1;
static const double BLAS_0F = 0.0;
static const double BLAS_1F = 1.0;

static const double LAPACK_EV_ABSTOL = NUMERIC_EPS;
static const double LAPACK_EV_RANGE_LOWER = LAPACK_EV_MIN;
static const double LAPACK_EV_RANGE_UPPER = DBL_MAX;

static BLAS_CHAR BLAS_SIDE_LEFT = "L";
static BLAS_CHAR BLAS_UPLO_UPPER = "U";
static BLAS_CHAR BLAS_TRANS_NO = "N";
static BLAS_CHAR BLAS_TRANS_TRANS = "T";
static BLAS_CHAR BLAS_DIAG_NO = "N";


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

PSC::PSC(bool nobsEigenvalues) : nobsEigenvalues(nobsEigenvalues)
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
    const int maxEVs = (this->nobsEigenvalues ? data.numObs() : data.numVar());
    const int compVal = (this->nobsEigenvalues ? this->data.numObs() : this->data.numVar());

    if (maxEVs > compVal) {
        delete[] evalues;
        delete[] evectors;
        delete[] evectorsSupport;

        this->evalues = new double[maxEVs];
        this->evectors = new double[maxEVs * maxEVs];
        this->evectorsSupport = new int[2 * maxEVs];
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

/***************************************************************************************************
 *
 * PSCs for OLS using identities from the paper
 *
 **************************************************************************************************/


PSC_OLS::PSC_OLS() : PSC(FALSE), XsqrtProvided(FALSE), initialized(FALSE)
{}

PSC_OLS::~PSC_OLS()
{
    if (this->initialized) {
        if (!this->XsqrtProvided) {
            delete[] this->Xsqrt;
        }
        delete[] this->XsqrtInvX;
        delete[] this->Z;
        delete[] this->Q;
        delete[] this->residuals;
    }
}

void PSC_OLS::setXsqrtMemory(double *RESTRICT Xsqrt)
{
    if (!this->XsqrtProvided && this->initialized) {
        delete[] this->Xsqrt;
    }
    this->Xsqrt = Xsqrt;
    this->XsqrtProvided = TRUE;
}


void PSC_OLS::setData(const Data &data) {
    if (!this->XsqrtProvided && (data.numVar() > this->data.numVar())) {
        if (this->initialized) {
            delete[] this->Xsqrt;
        }
        this->Xsqrt = new double[data.numVar() * data.numVar()];
    }

    if (data.numObs() > this->data.numObs()) {
        if (this->initialized) {
            delete[] this->Q;
            delete[] this->residuals;
        }
        this->Q = new double[data.numVar() * data.numVar()];
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
}

int PSC_OLS::computePSC() {
    BLAS_INT nvar = this->data.numVar();
    BLAS_INT nobs = this->data.numObs();
    int i, j, lapackInfo;
    int nevalues = 0;
    /*
     * G can point to the same address as XtXinvX, because we don't need
     * XtXinvX anymore when G is computed.
     * The same is true for Q and H, as well as Z and XtXinvX
     */
    double *RESTRICT XtXinvX = this->Z;
    double *RESTRICT G = this->Z;
    double *RESTRICT iterG;
    const double *RESTRICT iterX;
    const double *RESTRICT iterXtXinvX;
    const double *RESTRICT iterXsqrtInvX;
    double Wii;

    memcpy(XsqrtInvX, this->data.getXtr(), nobs * nvar * sizeof(double));

    /* Xsqrt = sqrt(Xtr . t(Xtr)) = sqrt(t(X) . X) */
    if (!this->XsqrtProvided) {
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_TRANS, nvar, nvar, nobs,
                   BLAS_1F, this->data.getXtrConst(), nvar, this->data.getXtr(), nvar,
                   BLAS_0F, this->Xsqrt, nvar);

        /* Xsqrt = chol(Xsqrt) */
        LAPACK_DPOTRF(BLAS_UPLO_UPPER, nvar, Xsqrt, nvar, lapackInfo);

        if (lapackInfo != 0) {
            throw LapackException("Could not compute Cholesky decomposition.", lapackInfo);
        }
    }
    /* We already have this one! */

    /* t(XsqrtInvX) = t(inv(Xsqrt)) %*% t(X) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, this->Xsqrt, nvar, this->XsqrtInvX, nvar);

    memcpy(XtXinvX, this->XsqrtInvX, nobs * nvar * sizeof(double));

    /* XtXinvX = inv(Xsqrt) %*% t(XsqrtInvX) */
    BLAS_DTRSM(BLAS_SIDE_LEFT, BLAS_UPLO_UPPER, BLAS_TRANS_NO, BLAS_DIAG_NO,
               nvar, nobs, BLAS_1F, this->Xsqrt, nvar, XtXinvX, nvar);


    /* t(G) = t(XsqrtInvX) %*% diag(W[i, i]) */
    /* NOTE: G and XtXinvX actually point to the same memory! */
    iterG = G;
    iterX = this->data.getXtrConst();
    iterXtXinvX = XtXinvX;
    iterXsqrtInvX = this->XsqrtInvX;

    for (i = 0; i < nobs; ++i) {
        /*
         * Set Wii to i-th diagonal element of H,
         * i.e., inner product of x_i and i-th column of XtXinvX
         */
        Wii = 0;
        for (j = 0; j < nvar; ++j, ++iterX, ++iterXtXinvX) {
            Wii += (*iterX) * (*iterXtXinvX);
        }

        /* W[i, i] = */
        Wii = this->residuals[i] / (1 - Wii);

        /* G[ , i] = W[i, i] * XsqrtInvX[, i] */
        for (j = 0; j < nvar; ++j, ++iterG, ++iterXsqrtInvX) {
            (*iterG) = Wii * (*iterXsqrtInvX);
        }
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

    return nevalues;
}



/***************************************************************************************************
 *
 * Manually compute PSCs for elastic net using the residuals and LOO-residuals
 *
 **************************************************************************************************/

PSC_EN::PSC_EN(ElasticNet &en) : PSC(TRUE), initialized(FALSE), en(en)
{}

PSC_EN::~PSC_EN()
{
    if (this->initialized) {
        delete[] this->buffer;
        delete[] this->Z;
        delete[] this->residMat;
        delete[] this->residuals;
    }
}


void PSC_EN::setData(const Data &data) {
    if (data.numVar() > this->data.numVar()) {
        if (this->initialized) {
            delete[] this->buffer;
        }

        this->buffer = new double[2 * data.numVar()];
    }

    if (data.numObs() > this->data.numObs()) {
        if (this->initialized) {
            delete[] this->residMat;
            delete[] this->Z;
            delete[] this->residuals;
        }
        this->residMat = new double[data.numObs() * data.numObs()];
        this->Z = new double[data.numObs() * data.numObs()];
        this->residuals = new double [data.numObs()];
    }

    PSC::setData(data);
    this->initialized = TRUE;
}

void PSC_EN::setResiduals(const double *RESTRICT residuals)
{
    memcpy(this->residuals, residuals, this->data.numObs() * sizeof(double));
}


int PSC_EN::computePSC() {
    const int nvar = this->data.numVar();
    const int nobs = this->data.numObs();
    int i, j;
    double *RESTRICT residualsOmitted;
    double *RESTRICT coefs = this->buffer;
    double *RESTRICT swapbuffer = this->buffer + nvar;
    double *RESTRICT iterX;
    double *RESTRICT iterY;
    double *RESTRICT lastCol = this->data.getXtr() + nvar * (nobs - 1);
    double *RESTRICT lastResponse = this->data.getY() + (nobs - 1);
    int nevalues = 0;
    int accuStatus = 0;

    memset(coefs, 0, nvar * sizeof(double));

    this->data.setNumObs(nobs - 1);

    residualsOmitted = this->residMat;
    iterX = this->data.getXtr();
    iterY = this->data.getY();

    for (i = 0; i < nobs - 1; ++i, residualsOmitted += nobs, ++iterY) {
        /* Swap i-th column with last column */
        memcpy(swapbuffer, iterX, nvar * sizeof(double));
        memcpy(iterX, lastCol, nvar * sizeof(double));
        memcpy(lastCol, swapbuffer, nvar * sizeof(double));

        /* Swap last response back with i-th response */
        *swapbuffer = (*iterY);
        *iterY = (*lastResponse);
        *lastResponse = *swapbuffer;

        this->en.setData(this->data);
        this->en.computeCoefs(coefs, residualsOmitted);
        accuStatus += this->en.getStatus();

        for (j = 0; j < i; ++j) {
            /* this actually calculates -r_i, but that's okay */
            residualsOmitted[j] -= this->residuals[j];
        }
        for (j = i + 1; j < nobs - 1; ++j) {
            /* this actually calculates -r_i, but that's okay */
            residualsOmitted[j] -= this->residuals[j];
        }

        /*
         * For the i'th value, we have to take the residual from the last observation!
         * This saves the need to swap the residuals as well.
         */
        residualsOmitted[i] -= this->residuals[nobs - 1];

        /* The last residual has to be computed by hand */
        residualsOmitted[nobs - 1] = (*lastResponse) - this->residuals[i];
        for (j = 0; j < nvar; ++j) {
            residualsOmitted[nobs - 1] -= coefs[j] * lastCol[j];
        }

        /* Swap last column back with i-th column */
        memcpy(swapbuffer, iterX, nvar * sizeof(double));
        memcpy(iterX, lastCol, nvar * sizeof(double));
        memcpy(lastCol, swapbuffer, nvar * sizeof(double));

        /* Swap last response back with i-th response */
        *swapbuffer = (*iterY);
        *iterY = (*lastResponse);
        *lastResponse = *swapbuffer;

        /* Also swap the two residuals back! */
        *swapbuffer = residualsOmitted[nobs - 1];
        residualsOmitted[nobs - 1] = residualsOmitted[i];
        residualsOmitted[i] = *swapbuffer;

        iterX += nvar;
    }

    /*
     * For the last observation, we don't need to swap anything
     */
    this->en.setData(this->data);
    this->en.computeCoefs(coefs, residualsOmitted);
    accuStatus += this->en.getStatus();

    for (j = 0; j < nobs - 1; ++j) {
        /* this actually calculates -r_i, but that's okay */
        residualsOmitted[j] -= this->residuals[j];
    }

    /* The last residual has to be computed by hand */
    residualsOmitted[nobs - 1] = (*lastResponse) - this->residuals[j];
    for (j = 0; j < nvar; ++j) {
        residualsOmitted[nobs - 1] -= coefs[j] * lastCol[j];
    }

    /* Reset the number of observations for data */
    this->data.setNumObs(nobs);

    /* Calculate matrix product Z = R'*R */
    BLAS_DSYRK(BLAS_UPLO_UPPER, BLAS_TRANS_TRANS, nobs, nobs,
               BLAS_1F, this->residMat, nobs, BLAS_0F, this->Z, nobs);

    nevalues = this->doEigenDecomposition(BLAS_UPLO_UPPER, this->Z, nobs);

    if (nevalues > 0) {
        BLAS_DGEMM(BLAS_TRANS_NO, BLAS_TRANS_NO, nobs, nevalues, nobs,
                   BLAS_1F, this->residMat, nobs, this->getEigenvectors(), nobs,
                   BLAS_0F, this->Z, nobs);
    }

    if (accuStatus > 0) {
        Rcpp::warning("Elastic net had non-successful status codes. PSCs may not be reliable!");
    }

    return nevalues;
}
