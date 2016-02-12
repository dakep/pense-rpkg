//
//  AuxMemory.c
//  penseinit
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#include <stdlib.h>

#include "AuxMemory.h"

void initAuxMemory(AuxMemory* auxmem)
{
    auxmem->dblWorkMemSize = 0;
    auxmem->intWorkMemSize = 0;
    auxmem->nvar = 0;
    auxmem->nobs = 0;
}

void freeAuxMemory(AuxMemory* auxmem)
{
    if (auxmem->nvar > 0) {
        free(auxmem->Xsqrt);
        free(auxmem->evalues);
        free(auxmem->evectorsSupport);
    }

    if (auxmem->nobs > 0) {
        free(auxmem->XsqrtInvX);
        free(auxmem->Q);
        free(auxmem->residuals);
    }

    if (auxmem->dblWorkMemSize > 0) {
        free(auxmem->dblWorkMem);
    }

    if (auxmem->intWorkMemSize > 0) {
        free(auxmem->intWorkMem);
    }
}

void resizeAuxMemory(AuxMemory* auxmem, const int nvar, const int nobs)
{
    if (nobs * nvar > auxmem->nobs * auxmem->nvar) {
        if (auxmem->nobs > 0 && auxmem->nvar > 0) {
            free(auxmem->XsqrtInvX);
            free(auxmem->Q);
        }

        auxmem->XsqrtInvX = (double*) malloc(nobs * nvar * sizeof(double));
        auxmem->Q = (double*) malloc(nobs * nvar * sizeof(double));
    }

    if (nvar > auxmem->nvar) {
        if (auxmem->nvar > 0) {
            free(auxmem->Xsqrt);
            free(auxmem->evalues);
            free(auxmem->evectorsSupport);
        }

        auxmem->nvar = nvar;

        auxmem->Xsqrt = (double*) malloc(nvar * nvar * sizeof(double));
        auxmem->evalues = (double*) malloc(nvar * sizeof(double));
        auxmem->evectorsSupport = (int*) malloc(2 * nvar * sizeof(int));
        
        auxmem->eigenvectors = auxmem->Xsqrt;
    }

    if (nobs > auxmem->nobs) {
        if (auxmem->nobs > 0) {
            free(auxmem->residuals);
        }

        auxmem->nobs = nobs;

        auxmem->residuals = (double*) malloc(nobs * sizeof(double));
    }

    if (auxmem->dblWorkMemSize == 0) {
        auxmem->dblWorkMem = (double*) malloc(sizeof(double));
        auxmem->dblWorkMemSize = 1;
    }

    if (auxmem->intWorkMemSize == 0) {
        auxmem->intWorkMem = (int*) malloc(sizeof(int));
        auxmem->intWorkMemSize = 1;
    }
}

void resizeDblWorkAuxMemory(AuxMemory* auxmem, const int newDblWorkMemSize)
{
    if (newDblWorkMemSize > auxmem->dblWorkMemSize) {
        if (auxmem->dblWorkMemSize > 0) {
            free(auxmem->dblWorkMem);
        }
        auxmem->dblWorkMemSize = newDblWorkMemSize;
        auxmem->dblWorkMem = (double*) malloc(auxmem->dblWorkMemSize * sizeof(double));
    }
}


void resizeIntWorkAuxMemory(AuxMemory* auxmem, const int newIntWorkMemSize)
{
    if (newIntWorkMemSize > auxmem->intWorkMemSize) {
        if (auxmem->intWorkMemSize > 0) {
            free(auxmem->intWorkMem);
        }
        auxmem->intWorkMemSize = newIntWorkMemSize;
        auxmem->intWorkMem = (int*) malloc(auxmem->intWorkMemSize * sizeof(int));
    }
}


