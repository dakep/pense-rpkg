//
//  AuxMemory.h
//  pense
//
//  Created by David Kepplinger on 2016-01-30.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef AuxMemory_h
#define AuxMemory_h

typedef struct _AuxMemoryTag {
    double *restrict Xsqrt;
    double *restrict evalues;
    double *restrict eigenvectors;
    int *restrict evectorsSupport;
    double *restrict dblWorkMem;
    int *restrict intWorkMem;

    double *restrict XsqrtInvX;
    double *restrict Q;
    double *restrict residuals;

    int dblWorkMemSize;
    int intWorkMemSize;
    int nvar;
    int nobs;
} AuxMemory;

void initAuxMemory(AuxMemory* auxmem);
void resizeAuxMemory(AuxMemory* auxmem, const int nvar, const int nobs);
void freeAuxMemory(AuxMemory* auxmem);
void resizeDblWorkAuxMemory(AuxMemory* auxmem, const int newDblWorkMemSize);
void resizeIntWorkAuxMemory(AuxMemory* auxmem, const int newIntWorkMemSize);


#endif /* AuxMemory_h */
