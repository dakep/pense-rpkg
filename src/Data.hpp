//
//  Data.hpp
//  pense
//
//  Created by David Kepplinger on 2016-01-24.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef Data_hpp
#define Data_hpp

#include <cstring>
#include "config.h"

class Data {
public:
    Data() : nobs(0), nvar(0)
    {}

    Data(double *Xtr, double * y, int nobs, int nvar) :
        Xtr(Xtr), y(y), nobs(nobs), nvar(nvar)
    {}

    Data(const Data &other) : Xtr(other.Xtr), y(other.y), nobs(other.nobs), nvar(other.nvar)
    {}

    Data &operator=(const Data &other) {
        if (this != &other) {
            this->Xtr = other.Xtr;
            this->y = other.y;
            this->nobs = other.nobs;
            this->nvar = other.nvar;
        }
        return *this;
    }

    int numObs() const
    {
        return this->nobs;
    }

    int numVar() const
    {
        return this->nvar;
    }

    void setNumObs(int numObs)
    {
        this->nobs = numObs;
    }

    void setNumVar(int numVar)
    {
        this->nvar = numVar;
    }

    double * getXtr()
    {
        return this->Xtr;
    }

    double * getY()
    {
        return this->y;
    }


    const double * getXtrConst() const
    {
        return this->Xtr;
    }

    const double * getYConst() const
    {
        return this->y;
    }

    /**
     * This function will not check if the memory is large enough
     * to hold the other's data.
     */
    void overrideMemory(const Data& other)
    {
        this->nobs = other.nobs;
        this->nvar = other.nvar;
        memcpy(this->Xtr, other.Xtr, this->nobs * this->nvar * sizeof(double));
        memcpy(this->y, other.y, this->nobs * sizeof(double));
    }

    /**
     * This function will allocate enough memory to hold the other's data.
     * However, it will NOT free the memory before!
     * Call free before if this is wanted!
     */
    void copy(const Data& other)
    {
        if (other.nobs > 0) {
            this->Xtr = new double[other.nobs * other.nvar];
            this->y = new double[other.nobs];

            this->overrideMemory(other);
        } else {
            this->nobs = 0;
            this->nvar = 0;
            this->Xtr = NULL;
            this->y = NULL;
        }
    }

    /**
     * Resize the memory so to be able to hold as many observations
     * and variables as set via setNumObs and setNumVar
     *
     * This function will NOT free the memory previously occupied
     * It is the caller's responsibility to clean up by calling free()!
     *
     */
    void resize()
    {
        this->Xtr = new double[this->nobs * this->nvar];
        this->y = new double[this->nobs];
    }

    void free()
    {
        if (this->nobs > 0) {
            delete[] this->Xtr;
            delete[] this->y;
        }

        this->Xtr = NULL;
        this->y = NULL;
        this->nobs = 0;
        this->nvar = 0;
    }

private:
    double *RESTRICT Xtr;
    double *RESTRICT y;
    int nobs;
    int nvar;
};


#endif /* Data_hpp */
