//
//  LapackException.hpp
//  pense
//
//  Created by David Kepplinger on 2016-01-28.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef LapackException_h
#define LapackException_h


#include <cstring>
#include <cstdio>
#include <exception>

class LapackException : public std::exception
{
public:
    LapackException(const char* what, const int lapackCode) :
            std::exception(), lapackCode(lapackCode)
    {
        /* Assuming that lapack error code < 999999999 */
        this->whatStr = new char[strlen(what) + 30];
        sprintf(this->whatStr, "%s (Lapack error code %d)", what, this->lapackCode);
    }

    virtual const char* what() const throw()
    {
        return this->whatStr;
    }

    int getLapackErrorCode() const
    {
        return this->lapackCode;
    }

    virtual ~LapackException() throw()
    {
        delete[] this->whatStr;
    }

private:
    const int lapackCode;
    char *whatStr;
};

#endif /* LapackException_h */
