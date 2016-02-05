//
//  IndexNumber.hpp
//  penseinit
//
//  Created by David Kepplinger on 2016-01-26.
//  Copyright © 2016 David Kepplinger. All rights reserved.
//

#ifndef IndexNumber_hpp
#define IndexNumber_hpp

#include <float.h>

class IndexNumber
{
public:
    enum SortOrder {
        INCREASING = 0,
        DECREASING = 1,
        ABSOLUTE_INCREASING = 2
    };

    IndexNumber(const double *const values, const int position = 0, const SortOrder sortOrder = INCREASING) :
        values(values), sortOrder(sortOrder), position(position)
    {}

    IndexNumber(const IndexNumber& other) : values(other.values), sortOrder(other.sortOrder), position(other.position)
    {}

    IndexNumber& operator=(IndexNumber& other)
    {
        this->position = other.position;
        return *this;
    }

    void setSortOrder(const SortOrder sortOrder)
    {
        this->sortOrder = sortOrder;
    }

    void setPosition(int pos)
    {
        this->position = pos;
    }

    int getPosition() const
    {
        return this->position;
    }

    double getValue() const
    {
        return this->values[this->position];
    }

    const double* getValues() const
    {
        return this->values;
    }

    bool operator<(const double b) const
    {
        switch (this->sortOrder) {
            case DECREASING:
                return b < this->values[this->position];
            case ABSOLUTE_INCREASING:
                return fabs(b) > fabs(this->values[this->position]);
            case INCREASING:
            default:
                return b > this->values[this->position];
        }
    }

    bool operator<=(const double b) const
    {
        return (*this < b) || (b == this->values[this->position]);
    }

    bool operator<(const IndexNumber &b) const
    {
        switch (this->sortOrder) {
            case DECREASING:
                return b.values[b.position] < this->values[this->position];
            case ABSOLUTE_INCREASING:
                return fabs(b.values[b.position]) > fabs(this->values[this->position]);
            case INCREASING:
            default:
                return b.values[b.position] > this->values[this->position];
        }

    }

private:
    const double *const values;
    SortOrder sortOrder;
    int position;
};


#endif /* IndexNumber_hpp */
