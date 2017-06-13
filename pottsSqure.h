#ifndef POTTSSQURE_H
#define POTTSSQURE_H

#include "potts.h"


class pottsSqure : public Potts
{
    public:
        pottsSqure();
        virtual ~pottsSqure();
    protected:
                virtual void constructNeighbor();
    virtual void i2coord(int i,int &xi,int &yi,int &ni);
    virtual int coord2i(int x,int y,int n);
    virtual void nei(int xi,int yi,int ni,int index,int &x,int &y,int &n);
    private:
};

#endif // POTTSSQURE_H
