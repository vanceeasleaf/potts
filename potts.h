#ifndef POTTS_H
#define POTTS_H
class Potts
{
public:

    double t_start,t_int,t_end;
    int MCS,Lx,Ly,Nq,method,Natom,nstart;
        int nNei;
    int *spin;
    double* energy;
    double *magnetic;
    double* energy2;
    double *magnetic2;
    int** neighbor;
    bool picture;
    bool init;
    Potts();
    virtual ~Potts();
    void run();
protected:
        double* x;
    double* y;
        virtual void constructNeighbor();
    void measure(int mcs);
    void drawImage(int mcs);
    virtual void i2coord(int i,int &xi,int &yi,int &ni);
    virtual int coord2i(int x,int y,int n);
    virtual void nei(int xi,int yi,int ni,int index,int &x,int &y,int &n);
private:

        double T;
        char meText[5][100];
    void mcising(double &en,double &ma,double &cv,double &kai);
    double mean(double *var,int start,int end);
    void SW();
    int** clusters;
    bool** bonds;

    void wolff();
    void metropolis();
    double siteEnergy(int i);

};

#endif // POTTS_H
