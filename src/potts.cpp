/*
* @Author: YangZhou
* @Date:   2017-06-13 20:13:20
* @Last Modified by:   YangZhou
* @Last Modified time: 2017-06-19 22:17:38
*/
#include "potts.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fstream>
#include <iostream>
using namespace std;
enum { METROPOLIS, SW0, WOLFF0 };
double rand01() {
  //    return (double)rand()/(double)(RAND_MAX+1);
  int i;
  double x = 0;
  for (i = 1; i < 6; i++) x = (double)(rand() % 10) * .1 + x * .1;
  return x;
}

Potts::Potts() {
  srand(unsigned(time(NULL)));
  MCS = 2000;
  nstart = 1000;
  picture = false;
  method = METROPOLIS;
  Lx = 32;
  Ly = 32;
  init = 0;
  Nq = 4;

  strcpy(meText[METROPOLIS], "METROPOLIS");
  strcpy(meText[SW0], "SW");
  strcpy(meText[WOLFF0], "WOLFF");
}
Potts::~Potts() {
  delete[] x;
  delete[] y;
}

void Potts::run() {
  constructNeighbor();
  fstream fout("result.txt", ios::out);
  fout << Lx << 'x' << Ly << " hex lattic " << meText[method] << " method "
       << endl;
  cout << Lx << 'x' << Ly << " hex lattic " << meText[method] << " method "
       << endl;
  fout << "Temperature\tEnergy\tHeat Capacity\t Magnetic\tPhi" << endl;
  cout << "Temperature\tEnergy\tHeat Capacity\t Magnetic\tPhi" << endl;
  spin = new int[Natom];

  energy = new double[MCS];
  magnetic = new double[MCS];
  energy2 = new double[MCS];
  magnetic2 = new double[MCS];
  clusters = new int *[Natom];
  for (int i = 0; i < Natom; i++) {
    clusters[i] = new int[Natom];
    for (int j = 0; j < Natom; j++) clusters[i][j] = 0;
  }
  bonds = new bool *[Natom];
  for (int i = 0; i < Natom; i++) {
    bonds[i] = new bool[Natom];
    for (int j = 0; j < Natom; j++) bonds[i][j] = false;
  }

  int t_n = (t_end - t_start) / t_int + 1;
  // double ens[t_n];
  // double mas[t_n];
  double en, ma;
  double cv, kai;
  for (int i = 0; i < t_n; i++) {
    T = i * t_int + t_start;
    mcising(en, ma, cv, kai);
    // ens[i]=en;mas[i]=ma;
    fout << T << '\t' << en << '\t' << cv << '\t' << ma << '\t' << kai << endl;
    cout << T << '\t' << en << '\t' << cv << '\t' << ma << '\t' << kai << endl;
    // emit output(T,en,cv,ma,kai);
    // qDebug()<<T<<'\t'<<en;
  }
  fout.close();
  for (int i = 0; i < Natom; i++) {
    delete[] bonds[i];
  }
  delete[] bonds;
  for (int i = 0; i < Natom; i++) {
    delete[] clusters[i];
  }
  delete[] clusters;
  for (int i = 0; i < Natom; i++) delete[] neighbor[i];
  delete[] neighbor;
  delete[] spin;
  delete[] energy;
  delete[] magnetic;
  delete[] energy2;
  delete[] magnetic2;
}

void Potts::mcising(double &en, double &ma, double &cv, double &kai) {
  if (init) {
    for (int i = 0; i < Natom; i++) {
      spin[i] = rand() % Nq;  // initial configuration
    }
  } else {
    int q = rand() % Nq;
    for (int i = 0; i < Natom; i++) {
      spin[i] = q;  // initial configuration
    }
  }
  for (int mcs = 0; mcs < MCS; mcs++) {
    switch (method) {
      case METROPOLIS:
        metropolis();
        break;
      case SW0:
        SW();
        break;
      case WOLFF0:
        wolff();
        break;
    }

    if (mcs >= nstart) {
      measure(mcs);
    }
  }
  en = mean(energy, nstart, MCS - 1);  // average
  ma = mean(magnetic, nstart, MCS - 1);
  cv = (mean(energy2, nstart, MCS - 1) - en * en) / T / T;
  kai = (mean(magnetic2, nstart, MCS - 1) - ma * ma) / T;
}
double Potts::mean(double *var, int start, int end) {
  double s = 0;
  for (int i = start; i <= end; i++) {
    s += var[i];
  }
  return s / (double)(MCS - start + 1);
}

void Potts::SW() {
  double beta = 1 / (T);
  bool occupied[Natom];
  int nclusters[Natom];
  memset(occupied, false, Natom * sizeof(bool));
  memset(nclusters, 0, Natom * sizeof(int));

  for (int i = 0; i < Natom; i++) {
    memset(clusters[i], 0, Natom * sizeof(int));
    memset(bonds[i], false, Natom * sizeof(bool));
  }

  for (int i = 0; i < Natom; i++)
    for (int j = i + 1; j < Natom; j++) {
      bonds[i][j] = rand01() < 1 - exp(-2 * beta);
    }
  // static int b=0;
  // qDebug()<<b;
  for (int swp = 0; swp < Natom; swp++) {
    int i = swp;  // choose a site
    if (occupied[i]) continue;
    occupied[i] = true;
    int cluster[Natom];  // record the current cluster
    memset(cluster, 0, Natom * sizeof(int));
    int ncluster =
        0;  // how many sites in the current cluster,1 represent the site i
    cluster[ncluster++] = i;  // used to flip the cluster
    int outmost[Natom];       // the outmost sites
    int noutmost = 0;
    memset(outmost, 0, Natom * sizeof(int));

    outmost[noutmost++] = i;  // stack push
    while (noutmost > 0)      // construst a graph,broad first search
    {
      int ii = outmost[noutmost - 1];
      noutmost--;  // stack pop
      for (int inei = 0; inei < nNei; inei++) {
        int jj = neighbor[ii][inei];
        if (spin[jj] != spin[ii] || occupied[jj]) continue;
        if (!bonds[i][jj] && !bonds[jj][i]) continue;
        outmost[noutmost++] = jj;
        cluster[ncluster++] = jj;
        occupied[jj] = true;
      }
    }

    for (int j = 0; j < ncluster; j++) clusters[i][j] = cluster[j];
    nclusters[i] = ncluster;
  }
  // qDebug()<<b++;
  for (int i = 0; i < Natom; i++) {
    int q = floor(rand01() * Nq);
    for (int j = 0; j < nclusters[i]; j++) {
      spin[clusters[i][j]] = q;
    }
  }
}
void Potts::wolff() {
  // static int b=0;
  // qDebug()<<b++;
  double beta = 1 / (T);
  for (int swp = 0; swp < Natom; swp++) {
    int i = floor(rand01() * Natom);  // choose a site
    int q = floor(rand01() * Nq);
    int cluster[Natom];  // record the current cluster
    int ncluster =
        1;  // how many sites in the current cluster,1 represent the site i
    cluster[ncluster - 1] = i;  // used to flip the cluster
    bool labeled[Natom];        // visited in the search
    for (int kk = 0; kk < Natom; kk++) labeled[kk] = false;
    labeled[i] = true;
    int outmost[Natom];  // the outmost sites
    for (int kk = 0; kk < Natom; kk++) outmost[kk] = 0;
    int noutmost = 0;
    outmost[noutmost++] = i;  // stack push
    while (noutmost > 0)      // construst a graph,broad first search
    {
      int ii = outmost[noutmost - 1];
      noutmost--;  // stack pop
      for (int inei = 0; inei < nNei; inei++) {
        int jj = neighbor[ii][inei];
        if (spin[jj] != spin[ii] || labeled[jj]) continue;
        double r = rand01();
        if (r < 1 - exp(-2 * beta)) {
          outmost[noutmost++] = jj;
          cluster[ncluster++] = jj;
        }
        labeled[jj] = true;
      }
    }
    for (int j = 0; j < ncluster; j++) spin[cluster[j]] = q;
  }
}

void Potts::metropolis() {
  double beta = 1.0 / (T);
  for (int swp = 0; swp < Natom; swp++) {
    // while(1){
    int i = rand() % Natom;  // choose a site
    int q = rand() % Nq;
    if (spin[i] != q) {
      double e1 = siteEnergy(i);
      int store = spin[i];
      spin[i] = q;
      double e2 = siteEnergy(i);
      spin[i] = store;
      double dH = e2 - e1;
      double w = 1.0;
      // if(dH>0)
      w = exp(-beta * dH);
      double zeta = rand01();
      if (zeta < w) spin[i] = q;
    }
    //   if(i==0)break;
  }
  //}
}
double Potts::siteEnergy(int i) {
  double e = 0.0;
  for (int inei = 0; inei < nNei; inei++) {
    int jj = neighbor[i][inei];
    if (spin[jj] == spin[i]) e--;
  }
  return e;
}
void Potts::constructNeighbor() {
  nNei = 3;
  Natom = Lx * Ly * 2;
  neighbor = new int *[Natom];
  x = new double[Natom];
  y = new double[Natom];
  for (int i = 0; i < Natom; i++) {
    neighbor[i] = new int[nNei];
    x[i] = y[i] = 0;
  }
  // fstream fneibor("neibor.txt",ios::out);
  int xx, yy, n;
  for (int i = 0; i < Natom; i++) {
    int xi, yi, ni;
    i2coord(i, xi, yi, ni);
    x[i] = xi + yi * .5;
    x[i] = x[i] > Lx - 0.5 ? x[i] - Lx : x[i];
    y[i] = yi * sqrt(3) / 2 - ni / sqrt(3);

    // fneibor<<i;
    for (int j = 0; j < nNei; j++) {
      nei(xi, yi, ni, j, xx, yy, n);
      neighbor[i][j] = coord2i(xx, yy, n);
      // fneibor<<'\t'<<neighbor[i][j];
    }
  }
  /*   for(int i=0;i<Natom;i++){
                   fneibor<<x[i]<<'\t'<<y[i];
                   for(int j=0;j<3;j++){
                           fneibor<<'\t'<<x[neighbor[i][j]]<<'\t'<<y[neighbor[i][j]];
                   }
                   fneibor<<endl;
     }
  fneibor.close();
  */
}
void Potts::measure(int mcs) {
  double m1 = 0;
  double e1 = 0.0;
  for (int i = 0; i < Natom; i++) e1 += siteEnergy(i);
  e1 /= Natom * 2.0;
  int Ns[Nq];
  for (int i = 0; i < Nq; i++) Ns[i] = 0;
  for (int i = 0; i < Natom; i++) {
    Ns[spin[i]]++;
  }
  int nmax = 0;
  for (int i = 0; i < Nq; i++) {
    nmax = nmax > Ns[i] ? nmax : Ns[i];
  }
  m1 = ((double)Nq * (double)nmax / (double)Natom - 1.0) / (double)(Nq - 1);
  energy[mcs] = e1;
  energy2[mcs] = e1 * e1;
  magnetic[mcs] = m1;
  magnetic2[mcs] = m1 * m1;
  if (picture) {
    /*        subplot(1,2,1);
    plot(1:mcs,energy(1:mcs));xlabel('mcs');ylabel('energy per site');
    text = sprintf('mcs = //d',mcs);
    plot_title = sprintf('//dx//d Potts model with //d
    spins\nTemperature=//g,//s',Lx,Ly,Nq,T,text );

    title(plot_title);
    drawnow;
    if(mcs%10==0){
            subplot(1,2,2);
            drawImage(spin,mcs);
    }
    */
  }
}
void Potts::drawImage(int mcs) {
  /* scatter(x,y,30,spin,'filled');
   text = sprintf('mcs = //d',mcs);
   plot_title = sprintf('//dx//d Potts model with //d
  spins\nTemperature=//g,//s',Lx,Ly,Nq,T,text );
   set(gcf,'Color','w');
   title(plot_title);

   //     xlabel(text,'Color','k');
   //     set(gca,'YTickLabel',[],'XTickLabel',[],'XTick',[],'YTick',[]);
   axis off
                   axis equal;
   drawnow;
  }
  */
}
void Potts::i2coord(int i, int &xi, int &yi, int &ni) {
  ni = floor((double)i / Lx / Ly);
  int a = i % (Lx * Ly);
  xi = a % Lx;
  yi = floor((double)a / Lx);
}
int Potts::coord2i(int x, int y, int n) {
  int i = y * Lx + x + n * Lx * Ly;
  return i;
}
void Potts::nei(int xi, int yi, int ni, int index, int &x, int &y, int &n) {
  n = 1 - ni;
  switch (index) {
    case 0:
      x = xi;
      y = yi + 1 - ni;
      break;
    case 1:
      x = xi - 1 + ni;
      y = yi + 1 - 2 * ni;
      break;
    case 2:

      x = xi + ni;
      y = yi - ni;
      break;
  }
  if (x == Lx) x = 0;
  if (x == -1) x = Lx - 1;
  if (y == Ly) y = 0;
  if (y == -1) y = Ly - 1;
}
