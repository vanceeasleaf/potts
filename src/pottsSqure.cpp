/*
* @Author: YangZhou
* @Date:   2017-06-13 20:13:20
* @Last Modified by:   YangZhou
* @Last Modified time: 2017-06-19 22:17:54
*/
#include "pottsSqure.h"
#include <math.h>
pottsSqure::pottsSqure() : Potts() {
  // ctor
}

pottsSqure::~pottsSqure() {
  // dtor
}
void pottsSqure::i2coord(int i, int &xi, int &yi, int &ni) {
  xi = i % Lx;
  yi = floor((double)i / Lx);
}
int pottsSqure::coord2i(int x, int y, int n) {
  int i = y * Lx + x;
  return i;
}
void pottsSqure::nei(int xi, int yi, int ni, int index, int &x, int &y,
                     int &n) {
  switch (index) {
    case 0:
      x = xi + 1;
      y = yi;
      break;
    case 1:
      x = xi;
      y = yi + 1;
      break;
    case 2:

      x = xi - 1;
      y = yi;
      break;
    case 3:

      x = xi;
      y = yi - 1;
      break;
  }
  if (x == Lx) x = 0;
  if (x == -1) x = Lx - 1;
  if (y == Ly) y = 0;
  if (y == -1) y = Ly - 1;
}

void pottsSqure::constructNeighbor() {
  nNei = 4;
  Natom = Lx * Ly;
  x = new double[Natom];
  y = new double[Natom];
  neighbor = new int *[Natom];
  for (int i = 0; i < Natom; i++) {
    neighbor[i] = new int[nNei];
    x[i] = y[i] = 0;
  }
  // fstream fneibor("neibor.txt",ios::out);
  int xx, yy, n;
  for (int i = 0; i < Natom; i++) {
    int xi, yi, ni;
    i2coord(i, xi, yi, ni);
    x[i] = xi;
    y[i] = yi;

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
