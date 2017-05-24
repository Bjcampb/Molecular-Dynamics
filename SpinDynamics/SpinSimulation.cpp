#include <cmath>
#include <math.h>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <iomanip>
#include <ctime>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <random>

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//           Define Constants for System
////////////////////////////////////////////////////////////////////////////////
const int N = 10;     // Number of spins

double randX, randY, randZ; // random values for setting initial positions
double mag;
double Sxo[N], Syo[N], Szo[N]; // t-1 spin components
double Sx[N], Sy[N], Sz[N]; // current spin comonents

////////////////////////////////////////////////////////////////////////////////
//           Initialize Functions
////////////////////////////////////////////////////////////////////////////////
double gasdev();     // Randomly generates numbers from a Gaussian Dist [-1,1]
void createLattice(); // Function that creates lattice
void interaction(); // Calculate exchange interaction
void integrate(double dt);   // Integrate equations of motion

////////////////////////////////////////////////////////////////////////////////
//           Write Functions
////////////////////////////////////////////////////////////////////////////////

double gasdev(){
  static bool available = false;
  static double gset;
  double fac, rsq, v1, v2;
  if (!available)
  {
    do {
      v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
      rsq = v1 * v1 + v2 * v2;
    } while (rsq >= 1.0 || rsq == 0.0);
    fac = sqrt(-2.0 * log(rsq) / rsq);
    gset = v1 * fac;
    available = true;
    return v2*fac;
  }
  else {
    available = false;
    return gset;
  }
}

void createLattice(){
  for (int i = 0; i < N; i++)
  {
    do
    {
      randX = gasdev();
      randY = gasdev();
      randZ = gasdev();
      mag = sqrt(randX * randX + randY * randY + randZ * randZ);
      //cout << mag << setw(5) << i << endl;
    }while(mag > 1);

    Sxo[i] = randX / mag;
    Syo[i] = randY / mag;
    Szo[i] = randZ / mag;
  }

}

void integrate(double dt){
    for(int i = 0; i < N; i++)
    {
        int a = i - 1;
        int b = i + 1;

        // Periodic boundary
        if(a == -1)
        {
            a = N - 1;
        }

        if(b == N)
        {
            b = 0;
        }

        double bx = -(Sxo[a] + Sxo[b]);
        double by = -(Syo[a] + Syo[b]);
        double bz = -(Szo[a] + Szo[b]);

        double dSx = Syo[i] * bz - Szo[i] * by;
        double dSy = Szo[i] * bx - Sxo[i] * bz;
        double dSz = Sxo[i] * by - Syo[i] * bx;

        double sx = Sxo[i] + 2.0 * dSx * dt;
        double sy = Syo[i] + 2.0 * dSy * dt;
        double sz = Szo[i] + 2.0 * dSz * dt;

        double dr = sqrt(sx * sx + sy * sy + sz * sz);

        Sx[i] = sx / dr;
        Sy[i] = sy / dr;
        Sz[i] = sz / dr;

        Sxo[i] = Sx[i];
        Syo[i] = Sy[i];
        Szo[i] = Sz[i];

    }

}

////////////////////////////////////////////////////////////////////////////////
//           Main Function
////////////////////////////////////////////////////////////////////////////////

int main(){


  srand(time(NULL)); // Call srand to reseed random number generator each time
  double dt = 0.01;

  createLattice();


  ofstream initpos;
  initpos.open("InitialPositions.txt");
  initpos << "sxo            syo            szo" << endl;
  for (int i = 0; i < N; i++){
    double length = sqrt(Sxo[i]*Sxo[i] + Syo[i] * Syo[i] + Szo[i] * Szo[i]);
    initpos << Sxo[i] << setw(21) << Syo[i] << setw(21) << Szo[i] << setw(21) << length << endl;
  }

  int NumberSteps = 10000;
  for(int i = 0; i < NumberSteps; i++){
    integrate(dt);
  }

  ofstream finpos;
  finpos.open("FinalPositions.txt");
  finpos << "sx            sy            sz" << endl;
  for (int i = 0; i < N; i++){
    double length = sqrt(Sx[i]*Sx[i] + Sy[i] * Sy[i] + Sz[i] * Sz[i]);
    finpos << Sx[i] << setw(21) << Sy[i] << setw(21) << Sz[i] << setw(21) << length << endl;
  }


}
