#include <math.h>
#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;



int main()
{
    float m, k, x0, v0, t0, h, tmax, x, v, xt, vt, t, energy;

    m = 1.0; // mass
    k = 1.0; // spring constant
    x0 = 1.0; // initial position
    v0 = 0.0; // initial velocity
    t0 = 0.0; // initial time
    tmax = 10.0; // max time
    energy = 0; // energy

    h = 0.001; // time-step

    xt = x0;
    vt = v0;
    t = t0;
    ofstream myfile;
    myfile.open ("SHM_07.txt");

    myfile << "Time Step = " << h << endl;
    myfile << "k = " << k << endl;
    myfile << "m = " << m << endl;
    myfile << "      Time                   x                   v                energy*2" << endl;

    for(t = t0; t <= tmax+h; t+=h){
        v = vt + h*(-k/m)*xt;
        x = xt + h*vt;
        xt = x;
        vt = v;
        energy = 0.5*((pow(v, 2.0)/m) + (k*pow(x, 2.0)));
        myfile << setw(10) << t << setw(20) << xt << setw(20) << vt << setw(20) << energy*2.0 << endl;

    }



}
