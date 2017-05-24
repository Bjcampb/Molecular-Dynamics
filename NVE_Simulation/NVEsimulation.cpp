#include <cmath>
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

//////////////////////////////////////////////////////////
//     Define Constants for System                    ////
//////////////////////////////////////////////////////////

const int N = 256;  // Defines the number of particles in box
double L;           // Defines one dimension of the box
double rho = 0.636;  // Defines the density of the atoms
double cutoff = 2.5; // Defines the cutoff values for the potential calculation
double TotEnergy = 101.79; // Defines desired initial total energy

double rxo[N], ryo[N], rzo[N]; // Array for old positions
double rxn[N], ryn[N], rzn[N]; // Array for current positions
double vx[N], vy[N], vz[N];    // Array for velocity of atoms
double ax[N], ay[N], az[N];    // Array for acceleartion of atoms

double ux[N], uy[N], uz[N];   // Array that helps with rescaling
double vxscale[N], vyscale[N], vzscale[N];

double momentum;      // Value for total momentum of system
double potential;    // Value for potential of current system
double kinetic;      // Value for kinetic energy of current system
double energy;       // Value for total energy of current system

///////////////////////////////////////////////////////////////
//    Initialize the Functions required for Program         ///
///////////////////////////////////////////////////////////////

void initializePositions();      // Places atoms in a FCC box
double gasdev();                 // Randomly generates numbers from a Gaussian Distribution -1 to 1
void rescaleVel();               // Rescales velocities that are initialized from gasdev to fit total energy
void initializeVelocities();     // Determines velocities for each atom to fit total energy
void calculateAcceleration();
void verletAlgorithm(double dt);
void initializeBox();



void calculateInitialPotential();       // Calculates initial potential of box for current state
void calculateKinetic();         // Calculates kinetic energy of box for current state
void calculateEnergy();          // Calculates total energy


///////////////////////////////////////////////////////////////
//     Write fucntions                                      ///
///////////////////////////////////////////////////////////////


////// Box initializaiton functions ///////////////////
//float RandomFloat() {
//	float random = ((float)rand()) / (float)RAND_MAX;
//
//	float range = L;
//	return (random*range);
//}

void initializePositions() {

	L = pow(N / rho, 1.0 / 3); // Length of one dimension of box
	int M = 1;
	while (4 * M * M * M < N) {
		++M;
	}

	double xFCC[4] = { 0.25, 0.75, 0.75, 0.25 }; // x position spacing for FCC lattice
	double yFCC[4] = { 0.25, 0.75, 0.25, 0.75 }; // y position spacing for FCC lattice
	double zFCC[4] = { 0.25, 0.25, 0.75, 0.75 }; // z position spacing for FCC lattice

	double d = L / M; // Lattice constant
	int n = 0; // placed molecule count
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < M; j++) {
			for (int k = 0; k < M; k++) {
				for (int F = 0; F < 4; F++) {
					if (n < N) {
						rxo[n] = (i + xFCC[F]) * d;
						ryo[n] = (j + yFCC[F]) * d;
						rzo[n] = (k + zFCC[F]) * d;
						++n;
					}
				}
			}
		}
	}


	ofstream initialpos;
	initialpos.open("InitialPositions.txt");
	initialpos << "rxo            ryo            rzo" << endl;
	for (int i = 0; i < N; i++) {
		initialpos << rxo[i] << setw(21) << ryo[i] << setw(21) << rzo[i] << endl;
	}

}

double gasdev() {
	static bool available = false;
	static double gset;
	double fac, rsq, v1, v2;
	if (!available) {
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

void initializeVelocities() {

	for (int i = 0; i < N; i++) {
		vx[i] = gasdev();
		vy[i] = gasdev();
		vz[i] = gasdev();
	}

	double vxCOM = 0;
	double vyCOM = 0;
	double vzCOM = 0;

	for (int i = 0; i < N; i++) {
		vxCOM += vx[i];
		vyCOM += vy[i];
		vzCOM += vz[i];
	}

	vxCOM /= N;
	vyCOM /= N;
	vzCOM /= N;

	for (int i = 0; i < N; i++) {
		vx[i] -= vxCOM;
		vy[i] -= vyCOM;
		vz[i] -= vzCOM;
	}
}

void rescaleVel() {

	double desiredKinetic = TotEnergy - potential;

	calculateKinetic();

	double lambda = sqrt(2.0 * desiredKinetic / kinetic);


	for (int i = 0; i < N; i++) {
		vx[i] *= lambda;
		vy[i] *= lambda;
		vz[i] *= lambda;
	}

	double multiplier = .7071067; // Multiplier to raise values to correct kinetic for total energy
	for (int i = 0; i < N; i++) {
		vx[i] *= multiplier;
		vy[i] *= multiplier;
		vz[i] *= multiplier;
	}


}

void calculateInitialPotential() {
	potential = 0; // Set potential value to zero before calulation

	for (int i = 0; i < (N - 1); i++) {
		for (int j = i + 1; j < N; j++) {
			double RXIJ = rxo[i] - rxo[j];
			double RYIJ = ryo[i] - ryo[j];
			double RZIJ = rzo[i] - rzo[j];

			if (fabs(RXIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RXIJ > 0) {
					RXIJ -= L;
				}
				else RXIJ += L;
			}

			if (fabs(RYIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RYIJ > 0) {
					RYIJ -= L;
				}
				else RYIJ += L;
			}

			if (fabs(RZIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RZIJ > 0) {
					RZIJ -= L;
				}
				else RZIJ += L;
			}


			double dr = sqrt((RXIJ * RXIJ) + (RYIJ * RYIJ) + (RZIJ * RZIJ)); // Calculate distance between atoms

			if (dr < cutoff) {
				double pot = pow(dr, -12) - pow(dr, -6);
				potential += pot;
				//initpot << RXIJ << setw(21) << RYIJ << setw(21) << RZIJ << setw(21) << dr << setw(21) << 4.0 * pot << endl;
			}
		}
	}

	potential = 4.0 * potential;
	//initpot << "Total Potential: " << potential << endl;

}

void initializeBox() {
	initializePositions();
	calculateInitialPotential();
	initializeVelocities();
	rescaleVel();

}


//////// Acceleartion and Verlet Functions ///////////////////
void calculateAcceleration() {

	for (int i = 0; i < N; i++) { // Set acceleration values equal to zero before calculating based on LJ potential
		ax[i] = 0;
		ay[i] = 0;
		az[i] = 0;
	}

	for (int i = 0; i < (N - 1); i++) {
		for (int j = i + 1; j < N; j++) {
			double dx = rxn[i] - rxn[j];
			double dy = ryn[i] - ryn[j];
			double dz = rzn[i] - rzn[j];

			if (fabs(dx) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (dx > 0) {
					dx -= L;
				}
				else dx += L;
			}

			if (fabs(dy) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (dy > 0) {
					dy -= L;
				}
				else dy += L;
			}

			if (fabs(dz) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (dz > 0) {
					dz -= L;
				}
				else dz += L;
			}

			double dr = sqrt((dx*dx) + (dy * dy) + (dz * dz));

			if (dr < cutoff) {
				double f = 48 * (pow(dr, -14) - (0.5*pow(dr, -8)));


				ax[i] += (dx) * f;
				ax[j] -= (dx) * f;

				ay[i] += (dy) * f;
				ay[j] -= (dy) * f;

				az[i] += (dz) * f;
				az[j] -= (dz) * f;
			}


		}
	}

}

void verletAlgorithm(double dt) {
	calculateAcceleration();

	double DTSQ = dt*dt;
	double DT2 = 2.0 * dt;

	for (int i = 0; i < N; i++) {
		double RXNEWI = 2.0 * rxn[i] - rxo[i] + DTSQ * ax[i];
		double RYNEWI = 2.0 * ryn[i] - ryo[i] + DTSQ * ay[i];
		double RZNEWI = 2.0 * rzn[i] - rzo[i] + DTSQ * az[i];
		vx[i] = (RXNEWI - rxo[i]) / DT2;
		vy[i] = (RYNEWI - ryo[i]) / DT2;
		vz[i] = (RZNEWI - rzo[i]) / DT2;

		// Periodic Boundary Conditions
		if (RXNEWI < 0) {
			RXNEWI += L;
		}
		if (RXNEWI >= L) {
			RXNEWI -= L;
		}
		if (RYNEWI < 0) {
			RYNEWI += L;
		}
		if (RYNEWI >= L) {
			RYNEWI -= L;
		}
		if (RZNEWI < 0) {
			RZNEWI += L;
		}
		if (RZNEWI >= L) {
			RZNEWI -= L;
		}

		rxo[i] = rxn[i];
		ryo[i] = ryn[i];
		rzo[i] = rzn[i];

		rxn[i] = RXNEWI;
		ryn[i] = RYNEWI;
		rzn[i] = RZNEWI;
	}

}



/////////// Value Calculations    ////////////////
void calculateMomentum() {
	momentum = 0;
	double VSUM = 0;
	for (int i = 0; i < N; i++) {
		VSUM += vx[i] + vy[i] + vz[i];
	}
	momentum = VSUM;
}

void calculatePotential() {
	potential = 0; // Set potential value to zero before calculation


	for (int i = 0; i < (N - 1); i++) {
		for (int j = i + 1; j < N; j++) {
			double RXIJ = rxn[i] - rxn[j];
			double RYIJ = ryn[i] - ryn[j];
			double RZIJ = rzn[i] - rzn[j];

			if (fabs(RXIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RXIJ > 0) {
					RXIJ -= L;
				}
				else RXIJ += L;
			}

			if (fabs(RYIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RYIJ > 0) {
					RYIJ -= L;
				}
				else RYIJ += L;
			}

			if (fabs(RZIJ) > 0.5 * L) { // If distance is greater than half L then use image from other box
				if (RZIJ > 0) {
					RZIJ -= L;
				}
				else RZIJ += L;
			}


			double dr = sqrt((RXIJ * RXIJ) + (RYIJ * RYIJ) + (RZIJ * RZIJ)); // Calculate distance between atoms

			if (dr < cutoff) {
				double pot = pow(dr, -12) - pow(dr, -6);
				potential += pot;
			}
		}
	}

	potential = 4.0 * potential;

}

void calculateKinetic() {

	kinetic = 0;

	for (int i = 0; i < N; i++) {
		double VSQDSUM = (vx[i] * vx[i]) + (vy[i] * vy[i]) + (vz[i] * vz[i]);
		kinetic += VSQDSUM;
	}

	kinetic = 0.5 * kinetic;
}

void calculateEnergy() {
	energy = 0;

	energy = kinetic + potential;
}

/////////////////////////////////////////////////////////////////
//                 Main Function                               //
/////////////////////////////////////////////////////////////////

int main() {
	srand(time(NULL));//srand(1);//srand(time(NULL)); // Call srand to reseed random number generator each time
	double dt = 0.01;

	initializeBox();
	/////////////////// Output Initial Info
	calculateMomentum();
	calculateInitialPotential();
	calculateKinetic();
	calculateEnergy();
	cout << "Initial Momentum: " << momentum << endl;
	cout << "Initial Potential: " << potential << endl;
	cout << "Initial Kinetic: " << kinetic << endl;
	cout << "Total Energy: " << energy << endl;
	//////////////////////////////////////////

	//////////// Euler Step //////////////////
	for (int i = 0; i < N; i++) {
		rxn[i] = rxo[i] + vx[i] * dt;
		ryn[i] = ryo[i] + vy[i] * dt;
		rzn[i] = rzo[i] + vz[i] * dt;

		if (rxn[i] < 0) {
			rxn[i] += L;
		}
		if (rxn[i] >= L) {
			rxn[i] -= L;
		}
		if (ryn[i] < 0) {
			ryn[i] += L;
		}
		if (ryn[i] >= L) {
			ryn[i] -= L;
		}
		if (rzn[i] < 0) {
			rzn[i] += L;
		}
		if (rzn[i] >= L) {
			rzn[i] -= L;
		}
	}

	ofstream euler;
	euler.open("EulerStepPos.txt");
	euler << "rx                        ry                         rz" << endl;
	for (int i = 0; i < N; i++) {
		euler << rxn[i] << setw(21) << ryn[i] << setw(21) << rzn[i] << endl;
	}

	calculateAcceleration();
	ofstream initAcc;
	initAcc.open("InitialAccelerations.txt");
	initAcc << "ax                        ay                              az" << endl;
	for (int i = 0; i < N; i++) {
		initAcc << ax[i] << setw(21) << ay[i] << setw(21) << az[i] << endl;
	}

	/////////////// Verlet Algorithm //////////////
	int NumberSteps = 12;
	ofstream verlet;
	verlet.open("Energy.txt");
	verlet << "Number of Steps: " << NumberSteps << endl;
	verlet << "Size of Step: " << dt << endl;
	verlet << "Momentum                         Kinetic                 Potential                      Energy" << endl;
	for (int i = 0; i < NumberSteps; i++) {
		verletAlgorithm(dt);
		calculateMomentum();
		calculateKinetic();
		calculatePotential();
		calculateEnergy();
		verlet << momentum << setw(21) << kinetic << setw(21) << potential << setw(21) << energy << endl;
	}

	// Debugging Outputs
	ofstream rxnpos;
	rxnpos.open("RXNpositions12.txt");
	rxnpos << "rx                   ry                rz" << endl;
	for (int i = 0; i < N; i++) {
		rxnpos << rxn[i] << setw(21) << ryn[i] << setw(21) << rzn[i] << endl;
	}

	ofstream vel;
	vel.open("Velocities12.txt");
	vel << "vx                   vy                vz" << endl;
	for (int i = 0; i < N; i++) {
		vel << vx[i] << setw(21) << vy[i] << setw(21) << vz[i] << endl;
	}

	ofstream accel;
	accel.open("Accelerations12.txt");
	accel << "ax                   ay                az" << endl;
	for (int i = 0; i < N; i++) {
		accel << ax[i] << setw(21) << ay[i] << setw(21) << az[i] << endl;
	}


}
