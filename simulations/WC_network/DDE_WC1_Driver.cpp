//// DDE_Example1 solves a DDE with a periodic solution using a fixed step third order Runge-Kutta solver ////

// Include relevant headers
#include "DDE_WC1.hpp"
#include "NFESOLVE.hpp"

#include <armadillo>
#include <cmath>
#include <iostream>
#include <string>


int main(int argc, char* argv[])
{
	if (argc <= 2)
	{
		printf("Insufficient input matrices!");
		return 0;
	}
	else if (argc > 3)
	{
		printf("Too many input matrices!");
		return 0;
	}
	else if (argc == 3) {
	// Set up solver parameters
	double initialTime = 0.0;
	double finalTime = 5;
	double ATol = 1e-6;
	double RTol = 1e-3;

	std::string outputFileName = "DDE_WC1.dat";
	int saveGap = 10;
	int printGap = 100;

	// Set up DDE parameters
	double tauE = 0.01;
	double tauI = 0.02;
	double PE= 0.34;
	double cee = 3.5;
	double cei = 3.75;
	double cie = -2.5;
	double C = 0.1;
	//sigmoid parameters
	double mu = 1.0;
	double sigm = 0.25;
	arma::vec P = {tauE, tauI, PE, cee, cei, cie, C, mu, sigm};

	double v = 4.0; //conductance velocity

	arma::mat W;
	W.load(argv[1]);
	arma::uword N = W.n_rows;

	arma::mat D;
	D.load(argv[2]);
	arma::uword ND = D.n_rows;
	
	if (N != ND) 
	{ 
		printf("Connectivity and delay matrices are different sizes!");
	}


	arma::uword numDelays = arma::accu(arma::trimatu(W,1) >= 1e-9);
	arma::uword numDelayValues = 2 * numDelays;
	arma::umat ZLocations(2, numDelayValues);
	arma::vec delays(numDelays);
	arma::uvec colIndex(numDelays);
	arma::uvec rowPointer(N);

	DelaySetup(v, W, D, delays, ZLocations, colIndex, rowPointer);

	// Set up initial condition
	arma::vec IC(2 * N, arma::fill::ones);
	//IC *= 0.1;
	IC.subvec(0,N-1) *= 0.5;
	IC.subvec(N,2*N-1) *= 0.1;

	// Set up DDE
	DDE_WC prob(P, W, IC.subvec(0, N - 1), colIndex, rowPointer);

	// Set up solver
	DelaySparseRungeKutta32Solver solver(prob, IC, delays, N, ZLocations, initialTime, finalTime, ATol, RTol, outputFileName, saveGap, printGap);

	// Solve
	solver.Solve();

	return 1;
	}
}


