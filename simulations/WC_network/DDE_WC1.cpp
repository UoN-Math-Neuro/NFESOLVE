// Include DDE_Example1 header
#include "DDE_WC1.hpp"

// Constructor
DDE_WC::DDE_WC(const arma::vec& parameters, arma::mat& connectivity, const arma::vec& constantHistory, const arma::uvec& colIndex, const arma::uvec& rowPointer)
{
	// Set private parameters vector and amplitudes
	mParameters = parameters;
	mpConnectivity = &connectivity;
	mHistoryVal = constantHistory;
	//mpDelays = &delayMatrix;
	mColIndex = colIndex;
	mRowPointer = rowPointer;

}
//Sigmoid function
arma::vec sigmoid(arma::vec x, const double mu, const double sigm) {
	arma::vec y;
	y.copy_size(x);

	y = pow((1.0 + exp(-(x - mu) / sigm)), -1.0);
	return y;
}

// ComputeF
void DDE_WC::ComputeF(const double t, const arma::vec& u, const SparseDelayMatrix& Z, arma::vec& F) const
{
	arma::uword n = mpConnectivity->n_rows; //number of nodes
	//Split vector into variables
	arma::vec uE = u.subvec(0,(n-1));
	arma::vec uI = u.subvec(n,(2*n-1));

	arma::arma_rng::set_seed_random();
	arma::vec xiE(n);
	arma::vec xiI(n);

	xiE.randn();
	xiI.randn();

	double tauE = mParameters(0);
	double tauI = mParameters(1);
	double PE = mParameters(2);
	double cee = mParameters(3);
	double cei = mParameters(4);
	double cie = mParameters(5);
	double C = mParameters(6);
	double mu = mParameters(7);
	double sigm = mParameters(8);

	arma::mat uEdelays(n, n, arma::fill::zeros);
	ComputeDelayedVars(uE, *mpConnectivity, Z, mRowPointer, mColIndex, uEdelays, 0);

	//arma::vec DuE = (-uE + sigmoid((cee * uE + cie * uI + PE + 0.01 * xiE + C * *mpConnectivity * uE), mu, sigm)) / tauE;
	//arma::vec DuE = (-uE + sigmoid((cee * uE + cie * uI + PE + 0.1 * xiE + C * *mpConnectivity * Z.col(0).subvec(0,n-1)), mu, sigm)) / tauE;
	arma::vec DuE = (-uE + sigmoid((cee * uE + cie * uI + PE + 0.1 * xiE + C * sum(*mpConnectivity % uEdelays, 1)), mu, sigm)) / tauE;
	arma::vec DuI = (-uI + sigmoid((cei * uE + 0.01 * xiI), mu, sigm)) / tauI;
	F = join_cols(DuE, DuI);

}


// ComputeHistory
void DDE_WC::ComputeHistory(const double t, arma::vec& history) const
{
	history = mHistoryVal;
}

void DDE_WC::SetConnectivityKernel(arma::mat& connectivity)
{
	mpConnectivity = &connectivity;
}

