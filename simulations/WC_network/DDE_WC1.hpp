#ifndef DDE_WCHEADERDEF
#define DDE_WCHEADERDEF

// Include NFESOLVE library and armadillo
#include "NFESOLVE.hpp"
#include <armadillo>

// Class name is DDE_Example1 and it inherits from DDEInterface
class DDE_WC : public SparseDDEInterface
{
public:
	DDE_WC(const arma::vec& parameters, arma::mat& connectivity, const arma::vec& constantHistory, const arma::uvec& colIndex, const arma::uvec& rowPointer);

	// ComputeF method deriving from DDEInterface
	void ComputeF(const double t, const arma::vec& u, const SparseDelayMatrix& Z, arma::vec& F) const;

	// ComputeHistory method deriving from DDEInterface
	void ComputeHistory(const double t, arma::vec& history) const;

	void SetConnectivityKernel(arma::mat& connectivity);

	void Compute_uE_Delay(const arma::vec& uE, const SparseDelayMatrix& Z, arma::mat& uEdelays) const;

private:
	// Private storage for parameters vector
	arma::vec mParameters;
	arma::mat* mpConnectivity;
	arma::vec mHistoryVal;
	//arma::mat* mpDelays;
	arma::uvec mColIndex;
	arma::uvec mRowPointer;

	arma::uword Compute_k(const arma::uword i, const arma::uword j) const;
};

#endif
