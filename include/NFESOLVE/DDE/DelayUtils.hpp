#ifndef DELAYUTILSHEADERDEF
#define DELAYUTILSHEADERDEF

/*
#include "SparseDelayMatrix.hpp"
#include "Debug.hpp"
*/
#include "../../NFESOLVE.hpp"
#include <armadillo>

void DelaySetup(const double v, const arma::mat& W, const arma::mat& D, arma::vec& delays, arma::umat& ZLocations, arma::uvec& colIndex, arma::uvec& rowPointer);

arma::uword Compute_k(const arma::uword i, const arma::uword j, const arma::uvec& mRowPointer, const arma::uvec& mColIndex);

void ComputeDelayedVars(const arma::vec& u, const arma::mat& mpConnectivity, const SparseDelayMatrix& Z, const arma::uvec& mRowPointer, const arma::uvec& mColIndex, arma::mat& uDelayed, const arma::uword varIndex);

#endif
