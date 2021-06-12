/*
 * Utility functions for constructing the history matrix Z and extracting variables without having to explicitly reference matrix coords.
 */

#include "DelayUtils.hpp"
#include <armadillo>

/*
 * Set up delay matrix from a distance matrix.
 * To be called from the "driver" file of a simulation.
 * Usage:
 * 	Takes conductance velocity v, connectivity matrix W, distance matrix D
 * 	Populates delays vector & CSR coordinates for the SparseDelayMatrix Z
 */
void DelaySetup(const double v, const arma::mat& W, const arma::mat& D, arma::vec& delays, arma::umat& ZLocations, arma::uvec& colIndex, arma::uvec& rowPointer) {

	arma::uword N = W.n_rows;

	int numDelayEqs = 2*N;
	arma::uword numDelays = arma::accu(arma::trimatu(W,1) >= 1e-9);
	arma::uword numDelayValues = 2 * numDelays;

	arma::uword delayIndex = 0;
	arma::uword colCounter = 0;
	for (arma::uword i = 0; i < N; ++i)
	{
		rowPointer(i) = colCounter;
		for (arma::uword j = 0; j < N; ++j)
		{
			if (W(i,j) == 0.0 || D(i,j) < 1e-9)
			{
				continue;
			}
			if (j > i)
			{
				delays(delayIndex) = D(i,j) / v;
				ZLocations(0, 2*delayIndex) = i;
				ZLocations(0, 2*delayIndex + 1) = j;
				ZLocations(1, 2*delayIndex) = delayIndex;
				ZLocations(1, 2*delayIndex + 1) = delayIndex;

				++delayIndex;
				colIndex(colCounter) = j;
				++colCounter;

			}
		}
	}

}
/*
 * Given CSR column & row index vectors, finds which column of SparseDelayMatrix Z corresponds to an (i,j) matrix entry
 */
arma::uword Compute_k(const arma::uword i, const arma::uword j, const arma::uvec& mRowPointer, const arma::uvec& mColIndex)
{
	for (arma::uword k = mRowPointer(i); k < mRowPointer(i+1); ++k)
	{
		if (mColIndex(k) == j)
		{
			return k;
		}
	}
	return 0;
}

/*
 * Extracts a matrix of delayed variables x^{k}_i(t - tau_{ij}).
 * x^{k} is the k-th state variable, i is the node index.
 * Usage:
 * 	Takes vector of variable of interest, connectivity matrix, variable index k
 * 	populates the NxN matrix uDelayed
 * 	assumes that all tau_{ii}s are zero
 * 	(TBD: constant offset delay?)
 */
void ComputeDelayedVars(const arma::vec& uVar, const arma::mat& mpConnectivity, const SparseDelayMatrix& Z, const arma::uvec& mRowPointer, const arma::uvec& mColIndex, arma::mat& uDelayed, const arma::uword varIndex)
{
	arma::uword n = mpConnectivity.n_rows;
	DebugCheck(uVar.n_elem != n, "uVar vector dimensions mismatch with connectivity matrix");
	DebugCheck( ( (uDelayed.n_rows != n ) || (uDelayed.n_cols != n) ), "Delayed variable matrix dimensions mismatch");


		for (arma::uword i = 0; i < n; ++i)
			{
				for (arma::uword j = 0; j < n; ++j)
				{
					if ((mpConnectivity)(i,j) == 0.0)
					{
						continue;
					}
					if (i == j)
					{
						uDelayed(i,j) = uVar(i);
					}
					else if (i > j)
					{
						arma::uword k = Compute_k(j, i, mRowPointer, mColIndex);
						uDelayed(i,j) = Z(i + n * varIndex, k);
					}
					else
					{
						arma::uword k = Compute_k(i,j, mRowPointer, mColIndex);
						uDelayed(i,j) = Z(i + n * varIndex, k);
					}
				}
			}

}
