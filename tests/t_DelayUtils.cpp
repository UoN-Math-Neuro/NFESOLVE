#include "NFEsolve.hpp"
#include <armadillo>
#include <assert.h>
#include <string>
#include <iostream>

void testUtils(int testNum)
{
	arma::arma_rng::set_seed_random();
	double v = 1.0;
	arma::mat W, D, uDelayed, uTest;
	arma::vec delays, delayTest, uVar, ZVal;
	arma::umat ZLocations, ZLTest;
	arma::uvec colIndex, rowPointer, colTest, rowTest;
	arma::uword numCols;

	if (testNum == 0)
	{
		W = { {1.0, 1.0}, {1.0, 1.0} };
		D = { {0.0, 2.0}, {2.0, 0.0} };
		delays.set_size(1);
		ZLocations.set_size(2, 2);
		colIndex.set_size(1);
		rowPointer.set_size(2);
		uDelayed.zeros(2, 2);
		uTest.zeros(2, 2);

		delayTest = {2.0};
		ZLTest = { {0, 1}, {0, 0} };
		colTest = {1};
		rowTest = {0, 1};
		uVar.randn(2);
		ZVal.randn(1);
		numCols = 1;
	}
	else if (testNum == 1)
	{
		W = { {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0} };
		D = { {0.0, 1.0, 2.0}, {1.0, 0.0, 3.0}, {2.0, 3.0, 0.0} };
		delays.set_size(2);
		ZLocations.set_size(2, 4);
		colIndex.set_size(2);
		rowPointer.set_size(3);
		uDelayed.zeros(3, 3);
		uTest.zeros(3, 3);

		delayTest = {2.0, 3.0};
		ZLTest = { {0, 2, 1, 2}, {0, 0, 1, 1} };

		colTest = {2, 2};
		rowTest = {0, 1, 2};
		uVar.randn(3);
		ZVal.randn(2);
		numCols = 2;
	}

	DelaySetup(v, W, D, delays, ZLocations, colIndex, rowPointer);

	assert((delays - delayTest).is_zero());

	assert((ZLocations - ZLTest).is_zero());

	assert((colIndex - colTest).is_zero());
	assert((rowPointer - rowTest).is_zero());

	std::cout << "DelaySetup test complete" << std::endl;


	arma::vec ZVals = arma::repelem(ZVal, 2, 1);



	SparseDelayMatrix Z(ZLocations, ZVals, W.n_rows, 2 * numCols);
/*
	arma::uvec ZRowIndex = Z.GetRowIndex();
	arma::uvec ZColPtr = Z.GetColPtr();
	arma::vec ZGetV = Z.GetValues();

	ZGetV.print();
	std::cout << "ColPtr" << std::endl;
	ZColPtr.print();
	std::cout << "RowIndex" << std::endl;
	ZRowIndex.print();


	arma::uword k0 = Compute_k(0, 2, rowPointer, colIndex);
	std::cout << k0 << std::endl;
	//std::cout << Z(0,k0) << std::endl;
	arma::uword k1 = Compute_k(1, 2, rowPointer, colIndex);
	std::cout << k1 << std::endl;
	//std::cout << Z(0,k1) <<  std::endl;

	try{
	std::cout << Z(0,1) <<  std::endl;
	} catch(...){}
	try{
		std::cout << Z(1,1) << std::endl;
	}catch(...){}
	try{
			std::cout << Z(1,2) << std::endl;
		}catch(...){}

	try{
			std::cout << Z(1,2) << std::endl;
		}catch(...){}
	try{
				std::cout << Z(0,2) << std::endl;
			}catch(...){}
*/
/*
 *
 *	std::cout << Z(1,0) << std::endl;
 *	std::cout << Z(0,1) << std::endl;
 *	std::cout << Z(1,1) << std::endl;
 */

	ComputeDelayedVars(uVar, W, Z, rowPointer, colIndex, uDelayed, 0);
//	std::cout << "test1" << std::endl;
//	uDelayed.print();

	uTest.diag() = uVar;
	if (testNum == 0)
	{
		uTest(0, 1) = ZVal(0);
		uTest(1, 0) = ZVal(0);
	}
	else if (testNum == 1)
	{
		uTest(2, 0) = ZVal(0);
		uTest(0, 2) = ZVal(0);
		uTest(2, 1) = ZVal(1);
		uTest(1, 2) = ZVal(1);
 	}

	assert(arma::approx_equal(uDelayed, uTest, "absdiff", 1e-3));

	std::cout << "ComputeDelayedVars test complete" << std::endl;
}


int main()
{
	testUtils(0);
	testUtils(1);
	return 1;
}



/*
void ComputeDelayedVars(const arma::vec& uVar, const arma::mat& mpConnectivity, const SparseDelayMatrix& Z, const arma::uvec& mRowPointer, const arma::uvec& mColIndex, arma::mat& uDelayed, const arma::uword varIndex)
{
	arma::uword n = mpConnectivity.n_rows;

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
*/
