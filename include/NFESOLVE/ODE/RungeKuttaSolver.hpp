#ifndef RUNGEKUTTASOLVERHEADERDEF
#define RUNGEKUTTASOLVERHEADERDEF

#include "AbstractODESolver.hpp"
#include <armadillo>

class RungeKuttaSolver : public AbstractODESolver
{
protected:
	void RungeKuttaAlgorithm(const arma::mat& a, const arma::vec& b, const arma::vec& c);
    void AdaptiveRungeKuttaAlgorithm(const arma::mat& a, const arma::vec& b, const arma::vec& c, const arma::vec& bHat, const double ATol, const double RTol);
    void SetStepSizeMin(const double stepSizeMin);
    void SetFacMin(const double facMin);
    void SetFac(const double fac);
    void SetFacMax(const double facMax);

    // Factors for stepsize adaptation
    double mStepSizeMin = 0.000001;
    double mFacMin = 0.1;
    double mFac = 0.8;
    double mFacMax = 3.0;
};

#endif
