#include "../include/core-financialmath/lewis.hpp"

namespace LewisEuropeanVanillaPrice {

    LewisEuropeanVanillaPrice::LewisEuropeanVanillaPrice(const GaussLaguerreQuadrature& gaussLaguerreQuadrature): 
    gaussLaguerreQuadrature_(gaussLaguerreQuadrature)
    {}

    void LewisEuropeanVanillaPrice::setQuadraturePoints(int value) { gaussLaguerreQuadrature_.setPoints(value); }

    double LewisEuropeanVanillaPrice::_getNormalizedPrice(double x, bool isCall) {

        int n = gaussLaguerreQuadrature_.getPoints();
        double sum = 0.0, u, w;
        std::complex<double> integrand;
        std::complex<double> i(0.0, 1.0);
        for (int j = 0; j < n; ++j) {
            
            u = gaussLaguerreQuadrature_.getRoot(j); w = gaussLaguerreQuadrature_.getWeight(j);
            integrand = std::exp(i * u * x) * _getCharacteriticFunction(std::complex<double>(u, -0.5));
            if (abs(integrand.real()) < INTEGRAND_REAL_PART_MIN or abs(w)< LAGUERRE_WEIGHT_MIN) {break;}
            else {sum += w * std::exp(u) * integrand.real() / (u * u + 0.25);}
            
        }
        double call = std::exp(x / 2.0) - sum / LetsBeRational::PI;
        return isCall ? call : call - LetsBeRational::getNormalizedIntrisicValue(x, true);
    }

    double LewisEuropeanVanillaPrice::_getNormalizedPrice2(double x, bool isCall) {

        int n = gaussLaguerreQuadrature_.getPoints();
        std::vector<double> laguerreNodes = gaussLaguerreQuadrature_.getRoots();
        std::vector<double> laguerreWeights = gaussLaguerreQuadrature_.getWeights();
        std::complex<double> i(0.0, 1.0);
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            std::cout << "U: " << laguerreNodes[j] << std::endl;
            std::cout << "W: " << laguerreWeights[j] << std::endl;
            double u = laguerreNodes[j];
            std::complex<double> integrand = std::exp(i * u * x) * _getCharacteriticFunction(std::complex<double>(u, -0.5));
            double rintegrand = abs(integrand.real()) < DBL_MIN ? 0.0 : integrand.real();
            if (rintegrand == 0.0 or abs(laguerreWeights[j])<1e-20) sum += 0.0;
            else sum += laguerreWeights[j] * std::exp(u) * rintegrand / (u * u + 0.25);
        }
        double call = std::exp(x / 2.0) - sum / LetsBeRational::PI;
        return isCall ? call : call - LetsBeRational::getNormalizedIntrisicValue(x, true);
    }

    double LewisEuropeanVanillaPrice::getNormalizedPrice(double x, bool isCall) {

        try { return _getNormalizedPrice(x, isCall); }
        catch (const std::exception& e) { return NAN; }
    }

    double LewisEuropeanVanillaPrice::getPrice(double F, double K,  bool isCall) {

        return std::sqrt(F*K) * _getNormalizedPrice(std::log(F/K), isCall);
    }

    BlackScholes::BlackScholes(double sigma, double timeToMaturity, const GaussLaguerreQuadrature& gaussLaguerreQuadrature): 
    LewisEuropeanVanillaPrice(gaussLaguerreQuadrature)
    , sigma_(sigma)
    , timeToMaturity_(timeToMaturity) 
    {}

    BlackScholes::BlackScholes(double sigma, double timeToMaturity, int n):
    LewisEuropeanVanillaPrice(GaussLaguerreQuadrature(n))
    , sigma_(sigma)
    , timeToMaturity_(timeToMaturity) {}

    BlackScholes::BlackScholes(double sigma, double timeToMaturity):
    LewisEuropeanVanillaPrice(GaussLaguerreQuadrature(20))
    , sigma_(sigma)
    , timeToMaturity_(timeToMaturity) {}

    double BlackScholes::getTimeToMaturity() const { return timeToMaturity_;}

    double BlackScholes::getSigma() const { return sigma_; }

    void BlackScholes::setTimeToMaturity(double value) { timeToMaturity_ = value; }

    void BlackScholes::setSigma(double value) { sigma_ = value; }

    std::complex<double> BlackScholes::_getCharacteriticFunction(std::complex<double> u) const {

        return std::exp(-0.5*u*(u+std::complex<double>(0.0,1.0))*timeToMaturity_*sigma_*sigma_);
    }

    Heston::Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity, const GaussLaguerreQuadrature& gaussLaguerreQuadrature):
    LewisEuropeanVanillaPrice(gaussLaguerreQuadrature)
    , kappa_(meanReversion)
    , theta_(longTermVariance)
    , eta_(varianceOfVariance)
    , rho_(correlation)
    , v0_(initialVariance)
    , T_(timeToMaturity)
    {}

    Heston::Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity, int n): 
    LewisEuropeanVanillaPrice(GaussLaguerreQuadrature(n))
    , kappa_(meanReversion)
    , theta_(longTermVariance)
    , eta_(varianceOfVariance)
    , rho_(correlation)
    , v0_(initialVariance)
    , T_(timeToMaturity)
    {}

    Heston::Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity): 
    LewisEuropeanVanillaPrice(GaussLaguerreQuadrature(20))
    , kappa_(meanReversion)
    , theta_(longTermVariance)
    , eta_(varianceOfVariance)
    , rho_(correlation)
    , v0_(initialVariance)
    , T_(timeToMaturity)
    {}

    double Heston::getTimeToMaturity() const { return T_; }

    double Heston::getMeanReversion() const { return kappa_; }

    double Heston::getInitialVariance() const { return v0_; }

    double Heston::getLongTermVariance() const { return theta_; }

    double Heston::getCorrelation() const { return rho_; }

    double Heston::getVarianceOfVariance() const { return eta_; }

    void Heston::setMeanReversion(double value) { kappa_ = value; }

    void Heston::setInitialVariance(double value) { v0_ = value; }

    void Heston::setLongTermVariance(double value) { theta_ = value; }

    void Heston::setCorrelation(double value) { rho_ = value; }

    void Heston::setVarianceOfVariance(double value) { eta_ = value; }

    void Heston::setTimeToMaturity(double value) { T_ = value; }

    bool Heston::isFellerConditionSatisfied() const { return (2.0 * kappa_ * theta_) > (eta_ * eta_); }

    std::complex<double> Heston::_getCharacteriticFunction(std::complex<double> u) const {

        std::complex<double> i(0.0, 1.0);
        std::complex<double> beta = kappa_ - rho_*eta_*i*u;
        std::complex<double> d = std::sqrt(beta*beta + eta_*eta_*(i*u + u*u));
        std::complex<double> g = (beta - d)/(beta + d); 
        std::complex<double> D = (beta - d)*(1.0-std::exp(-d*T_))/(eta_*eta_*(1.0-g*std::exp(-d*T_)));
        std::complex<double> C = kappa_*(T_*(beta-d) - 2.0*std::log((1.0-g*std::exp(-d*T_))/(1.0-g)))/(eta_*eta_); 
        return std::exp(C*theta_+D*v0_); 
    }

}