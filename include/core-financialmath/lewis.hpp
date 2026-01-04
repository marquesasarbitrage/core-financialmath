#pragma once 
#include <iostream>
#include <complex>
#include "letsberational.hpp"
#include "core-math/quadratures.hpp"

// Model references 
// Option valuation under stochastic volatility - Lewis (2000) : https://econpapers.repec.org/bookchap/vsvvbooks/ovsv.htm
// The volatility surface, a practitioner’s guide - Gatheral (2006) : https://link.springer.com/article/10.1007/s11408-007-0072-4
// The Heston Model and its Extensions in Matlab and C# - Fabrice Douglas Rouah (2013): https://onlinelibrary.wiley.com/doi/book/10.1002/9781118656471

namespace LewisEuropeanVanillaPrice {

    class LewisEuropeanVanillaPrice {

        public: 
            LewisEuropeanVanillaPrice(const GaussLaguerreQuadrature& gaussLaguerreQuadrature); 
            virtual ~LewisEuropeanVanillaPrice() = default; 

            double getNormalizedPrice(double x, bool isCall); 
            double _getNormalizedPrice2(double x, bool isCall);
            double getPrice(double F, double K,  bool isCall); 
            void setQuadraturePoints(int value);

        protected: 
            virtual std::complex<double> _getCharacteriticFunction(std::complex<double> u) const = 0;
        
        private: 
            GaussLaguerreQuadrature gaussLaguerreQuadrature_;
            double _getNormalizedPrice(double x, bool isCall);
            static constexpr double LAGUERRE_WEIGHT_MIN = 1e-20;
            static constexpr double INTEGRAND_REAL_PART_MIN = 1e-20;
        
    };

    class BlackScholes final: public LewisEuropeanVanillaPrice {

        public: 
            BlackScholes(double sigma, double timeToMaturity, const GaussLaguerreQuadrature& gaussLaguerreQuadrature); 
            BlackScholes(double sigma, double timeToMaturity, int n);
            BlackScholes(double sigma, double timeToMaturity);
            ~BlackScholes() = default; 

            double getTimeToMaturity() const; 
            double getSigma() const; 
            void setTimeToMaturity(double value); 
            void setSigma(double value); 

        
        protected: 
            std::complex<double> _getCharacteriticFunction(std::complex<double> u) const override final;
        
        private: 
            double sigma_; 
            double timeToMaturity_; 
            
    };

    class Heston final : public LewisEuropeanVanillaPrice {

        public: 
            Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity, const GaussLaguerreQuadrature& gaussLaguerreQuadrature); 
            Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity, int n);
            Heston(double meanReversion, double longTermVariance, double varianceOfVariance, double correlation, double initialVariance, double timeToMaturity);
            ~Heston() = default; 
        
            double getTimeToMaturity() const; 
            double getMeanReversion() const; 
            double getInitialVariance() const; 
            double getLongTermVariance() const; 
            double getCorrelation() const; 
            double getVarianceOfVariance() const; 

            void setMeanReversion(double value); 
            void setInitialVariance(double value); 
            void setLongTermVariance(double value); 
            void setCorrelation(double value); 
            void setVarianceOfVariance(double value); 
            void setTimeToMaturity(double value); 

            bool isFellerConditionSatisfied() const;

        
        protected: 
            std::complex<double> _getCharacteriticFunction(std::complex<double> u) const override final;
        
        private: 
            double kappa_;
            double theta_; 
            double eta_;
            double rho_; 
            double v0_;
            double T_; 

    };

}




