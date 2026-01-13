#pragma once 
#include <iostream>
#include <tuple>
#include "letsberational.hpp"

// Model references 
// Arbitrage-free SVI volatility surfaces - Jim Gatheral and Antoine Jacquier (2013) : https://doi.org/10.48550/arXiv.1204.0646

class SVI {

    public: 
        SVI(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity); 
        SVI(std::tuple<double, double, double, double, double> params, double timeToMaturity); 
        ~SVI() = default; 
        static SVI getSVIFromRawParameters(double a, double b, double m, double s, double p, double timeToMaturity);
        std::tuple<double, double, double, double, double> getRawParameters() const;

        double getImpliedVarianceATM() const; 
        double getImpliedVolatilityATM() const; 
        double getTotalVarianceATM() const; 
        double getImpliedVariance(double x) const; 
        double getImpliedVolatility(double x) const; 
        double getTotalVariance(double x) const; 
        double getMinimumVariance() const; 
        double getMinimumTotalVariance() const; 
        double getMinimumImpliedVolatility() const; 
        double getRiskNeutralDensity(double x) const;
        double getLocalVariance(double x) const;
        double getLocalVolatility(double x) const;    
        double getNormalizedBlackPrice(double x, bool isCall) const;
        double getSkewATM() const; 
        double getCallSlope() const; 
        double getPutSlope() const; 
        double getTimeToMaturity() const;
        bool isFreeOfButterflyArbitrage() const;
        void checkFreeOfButterflyArbitrage() const;
        double getLeftLogMoneynessAsymptoticBehavior(double epsilon) const;
        double getRightLogMoneynessAsymptoticBehavior(double epsilon) const;
        
    private: 
        const double a_; 
        const double b_;
        const double p_; 
        const double m_; 
        const double s_; 
        const double dbdt_; 
        const double dmdt_; 
        const double dsdt_;
        const double dadt_;   
        const double timeToMaturity_;

        double getParameterB(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterRho(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterBeta(double atmVariance, double atmSkew, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterAlpha(double atmVariance, double atmSkew, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterM(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterA(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterS(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterdBdT(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterdMdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterdSdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getParameterdAdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const;
        double getRawParameterdBdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const;
        double getRawParameterdMdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const;
        double getRawParameterdSdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const;
        double getRawParameterdAdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const;
        static std::tuple<double, double, double, double, double, double> processRawParameters(std::tuple<double, double, double, double, double> params, double timeToMaturity);
        
        
};