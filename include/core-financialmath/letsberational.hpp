#pragma once 
#include <iostream>
#include <cfloat>
#include "errors/main.hpp"
#include "core-math/optim/newtonraphson.hpp"
#include "core-math/probability/univariate.hpp"

// Model references 
// The pricing of options and corporate liabilities - Black and Scholes (1973) : https://www.jstor.org/stable/1831029
// The pricing of commodity contract - Black (1976) : https://www.sciencedirect.com/science/article/abs/pii/0304405X76900246
// Let's be rational - Jäckel (2015) : http://www.jaeckel.org/

class LetsBeRational {

    public: 

        static constexpr double PI = 3.14159265358979323846;

        struct ImpliedVolatilityResult {

            double value_; 
            std::exception_ptr error_;
            int iterations_;
            double targetValue_;
        };

        static double getNormalizedPrice(double x, double normalizedSigma, bool isCall);
        static double getNormalizedPrice(double x, double timeToMaturity, double sigma, bool isCall);
        static double getPrice(double futurePrice, double strike, double timeToMaturity, double sigma, bool isCall);
        static double getNormalizedIntrisicValue(double x, bool isCall);
        static ImpliedVolatilityResult getImpliedVolatility(double normalizedPrice, double x, double timeToMaturity, bool isCall);
        static ImpliedVolatilityResult getImpliedVolatility(double price, double futurePrice, double strike, double timeToMaturity, bool isCall); 
        static ImpliedVolatilityResult getImpliedNormalizedVolatility(double normalizedPrice, double x, bool isCall);


    private: 
        static constexpr double H_LARGE = -10.0;
        static constexpr double T_SMALL = 0.21;
        static constexpr double SQRT_TWO_PI = 2.50662827463100050242;
        static constexpr double ONE_OVER_SQRT_TWO = 0.7071067811865475244008443621048490392848359376887;
        static constexpr double ONE_OVER_SQRT_TWO_PI = 0.3989422804014326779399460599343818684758586311649;
        // Some of the code is sourced from https://github.com/vollib/lets_be_rational/
        static double _getCallPriceRegion1(double h, double t); 
        static double _getCallPriceRegion2(double h, double t); 
        static double _getCallPriceRegion3(double h, double t); 
        static double _getCallPriceRegion4(double h, double t); 
        static double _getCallPrice(double x, double normalizedSigma);
        static double _getVega(double x, double normalizedSigma);
        static double _getVolga(double x, double normalizedSigma);
        static double _getNormalizedPrice(double x, double normalizedSigma, bool isCall);

        static double _getRationalCubicInterpolate(double x, double x0, double x1,double y0, double y1,double dy0, double dy1,double r); 

        static double _getRightAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy0);
        static double _getRightAsymptotic(double sigma); 
        static double _getFirstDerivativeRightAsymptotic(double x, double sigma);
        static double _getSecondDerivativeRightAsymptotic(double x, double sigma);

        static double _getLeftAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy1);
        static double _getLeftAsymptotic(double x, double sigma); 
        static double _getFirstDerivativeLeftAsymptotic(double x, double sigma); 
        static double _getSecondDerivativeLeftAsymptotic(double x, double sigma); 
        static double _getLeftAsymptoticZ(double x, double sigma);

        static std::tuple<double, double, double, double, double, double, double, double, double, double, bool> _getInitialData(double beta, double x, double isCall);
        static double _getInitialGuess(double beta, double x, bool isCall); 
        static std::function<double(double)> _getTarget(double beta, double x, double bLower, double bUpper);
        static std::function<double(double)> _getTargetFirstDerivative(double beta, double x, double bLower, double bUpper); 

        static std::tuple<double, double, int> _getNewtonNormalizedVolatilityUnsafe(double beta, double x, bool isCall);
        
};

