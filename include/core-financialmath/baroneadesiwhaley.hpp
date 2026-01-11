#pragma once 
#include <iostream>
#include "letsberational.hpp"
#include "core-math/optim/neldermead.hpp"

// Model references 
// Efficient analytic approximation of American Option Values - Barone-Adesi and Whaley (1987) : https://doi.org/10.1111/j.1540-6261.1987.tb02569.x

class BaroneAdesiWhaley {

    public: 
        static double getEuroPrice(double futurePrice, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getEuroDelta(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getEuroGamma(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getEuroVega(double futurePrice, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getInitialGuessOptimalExercisePrice(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getOptimalExercisePrice(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getExercisePremium(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getPrice(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static LetsBeRational::ImpliedVolatilityResult getImpliedVolatility(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall);
    
    private: 
        static double getDiscountFactor(double timeToMaturity, double interestRate);
        static std::function<double(double)> _getTarget(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static std::function<double(double)> _getTargetFirstDerivative(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall);
        static double getQ(double drift, double sigma, double interestRate, double timeToMaturity, bool isCall); 
        static double getQInf(double drift, double sigma, double interestRate, bool isCall); 
        static std::function<double(double)> _getTargetImpliedVolatiltiy(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall);
        static std::function<double(double)> _getTargetFirstDerivativeImpliedVolatiltiy(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall);
        static double _getExercisePremium(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice);
        static double _getPrice(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice);
        static double _getVega(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice1, double exercisePrice2, double epsilon);
        static LetsBeRational::ImpliedVolatilityResult _getImpliedVolatilityNewtonRaphson(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall);
};