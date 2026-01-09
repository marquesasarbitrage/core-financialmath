#include "../include/core-financialmath/baroneadesiwhaley.hpp"

double BaroneAdesiWhaley::getQ(double drift, double sigma, double interestRate, double timeToMaturity, bool isCall) {

    double putCallFlag = isCall ? 1.0 : -1.0; 
    double N = 2*drift/(sigma*sigma);
    double M = 2*interestRate/(sigma*sigma);
    return .5*(1-N+putCallFlag*sqrt((N-1)*(N-1)+4*M/(1-getDiscountFactor(timeToMaturity,interestRate))));
}

double BaroneAdesiWhaley::getQInf(double drift, double sigma, double interestRate, bool isCall) {

    double putCallFlag = isCall ? 1.0 : -1.0; 
    double N = 2*drift/(sigma*sigma);
    double M = 2*interestRate/(sigma*sigma);
    return .5*(1-N+putCallFlag*sqrt((N-1)*(N-1)+4*M));
}

double BaroneAdesiWhaley::getDiscountFactor(double timeToMaturity, double interestRate) {

    return std::exp(-interestRate*timeToMaturity);
}

double BaroneAdesiWhaley::getEuroPrice(double futurePrice, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    return getDiscountFactor(timeToMaturity,interestRate)*LetsBeRational::getPrice(futurePrice,strike,timeToMaturity,sigma,isCall);
}

double BaroneAdesiWhaley::getEuroDelta(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    int putCallFlag = isCall ? 1 : -1; 
    double d1 = (std::log(futurePrice / strike) + 0.5 * sigma * sigma * timeToMaturity) / (sigma * std::sqrt(timeToMaturity)); 
    return getDiscountFactor(timeToMaturity,interestRate)*putCallFlag*exp(drift*timeToMaturity)*Gaussian().cdf(putCallFlag*d1);
}

double BaroneAdesiWhaley::getEuroGamma(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    double d1 = (std::log(futurePrice / strike) + 0.5 * sigma * sigma * timeToMaturity) / (sigma * std::sqrt(timeToMaturity)); 
    double driftSq = exp(drift*timeToMaturity)*exp(drift*timeToMaturity);
    return getDiscountFactor(timeToMaturity,interestRate)*driftSq*Gaussian().pdf(d1)/(futurePrice*sigma*std::sqrt(timeToMaturity));
}

double BaroneAdesiWhaley::getInitialGuessOptimalExercisePrice(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    double putCallFlag = isCall ? 1.0 : -1.0; 
    double qq_inf_ = getQInf(drift,sigma,interestRate,isCall);
    double Sinf = strike/(1 - 1/qq_inf_);
    double payoff = putCallFlag*(Sinf - strike); 
    double h = -putCallFlag*(drift*timeToMaturity+putCallFlag*2*sigma*std::sqrt(timeToMaturity))*strike/payoff; 
    return Sinf-putCallFlag*std::exp(h)*payoff;
}

std::function<double(double)> BaroneAdesiWhaley::_getTarget(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    return [drift, strike, timeToMaturity, sigma, interestRate, isCall](double F) {

        double putCallFlag = isCall ? 1.0 : -1.0;
        double q = getQ(drift,sigma,interestRate,timeToMaturity,isCall);
        double price = getEuroPrice(F*std::exp(drift*timeToMaturity),strike,timeToMaturity,sigma,interestRate,isCall);
        double delta = getEuroDelta(F*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
        return putCallFlag*(F-strike) - price - putCallFlag*F*(1-putCallFlag*delta)/q;

    };
}

std::function<double(double)> BaroneAdesiWhaley::_getTargetFirstDerivative(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    return [drift, strike, timeToMaturity, sigma, interestRate, isCall](double F) {

        double putCallFlag = isCall ? 1.0 : -1.0;
        double q = getQ(drift,sigma,interestRate,timeToMaturity,isCall);
        double delta = getEuroDelta(F*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
        double gamma = getEuroGamma(F*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
        return putCallFlag - delta*(1-1/q) - (putCallFlag-gamma*F)/q;

    };
}

double BaroneAdesiWhaley::_getOptimalExercisePrice(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    if (drift >= interestRate and isCall){return DBL_MAX;}
    NewtonRaphson newton(
        getInitialGuessOptimalExercisePrice(drift, strike,timeToMaturity,sigma,interestRate,isCall),
        _getTarget(drift, strike,timeToMaturity,sigma,interestRate,isCall),
        _getTargetFirstDerivative(drift, strike,timeToMaturity,sigma,interestRate,isCall)
    );
    newton.setToleranceThreshold(1e-20);
    newton.setMaximumIterations(20);
    newton.optimize(); 
    return newton.getResult();
}

double BaroneAdesiWhaley::getExercisePremium(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    double putCallFlag = isCall ? 1.0 : -1.0;
    double exercisePrice = _getOptimalExercisePrice(drift,strike,timeToMaturity,sigma,interestRate,isCall);
    double qq_ = getQ(drift,sigma,interestRate,timeToMaturity,isCall);
    double delta = getEuroDelta(exercisePrice*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
    return putCallFlag * exercisePrice * std::pow(futurePrice / exercisePrice, qq_) * (1 - putCallFlag * delta) / qq_;
}

double BaroneAdesiWhaley::getPrice(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    double putCallFlag = isCall ? 1.0 : -1.0;
    double qq_ = getQ(drift,sigma,interestRate,timeToMaturity,isCall);
    double exercisePrice = _getOptimalExercisePrice(drift, strike,timeToMaturity,sigma,interestRate,isCall);
    double delta = getEuroDelta(exercisePrice*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
    double exercisePremium = putCallFlag * exercisePrice * std::pow(futurePrice / exercisePrice, qq_) * (1 - putCallFlag * delta) / qq_;
    double euroPrice = getEuroPrice(futurePrice*std::exp(drift*timeToMaturity),strike,timeToMaturity,sigma,interestRate,isCall);
    return (putCallFlag*futurePrice >= putCallFlag*exercisePrice) ? putCallFlag*(futurePrice-strike) : euroPrice + exercisePremium;
}


