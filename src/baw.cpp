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

double BaroneAdesiWhaley::getEuroVega(double futurePrice, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    double d1 = (std::log(futurePrice / strike) + 0.5 * sigma * sigma * timeToMaturity) / (sigma * std::sqrt(timeToMaturity)); 
    return futurePrice*getDiscountFactor(timeToMaturity,interestRate)*Gaussian().pdf(d1)*sqrt(timeToMaturity);
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

double BaroneAdesiWhaley::getOptimalExercisePrice(double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    if (drift >= interestRate and isCall){return DBL_MAX;}
    NewtonRaphson newton(
        getInitialGuessOptimalExercisePrice(drift, strike,timeToMaturity,sigma,interestRate,isCall),
        _getTarget(drift, strike,timeToMaturity,sigma,interestRate,isCall),
        _getTargetFirstDerivative(drift, strike,timeToMaturity,sigma,interestRate,isCall)
    );
    newton.setToleranceThreshold(1e-20);
    newton.setMaximumIterations(5);
    newton.optimize(); 
    return newton.getResult();
}

double BaroneAdesiWhaley::getExercisePremium(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    return _getExercisePremium(
        futurePrice
        ,drift
        ,strike
        ,timeToMaturity
        ,sigma
        ,interestRate
        ,isCall
        ,getOptimalExercisePrice(drift,strike,timeToMaturity,sigma,interestRate,isCall));
}

double BaroneAdesiWhaley::getPrice(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall) {

    return _getPrice(
        futurePrice
        ,drift
        ,strike
        ,timeToMaturity
        ,sigma
        ,interestRate
        ,isCall
        ,getOptimalExercisePrice(drift,strike,timeToMaturity,sigma,interestRate,isCall));
}

double BaroneAdesiWhaley::_getExercisePremium(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice) {

    if (drift >= interestRate and isCall){return 0.0;}
    double putCallFlag = isCall ? 1.0 : -1.0;
    double qq_ = getQ(drift,sigma,interestRate,timeToMaturity,isCall);
    double delta = getEuroDelta(exercisePrice*std::exp(drift*timeToMaturity),drift,strike,timeToMaturity,sigma,interestRate,isCall);
    return putCallFlag * exercisePrice * std::pow(futurePrice / exercisePrice, qq_) * (1 - putCallFlag * delta) / qq_;
}

double BaroneAdesiWhaley::_getPrice(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice) {

    double putCallFlag = isCall ? 1.0 : -1.0;
    double euroPrice = getEuroPrice(futurePrice*std::exp(drift*timeToMaturity),strike,timeToMaturity,sigma,interestRate,isCall);

    if (drift >= interestRate and isCall){ 

        return euroPrice; 

    } else if (putCallFlag*futurePrice >= putCallFlag*exercisePrice) { 

        return putCallFlag*(futurePrice-strike); 

    } else {

        return euroPrice + _getExercisePremium(futurePrice,drift,strike,timeToMaturity,sigma,interestRate,isCall,exercisePrice);

    }
}

double BaroneAdesiWhaley::_getVega(double futurePrice, double drift, double strike, double timeToMaturity, double sigma, double interestRate, bool isCall, double exercisePrice1, double exercisePrice2, double epsilon) {

    double putCallFlag = isCall ? 1.0 : -1.0;
    double euroVega = getEuroVega(futurePrice*std::exp(drift*timeToMaturity),strike,timeToMaturity,sigma,interestRate,isCall);

    if (drift >= interestRate and isCall){ 

        return euroVega; 

    } else if (putCallFlag*futurePrice >= putCallFlag*exercisePrice1) {

        return 0.0; 
    
    } else {

        double premium1 = _getExercisePremium(futurePrice,drift,strike,timeToMaturity,sigma,interestRate,isCall,exercisePrice1);
        double premium2 = _getExercisePremium(futurePrice,drift,strike,timeToMaturity,sigma+epsilon,interestRate,isCall,exercisePrice2);
        return euroVega + std::max((premium2-premium1)/epsilon,0.0);
    }
    
}

LetsBeRational::ImpliedVolatilityResult BaroneAdesiWhaley::_getImpliedVolatilityNewtonRaphson(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall) {

    LetsBeRational::ImpliedVolatilityResult result;
    double x = LetsBeRational::getInitialGuessImpliedVolatility(
        price/getDiscountFactor(timeToMaturity,interestRate)
        , futurePrice
        , strike
        , timeToMaturity
        , isCall
    );
    double exercisePrice1, exercisePrice2, vega, calculatedPrice, epsilon = 0.001, fX, dfX, xStep = 0.0, threshold = 1e-5;
    int numberIterations_ = 0, maxIterations = 10;
    for (int i = 1; i <= maxIterations; ++i) {

        numberIterations_ += 1;
        exercisePrice1 = getOptimalExercisePrice(0.0,strike,timeToMaturity,x,interestRate,isCall);
        exercisePrice2 = getOptimalExercisePrice(0.0,strike,timeToMaturity,x+epsilon,interestRate,isCall);
        calculatedPrice = _getPrice(futurePrice,0.0,strike,timeToMaturity,x,interestRate,isCall,exercisePrice1); 
        vega = _getVega(futurePrice,0.0,strike,timeToMaturity,x,interestRate,isCall,exercisePrice1,exercisePrice2,epsilon);
        if (vega < __DBL_MIN__) { x *= 3; maxIterations += 10;  continue; }
        fX = std::log(1/calculatedPrice) - std::log(1/price);
        dfX = -vega/calculatedPrice;

        if (std::abs(fX) < threshold) break;
        
        if (std::abs(dfX) < std::sqrt(__DBL_MIN__)) {

            result.value_ = 0.0; 
            result.error_ = nullptr; 
            result.iterations_= numberIterations_;
            result.targetValue_ = fX;
        }
        
        double xNew = x - fX / dfX;
        xStep = std::abs(xNew - x);
        x = xNew;

        if (xStep < threshold) break;

        if (i==maxIterations) break;
        
    }
    result.value_ = x; 
    result.error_ = nullptr; 
    result.iterations_= numberIterations_;
    result.targetValue_ = fX;

    return result;
}

LetsBeRational::ImpliedVolatilityResult BaroneAdesiWhaley::getImpliedVolatility(double price, double futurePrice, double strike, double timeToMaturity, double interestRate, bool isCall) {

    try {
        
        if (price<=0.0) throw FinMathErrorRegistry::LetsBeRational::NegativeOptionPriceError();
        if (!isCall and std::abs((strike - futurePrice)-price) < DBL_MIN) throw FinMathErrorRegistry::LetsBeRational::OptionHasNoTimeValueError();
        FinmathToolbox::checkSpotPrice(futurePrice);
        FinmathToolbox::checkStrikePrice(strike);
        FinmathToolbox::checkYearFraction(timeToMaturity);

        double putCallFlag = (isCall ? 1.0 : -1.0);
        if (0.0 >= interestRate and isCall) { 

            return LetsBeRational::getImpliedVolatility(
                price/getDiscountFactor(timeToMaturity,interestRate)
                , futurePrice
                , strike
                , timeToMaturity
                , isCall);

        } else {

            return _getImpliedVolatilityNewtonRaphson(price,futurePrice,strike,timeToMaturity,interestRate,isCall);
        }

    } catch (const std::exception& e) {

        return {0.0,std::current_exception(),0,NAN};
    }
    
}



