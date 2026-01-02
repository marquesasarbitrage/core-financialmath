#include "../include/core-financialmath/letsberational.hpp"

double LetsBeRational::getPutPrice(double logMoneyness, double timeToMaturity, double sigma) {

    return _getPrice(logMoneyness,sigma*std::sqrt(timeToMaturity), false);
}

double LetsBeRational::getCallPrice(double logMoneyness, double timeToMaturity, double sigma) {

    return _getPrice(logMoneyness,sigma*std::sqrt(timeToMaturity), true);
}

LetsBeRational::ImpliedVolatilityResult LetsBeRational::getNewtonNormalizedVolatility(double beta, double x, bool isCall) {

    ImpliedVolatilityResult result; 
    try {
        std::tuple<double, double, int> tupleResult = _getNewtonNormalizedVolatility(beta,x,isCall);
        
        result.value_ = std::get<0>(tupleResult);
        result.error_ = nullptr;
        result.iterations_ = std::get<2>(tupleResult);
        result.targetValue_ = std::get<1>(tupleResult);
    } catch (const std::exception& e) {
        ImpliedVolatilityResult result; 
        result.value_ = 0.0;
        result.error_ = std::current_exception();
        result.targetValue_ = NAN; 
        result.iterations_ = 0;
    }
    return result;
}

LetsBeRational::ImpliedVolatilityResult LetsBeRational::getCallImpliedVolatility(double logMoneyness, double timeToMaturity, double normalizedUndiscountedPrice) {

    ImpliedVolatilityResult result = getNewtonNormalizedVolatility(normalizedUndiscountedPrice,logMoneyness,true);
    result.value_ = result.value_/std::sqrt(timeToMaturity);
    return result;
}

LetsBeRational::ImpliedVolatilityResult LetsBeRational::getPutImpliedVolatility(double logMoneyness, double timeToMaturity, double normalizedUndiscountedPrice) {

    ImpliedVolatilityResult result = getNewtonNormalizedVolatility(normalizedUndiscountedPrice,logMoneyness,false);
    result.value_ = result.value_/std::sqrt(timeToMaturity);
    return result;
}

double LetsBeRational::_getCallPriceRegion1(double h, double t) {

    const double e=(t/h)*(t/h), r=((h+t)*(h-t)), q=(h/r)*(h/r);
    const double asymptotic_expansion_sum = (2.0+q*(-6.0E0-2.0*e+3.0*q*(1.0E1+e*(2.0E1+2.0*e)+5.0*q*(-1.4E1+e*(-7.0E1+e*(-4.2E1-2.0*e))+7.0*q*(1.8E1+e*(1.68E2+e*(2.52E2+e*(7.2E1+2.0*e)))+9.0*q*(-2.2E1+e*(-3.3E2+e*(-9.24E2+e*(-6.6E2+e*(-1.1E2-2.0*e))))+1.1E1*q*(2.6E1+e*(5.72E2+e*(2.574E3+e*(3.432E3+e*(1.43E3+e*(1.56E2+2.0*e)))))+1.3E1*q*(-3.0E1+e*(-9.1E2+e*(-6.006E3+e*(-1.287E4+e*(-1.001E4+e*(-2.73E3+e*(-2.1E2-2.0*e))))))+1.5E1*q*(3.4E1+e*(1.36E3+e*(1.2376E4+e*(3.8896E4+e*(4.862E4+e*(2.4752E4+e*(4.76E3+e*(2.72E2+2.0*e)))))))+1.7E1*q*(-3.8E1+e*(-1.938E3+e*(-2.3256E4+e*(-1.00776E5+e*(-1.84756E5+e*(-1.51164E5+e*(-5.4264E4+e*(-7.752E3+e*(-3.42E2-2.0*e))))))))+1.9E1*q*(4.2E1+e*(2.66E3+e*(4.0698E4+e*(2.3256E5+e*(5.8786E5+e*(7.05432E5+e*(4.0698E5+e*(1.08528E5+e*(1.197E4+e*(4.2E2+2.0*e)))))))))+2.1E1*q*(-4.6E1+e*(-3.542E3+e*(-6.7298E4+e*(-4.90314E5+e*(-1.63438E6+e*(-2.704156E6+e*(-2.288132E6+e*(-9.80628E5+e*(-2.01894E5+e*(-1.771E4+e*(-5.06E2-2.0*e))))))))))+2.3E1*q*(5.0E1+e*(4.6E3+e*(1.0626E5+e*(9.614E5+e*(4.08595E6+e*(8.9148E6+e*(1.04006E7+e*(6.53752E6+e*(2.16315E6+e*(3.542E5+e*(2.53E4+e*(6.0E2+2.0*e)))))))))))+2.5E1*q*(-5.4E1+e*(-5.85E3+e*(-1.6146E5+e*(-1.77606E6+e*(-9.37365E6+e*(-2.607579E7+e*(-4.01166E7+e*(-3.476772E7+e*(-1.687257E7+e*(-4.44015E6+e*(-5.9202E5+e*(-3.51E4+e*(-7.02E2-2.0*e))))))))))))+2.7E1*q*(5.8E1+e*(7.308E3+e*(2.3751E5+e*(3.12156E6+e*(2.003001E7+e*(6.919458E7+e*(1.3572783E8+e*(1.5511752E8+e*(1.0379187E8+e*(4.006002E7+e*(8.58429E6+e*(9.5004E5+e*(4.7502E4+e*(8.12E2+2.0*e)))))))))))))+2.9E1*q*(-6.2E1+e*(-8.99E3+e*(-3.39822E5+e*(-5.25915E6+e*(-4.032015E7+e*(-1.6934463E8+e*(-4.1250615E8+e*(-6.0108039E8+e*(-5.3036505E8+e*(-2.8224105E8+e*(-8.870433E7+e*(-1.577745E7+e*(-1.472562E6+e*(-6.293E4+e*(-9.3E2-2.0*e))))))))))))))+3.1E1*q*(6.6E1+e*(1.0912E4+e*(4.74672E5+e*(8.544096E6+e*(7.71342E7+e*(3.8707344E8+e*(1.14633288E9+e*(2.07431664E9+e*(2.33360622E9+e*(1.6376184E9+e*(7.0963464E8+e*(1.8512208E8+e*(2.7768312E7+e*(2.215136E6+e*(8.184E4+e*(1.056E3+2.0*e)))))))))))))))+3.3E1*(-7.0E1+e*(-1.309E4+e*(-6.49264E5+e*(-1.344904E7+e*(-1.4121492E8+e*(-8.344518E8+e*(-2.9526756E9+e*(-6.49588632E9+e*(-9.0751353E9+e*(-8.1198579E9+e*(-4.6399188E9+e*(-1.6689036E9+e*(-3.67158792E8+e*(-4.707164E7+e*(-3.24632E6+e*(-1.0472E5+e*(-1.19E3-2.0*e)))))))))))))))))*q)))))))))))))))));
    return Gaussian::ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*(t/r)*asymptotic_expansion_sum;
}

double LetsBeRational::_getCallPriceRegion2(double h, double t) {

    const double a = 1+h*(0.5*SQRT_TWO_PI)*Gaussian::erfcxCody(-ONE_OVER_SQRT_TWO*h), w=t*t, h2=h*h;
    const double expansion = 2*t*(a+w*((-1+3*a+a*h2)/6+w*((-7+15*a+h2*(-1+10*a+a*h2))/120+w*((-57+105*a+h2*(-18+105*a+h2*(-1+21*a+a*h2)))/5040+w*((-561+945*a+h2*(-285+1260*a+h2*(-33+378*a+h2*(-1+36*a+a*h2))))/362880+w*((-6555+10395*a+h2*(-4680+17325*a+h2*(-840+6930*a+h2*(-52+990*a+h2*(-1+55*a+a*h2)))))/39916800+((-89055+135135*a+h2*(-82845+270270*a+h2*(-20370+135135*a+h2*(-1926+25740*a+h2*(-75+2145*a+h2*(-1+78*a+a*h2))))))*w)/6227020800.0))))));
    return Gaussian::ONE_OVER_SQRT_TWO_PI*exp((-0.5*(h*h+t*t)))*expansion;
} 

double LetsBeRational::_getCallPriceRegion3(double h, double t) {

    Gaussian norm = Gaussian();
    return norm.cdf(h+t)*exp(h*t) - norm.cdf(h-t)/exp(h*t);
}

double LetsBeRational::_getCallPriceRegion4(double h, double t) {
    
    return 0.5 * exp(-0.5*(h*h+t*t)) * (Gaussian::erfcxCody(-ONE_OVER_SQRT_TWO*(h+t)) - Gaussian::erfcxCody(-ONE_OVER_SQRT_TWO*(h-t)));
}

double LetsBeRational::_getCallPrice(double x, double normalizedSigma) {

    if (x>0) { 
        
        return _getNormalizedIntrisicValue(x, true)+_getCallPrice(-x,normalizedSigma);
    }

    if (normalizedSigma<=0) {
        
        return _getNormalizedIntrisicValue(x, true);
    }

    if ( x < normalizedSigma*H_LARGE  &&  0.5*normalizedSigma*normalizedSigma+x < normalizedSigma*(T_SMALL+H_LARGE)) {
        
        return _getCallPriceRegion1(x/normalizedSigma, .5*normalizedSigma);
    }

    if ( 0.5*normalizedSigma < T_SMALL) {

        return _getCallPriceRegion2(x/normalizedSigma, .5*normalizedSigma);
    }
    
    if (x+0.5*normalizedSigma*normalizedSigma > normalizedSigma*0.85) {
        
        return _getCallPriceRegion3(x/normalizedSigma, .5*normalizedSigma);
    }

    return _getCallPriceRegion4(x/normalizedSigma, .5*normalizedSigma);
}

double LetsBeRational::_getVega(double x, double normalizedSigma) {

    double ax = fabs(x);
    double nsSq = normalizedSigma*normalizedSigma;
    return (ax<=0) ? ONE_OVER_SQRT_TWO_PI*exp(-0.125*nsSq) : ( (normalizedSigma<=0 || normalizedSigma<=ax*std::sqrt(DBL_MIN)) ? 0 : ONE_OVER_SQRT_TWO_PI*exp(-0.5*(x*x/nsSq+.25*nsSq)));
}

double LetsBeRational::_getVolga(double x, double normalizedSigma) {

    if (normalizedSigma<=DBL_MIN) return 0;
    return (x/(normalizedSigma*normalizedSigma*normalizedSigma) - .125)*_getVega(x, normalizedSigma);
}

double LetsBeRational::_getPrice(double x, double normalizedSigma, bool isCall) {

    if (normalizedSigma<=DBL_MIN) return _getNormalizedIntrisicValue(x, isCall);
    if (!isCall) x = -x;
    return _getCallPrice(x,normalizedSigma);
}

double LetsBeRational::_getNormalizedIntrisicValue(double x, bool isCall) {

    int putCallFlag = isCall ? 1:-1;
    if (putCallFlag*x <= 0) return 0;
    double bm = exp(.5*x);
    return std::max(putCallFlag*(bm - 1/bm), 0.0);
}

double LetsBeRational::_getRationalCubicInterpolate(double x, double x0, double x1,double y0, double y1,double dy0, double dy1,double r) {

    double h = x1 - x0;
    double s = (x - x0) / h;

    double numerator =
        y1 * s * s * s + 
        (r * y1 - h * dy1) * s * s * (1 - s) +
        (r * y0 + h * dy0) * s * (1 - s) * (1 - s) +
        y0 * (1 - s) * (1 - s) * (1 - s);

    double denominator = 1 + (r - 3) * s * (1 - s);

    return numerator / denominator;
} 

double LetsBeRational::_getRightAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy0) {

    double denominator = ((y1 - y0) / (x1 - x0) - dy0);
    double numerator =  (.5*(x1 - x0)*ddy0 + dy1 - dy0);
    if (abs(denominator) < DBL_MIN) return (numerator > 0) ? 2 / (DBL_EPSILON * DBL_EPSILON) : std::sqrt(DBL_EPSILON) - 1.0;
    else return numerator / denominator;
}

double LetsBeRational::_getRightAsymptotic(double sigma) { return Gaussian().cdf(-sigma/2); }

double LetsBeRational::_getFirstDerivativeRightAsymptotic(double x, double sigma) {

    return (abs(x)<DBL_MIN) ? -.5 : -.5 * exp(.5*x*x/(sigma*sigma));
}

double LetsBeRational::_getSecondDerivativeRightAsymptotic(double x, double sigma) {

    return (abs(x)<DBL_MIN) ? 0 : exp(sigma*sigma / 8 +  x*x/sigma*sigma)*sqrt(PI/2) * x * x / (sigma*sigma*sigma);
}

double LetsBeRational::_getLeftAsymptoticR(double x0, double x1,double y0, double y1,double dy0, double dy1, double ddy1) {

    double denominator = (dy1 - (y1 - y0) / (x1 - x0));
    double numerator =  (.5*(x1 - x0)*ddy1 + dy1 - dy0);
    if (abs(denominator) < DBL_MIN) return (numerator > 0) ? 2 / (DBL_EPSILON * DBL_EPSILON) : std::sqrt(DBL_EPSILON) - 1.0;
    else return numerator / denominator;
}

double LetsBeRational::_getLeftAsymptotic(double x, double sigma) {

    double cdfz = Gaussian().cdf(_getLeftAsymptoticZ(x, sigma));
    return (2*PI*cdfz*cdfz*cdfz*abs(x))/(3*sqrt(3));
}

double LetsBeRational::_getFirstDerivativeLeftAsymptotic(double x, double sigma) {

    double z = _getLeftAsymptoticZ(x, sigma); double cdfz = Gaussian().cdf(z);
    return 2*PI*z*z*cdfz*cdfz*exp(z*z+sigma*sigma/8.0);
}

double LetsBeRational::_getSecondDerivativeLeftAsymptotic(double x, double sigma) {

    double z = _getLeftAsymptoticZ(x, sigma); double cdfz = Gaussian().cdf(z); double pdfz = Gaussian().pdf(z);
    double sigmasq = sigma*sigma;
    double term1 = (PI/6)*(z*z/(sigma*sigma*sigma))*cdfz*exp(2*z*z+sigma*sigma/4.0);
    double term2 = 8*sqrt(3)*sigma*abs(x)+(3*sigma*sigma*(sigma*sigma-8)-8*x*x)*cdfz/pdfz;
    return term1*term2;
}

double LetsBeRational::_getLeftAsymptoticZ(double x, double sigma) { return -abs(x)/(sqrt(3)*sigma); }

std::tuple<double, double, double, double, double, double, double, double, double, double, bool> LetsBeRational::_getInitialData(double beta, double x, double isCall) {
 
    int c = isCall ? 1 : -1;
    if (c*x>0) {
      beta = fabs(std::max(beta-_getNormalizedIntrisicValue(x, isCall),0.));
      c = -c; isCall = !isCall;
    }
    if (c<0) { x = -x; c = -c; isCall = !isCall; }
    
    double b_max = exp(.5*x);
    double sigma_c = sqrt(2*abs(x));
    double b_c = _getPrice(x, sigma_c, isCall);
    double v_c = _getVega(x, sigma_c);
    double sigma_u = (v_c>DBL_MIN) ? sigma_c + (b_max - b_c)/v_c : sigma_c;
    double sigma_l = (v_c>DBL_MIN) ? sigma_c - b_c/v_c : sigma_c;
    double b_u = _getPrice(x, sigma_u, isCall);
    double b_l = _getPrice(x, sigma_l, isCall);
    return std::make_tuple(beta,x,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, isCall); 
}

double LetsBeRational::_getInitialGuess(double beta, double x, bool isCall){

    auto [beta_,x_,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, is_call_] = _getInitialData(beta, x, isCall);

    if (beta_<=0 or beta>=b_max) throw FinMathErrorRegistry::LetsBeRational::InvalidNormalizedPriceRangeError();

    if (beta_ < b_l && beta_ > 0) {

        double fl = _getLeftAsymptotic(x_, sigma_l);
        double dfl = _getFirstDerivativeLeftAsymptotic(x_, sigma_l);
        double ddfl = _getSecondDerivativeLeftAsymptotic(x_, sigma_l);
        double r = _getLeftAsymptoticR(0, b_l, 0, fl, 1, dfl, ddfl);
        double frc = _getRationalCubicInterpolate(beta_, 0, b_l, 0, fl, 1, dfl, r);
        double sq3 = sqrt(3), a = frc/(2*PI*abs(x_));
        return abs(x_/(sq3*Gaussian().quantile(sq3*std::pow(a, 1.0/3.0))));
    }
    else if (beta_ <= b_c && beta_ >= b_l){

        double v_l = _getVega(x_, sigma_l);
        double r = _getLeftAsymptoticR(b_l, b_c, sigma_l, sigma_c, 1/v_l, 1/v_c, 0);
        return _getRationalCubicInterpolate(beta_, b_l, b_c, sigma_l, sigma_c, 1/v_l, 1/v_c, r);
    }
    else if (beta_ <= b_u && beta_ > b_c){

        double v_u = _getVega(x_, sigma_u);
        double r = _getRightAsymptoticR(b_c, b_u, sigma_c, sigma_u, 1/v_c, 1/v_u, 0);
        return _getRationalCubicInterpolate(beta_, b_c, b_u, sigma_c, sigma_u, 1/v_c, 1/v_u, r);
    }
    else {
        double fu =_getRightAsymptotic(sigma_u);
        double dfu = _getFirstDerivativeRightAsymptotic(x_, sigma_u);
        double ddfu = _getSecondDerivativeRightAsymptotic(x_, sigma_u);
        double r = _getRightAsymptoticR(b_u, b_max, fu, 0, dfu, -.5, ddfu);
        double frc = _getRationalCubicInterpolate(beta_, b_u, b_max, fu, 0, dfu, -.5, r);
        return -2.0*Gaussian().quantile(frc);
    }
} 

std::function<double(double)> LetsBeRational::_getTarget(double beta, double x, double bLower, double bUpper) {

    double bMax = exp(.5*x);
    double bUpper2 = std::max(bUpper, bMax/2);
    if (beta > 0 && beta < bLower) return [beta, x](double s) {return log(1/_getCallPrice(x, s)) - log(1/beta);};
    else if (beta >= bLower && beta <= bUpper2) return [beta, x](double s) {return _getCallPrice(x, s) - beta;};
    else return [bMax, beta, x](double s) {return log((bMax-beta)/(bMax-_getCallPrice(x, s)));};
}

std::function<double(double)> LetsBeRational::_getTargetFirstDerivative(double beta, double x, double bLower, double bUpper) {

    double bMax = exp(.5*x);
    double bUpper2 = std::max(bUpper, bMax/2);
    if (beta > 0 && beta < bLower) return [x](double s) {return -_getVega(x, s)/_getCallPrice(x, s);};
    else if (beta >= bLower && beta <= bUpper2) return [x](double s) {return _getVega(x, s);};
    else return [bMax, x](double s) {return _getVega(x, s)/(bMax - _getCallPrice(x, s));};
}

std::tuple<double, double, int> LetsBeRational::_getNewtonNormalizedVolatility(double beta, double x, bool isCall) {

    auto [beta_,x_,b_max,sigma_l, sigma_c, sigma_u, b_l, b_c, b_u, v_c, is_call_] = _getInitialData(beta, x, isCall);
    
    if (beta_<=0 or beta>=b_max) throw FinMathErrorRegistry::LetsBeRational::InvalidNormalizedPriceRangeError();
    
    NewtonRaphson newton(
        _getInitialGuess(beta, x, isCall),
        _getTarget(beta_, x_, b_l, b_u),
        _getTargetFirstDerivative(beta_, x_, b_l, b_u)
    );
    newton.setToleranceThreshold(1e-20);
    newton.setMaximumIterations(20);
    newton.optimize(); 
    
    if (!newton.getError()) return {newton.getResult(), newton.getFunctionResult(), newton.getNumberIterations()};
    else std::rethrow_exception(newton.getError());
}

