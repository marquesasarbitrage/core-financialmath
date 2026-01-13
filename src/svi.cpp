#include "../include/core-financialmath/svi.hpp"

double SVI::getImpliedVarianceATM() const { return getTotalVarianceATM()/timeToMaturity_; }

double SVI::getImpliedVolatilityATM() const { return std::sqrt(getImpliedVarianceATM()); }

double SVI::getImpliedVariance(double x) const { return getTotalVariance(x)/timeToMaturity_; }

double SVI::getImpliedVolatility(double x) const { return std::sqrt(getImpliedVariance(x)); }

double SVI::getTotalVariance(double x) const { return a_ + b_*(p_*(x-m_)+std::sqrt((x-m_)*(x-m_)+s_*s_)); }

double SVI::getTotalVarianceATM() const { return getTotalVariance(0.0); }

double SVI::getMinimumVariance() const { return (a_ + b_*s_*std::sqrt(1-p_*p_))/timeToMaturity_ ; }

double SVI::getMinimumTotalVariance() const { return a_ + b_*s_*std::sqrt(1-p_*p_); }

double SVI::getMinimumImpliedVolatility() const { return std::sqrt(getMinimumVariance()); }

double SVI::getSkewATM() const { return b_*(p_ - m_ / std::sqrt(m_*m_ + s_*s_))/(2.0*std::sqrt(getTotalVarianceATM())); }

double SVI::getCallSlope() const { return b_*(1+p_) /std::sqrt(getTotalVarianceATM()); }

double SVI::getPutSlope() const { return b_*(1-p_) /std::sqrt(getTotalVarianceATM()); }

bool SVI::isFreeOfButterflyArbitrage() const {

    double sqwt = std::sqrt(getTotalVarianceATM()),ct = getCallSlope(),pt = getPutSlope();
    return (sqwt*std::max(ct,pt)<2 and ((ct+pt)*std::max(ct,pt)) <=2);
}

void SVI::checkFreeOfButterflyArbitrage() const {

    if (!isFreeOfButterflyArbitrage()) throw FinMathErrorRegistry::SVI::ButterflyArbitrageError();
}

double SVI::getTimeToMaturity() const { return timeToMaturity_; }

double SVI::getRiskNeutralDensity(double x) const {

    double w = getTotalVariance(x);
    double dwdk = b_*(p_+(x-m_)/(sqrt((x-m_)*(x-m_) + s_*s_)));
    double dw2dk2 = b_*s_*s_/std::pow((x-m_)*(x-m_) + s_*s_,1.5);
    return (1 - .5*x*dwdk/w) * (1 - .5*x*dwdk/w) - .25 * dwdk*dwdk * (0.25 + 1/w) + .5 * dw2dk2;
}

double SVI::getLocalVariance(double x) const {
    
    double g = p_*(x-m_)+sqrt((x-m_)*(x-m_)+s_*s_);
    double dgdt = -p_*dmdt_ + (dsdt_*s_-dmdt_*(x-m_))/sqrt((x-m_)*(x-m_)+s_*s_);
    return (dadt_ + b_*dgdt + dbdt_*g)/getRiskNeutralDensity(x);

}

double SVI::getLocalVolatility(double x) const { return std::sqrt(getLocalVariance(x)); }
        
double SVI::getNormalizedBlackPrice(double x, bool isCall) const {
    
    return LetsBeRational::getNormalizedPrice(x,getImpliedVolatility(x)*std::sqrt(timeToMaturity_),isCall);
}

double SVI::getLeftLogMoneynessAsymptoticBehavior(double epsilon) const {

    return m_ - s_ * (1 - epsilon*(1-p_)) / sqrt(1 - (1 - epsilon*(1-p_))*(1 - epsilon*(1-p_)));
}

double SVI::getRightLogMoneynessAsymptoticBehavior(double epsilon) const {

    return m_ + s_ * (1 - epsilon*(1 + p_)) / sqrt(1 - (1 - epsilon*(1 + p_))*(1 - epsilon*(1 + p_)));
}

double SVI::getParameterB(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double b = std::sqrt(atmVariance*timeToMaturity)*(callSlope+putSlope)/2.0;
    if (b_<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    return b;
}

double SVI::getParameterRho(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
    double p = (b == 0) ? 0.0 : 1 - putSlope*sqrt(atmVariance*timeToMaturity)/b;
    if (abs(p)>1.0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    return p;
}

double SVI::getParameterBeta(double atmVariance, double atmSkew, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
    double p = getParameterRho(atmVariance,callSlope,putSlope,timeToMaturity);
    double beta = b == 0.0 ? .5 : p - 2*atmSkew*sqrt(atmVariance*timeToMaturity)/b;
    if (abs(beta)>1.0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    return beta;
}

double SVI::getParameterAlpha(double atmVariance, double atmSkew, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double beta = getParameterBeta(atmVariance, atmSkew, callSlope, putSlope, timeToMaturity);
    return beta == 0.0 ? 0.0 : beta > 0 ? std::sqrt(1/(beta*beta) - 1) : -std::sqrt(1/(beta*beta) - 1);
}

double SVI::getParameterM(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
    if (b==0.0) return 0.0; 
    double p = getParameterRho(atmVariance,callSlope,putSlope,timeToMaturity);
    double alpha = getParameterAlpha(atmVariance,atmSkew,callSlope,putSlope,timeToMaturity);
    double sqAlpha = alpha < 0 ? -std::sqrt(1+alpha*alpha) : std::sqrt(1+alpha*alpha); 
    return timeToMaturity*(atmVariance-minimumVariance)/(b*(-p+sqAlpha-alpha*std::sqrt(1-p*p)));
}

double SVI::getParameterA(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double m = getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity);
    double p = getParameterRho(atmVariance,callSlope,putSlope,timeToMaturity);
    if (m==0.0)  {

        return timeToMaturity*(minimumVariance+atmVariance*std::sqrt(1-p*p))/(1-std::sqrt(1-p*p));
    }
    else {

        double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
        double alpha = getParameterAlpha(atmVariance,atmSkew,callSlope,putSlope,timeToMaturity);
        return timeToMaturity*minimumVariance-b*alpha*m*sqrt(1-p*p);
    }
}

double SVI::getParameterS(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double s;
    double m = getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity); 
    double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity); 
    if ( m==0.0 ) {

        double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity); 
        s = (b == 0.0) ? 1.0 : (atmVariance*timeToMaturity - getParameterA(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))/b ;

    } else s = m * getParameterAlpha(atmVariance,atmSkew,callSlope,putSlope ,timeToMaturity);
    if (s<=0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    return s;
}

double SVI::getParameterdBdT(double atmVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    return std::sqrt(atmVariance)*(callSlope+putSlope)/(4*std::sqrt(timeToMaturity));
}

double SVI::getParameterdMdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
    
    if ( b_==0 ) { return 0.0; } 
    else {

        double m = getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity); 
        double dbdt = getParameterdBdT(atmVariance,callSlope,putSlope,timeToMaturity);
        return m * (b - timeToMaturity*dbdt)/(timeToMaturity*b);
    }
}

double SVI::getParameterdSdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double m = getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity); 
    if ( m_==0.0 ){

        double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
        double a = getParameterA(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity);
        double dbdt = getParameterdBdT(atmVariance,callSlope,putSlope,timeToMaturity);
        return b==0.0 ? 0.0 : (b*(atmVariance-a/timeToMaturity)+dbdt*(atmVariance*timeToMaturity-a))/(b*b);

    } else {

        double alpha = getParameterAlpha(atmVariance,atmSkew,callSlope,putSlope ,timeToMaturity);
        double dmdt = getParameterdMdT(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity);
        return alpha*dmdt;
    }
}

double SVI::getParameterdAdT(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity) const {

    FinmathToolbox::checkImpliedVolatility(atmVariance);
    FinmathToolbox::checkImpliedVolatility(minimumVariance);
    FinmathToolbox::checkYearFraction(timeToMaturity);
    double m = getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity); 
    if ( m==0.0 ) return  getParameterA(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity)/timeToMaturity;
    else {
        double b = getParameterB(atmVariance,callSlope,putSlope,timeToMaturity);
        double s = getParameterS(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity);
        double p = getParameterRho(atmVariance,callSlope,putSlope,timeToMaturity);
        double dbdt = getParameterdBdT(atmVariance,callSlope,putSlope,timeToMaturity);
        double dsdt = getParameterdSdT(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity);
        return minimumVariance-sqrt(1-p*p)*(dbdt*s+dsdt*b);
    }
}

double SVI::getRawParameterdBdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const {

    auto [a,b,m,s,p,t] = processRawParameters(params,timeToMaturity);
    double atmVariance = (a + b * (-p*m + std::sqrt(m*m + s*s)))/t; 
    double callSlope = b*(1+p) /std::sqrt(atmVariance); 
    double putSlope = b*(1-p) /std::sqrt(atmVariance); 
    return std::sqrt(atmVariance)*(callSlope+putSlope)/(4*std::sqrt(t));
}

double SVI::getRawParameterdMdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const {

    auto [a,b,m,s,p,t] = processRawParameters(params,timeToMaturity);
    
    if ( b_==0 ) { return 0.0; } 
    else {

        double dbdt = getRawParameterdBdT(params,t);
        return m * (b - t*dbdt)/(t*b);
    }
}

double SVI::getRawParameterdSdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const {

    auto [a,b,m,s,p,t] = processRawParameters(params,timeToMaturity);
    double atmVariance = (a + b * (-p*m + std::sqrt(m*m + s*s)))/t; 
    if ( m_==0.0 ){

        double dbdt = getRawParameterdBdT(params,timeToMaturity);
        return b==0.0 ? 0.0 : (b*(atmVariance-a/t)+dbdt*(atmVariance*t-a))/(b*b);

    } else {

        double skewATM = b*(p - m / std::sqrt(m*m + s*s))/(2.0*std::sqrt(atmVariance));
        double beta = b == 0.0 ? .5 : p - 2*skewATM*sqrt(atmVariance*t)/b;
        if (abs(beta)>1.0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
        double alpha = (beta == 0.0) ? 0.0 : beta > 0 ? std::sqrt(1/(beta*beta) - 1) : -std::sqrt(1/(beta*beta) - 1);
        double dmdt = getRawParameterdMdT(params, timeToMaturity);
        return alpha*dmdt;
    }
}

double SVI::getRawParameterdAdT(std::tuple<double, double, double, double, double> params, double timeToMaturity) const {

    auto [a,b,m,s,p,t] = processRawParameters(params,timeToMaturity);

    if ( m==0.0 ) return  a/timeToMaturity;
    else {
        double dbdt = getRawParameterdBdT(params,timeToMaturity);
        double dsdt = getRawParameterdSdT(params, timeToMaturity);
        return (a + b*s*std::sqrt(1-p*p))/timeToMaturity-sqrt(1-p*p)*(dbdt*s+dsdt*b);
    }
}

std::tuple<double, double, double, double, double, double> SVI::processRawParameters(std::tuple<double, double, double, double, double> params, double timeToMaturity) {

    auto [a,b,m,s,p] = params;
    FinmathToolbox::checkYearFraction(timeToMaturity);
    if (s<=0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if (b<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if (abs(p)>1.0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if(a+b*s*sqrt(1-p*p)<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    return {a,b,m,s,p,timeToMaturity};
}

SVI::SVI(double atmVariance, double atmSkew, double minimumVariance, double callSlope, double putSlope, double timeToMaturity): 
a_(getParameterA(atmVariance, atmSkew,minimumVariance,callSlope,putSlope,timeToMaturity))
, b_(getParameterB(atmVariance,callSlope,putSlope,timeToMaturity))
, p_(getParameterRho(atmVariance,callSlope,putSlope,timeToMaturity))
, m_(getParameterM(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))
, s_(getParameterS(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))
, dbdt_(getParameterdBdT(atmVariance,callSlope,putSlope,timeToMaturity))
, dmdt_(getParameterdMdT(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))
, dsdt_(getParameterdSdT(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))
, dadt_(getParameterdAdT(atmVariance, atmSkew, minimumVariance, callSlope, putSlope, timeToMaturity))
, timeToMaturity_(timeToMaturity) {

    if(a_+b_*s_*sqrt(1-p_*p_)<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
}

SVI SVI::getSVIFromRawParameters(double a, double b, double m, double s, double p, double timeToMaturity) {

    FinmathToolbox::checkYearFraction(timeToMaturity);
    if (s<=0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if (b<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if (abs(p)>1.0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    if(a+b*s*sqrt(1-p*p)<0) throw FinMathErrorRegistry::SVI::InvalidSVIParamError();
    std::cout << "ok" << std::endl;
    double atmVariance = (a + b * (-p*m + std::sqrt(m*m + s*s)))/timeToMaturity; 
    double minVariance = (a + b*s*std::sqrt(1-p*p))/timeToMaturity;
    double skewATM = b*(p - m / std::sqrt(m*m + s*s))/(2.0*std::sqrt(atmVariance));
    double callSlope = b*(1+p) /std::sqrt(atmVariance); 
    double putSlope = b*(1-p) /std::sqrt(atmVariance); 
    return SVI(atmVariance,skewATM,minVariance,callSlope,putSlope,timeToMaturity);
}

SVI::SVI(std::tuple<double, double, double, double, double> params, double timeToMaturity): 
a_(std::get<0>(params))
, b_(std::get<1>(params))
, p_(std::get<4>(params))
, m_(std::get<2>(params))
, s_(std::get<3>(params))
, dbdt_(getRawParameterdBdT(params,timeToMaturity))
, dmdt_(getRawParameterdMdT(params, timeToMaturity))
, dsdt_(getRawParameterdSdT(params, timeToMaturity))
, dadt_(getRawParameterdAdT(params, timeToMaturity))
, timeToMaturity_(timeToMaturity) {}

std::tuple<double, double, double, double, double> SVI::getRawParameters() const { return {a_,b_,m_,s_,p_}; }

