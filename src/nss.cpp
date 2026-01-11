#include "../include/core-financialmath/nelsonsiegelsvensson.hpp"

double NelsonSiegelFamily::rateFuntion1(double t, double tau){ return tau*(1-std::exp(-t/tau))/t; }
double NelsonSiegelFamily::rateFuntion2(double t, double tau){ return tau*(1-std::exp(-t/tau))/t - std::exp(-t/tau); }
double NelsonSiegelFamily::forwardRateFuntion1(double t, double tau){ return std::exp(-t/tau); }
double NelsonSiegelFamily::forwardRateFuntion2(double t, double tau){ return t*std::exp(-t/tau)/tau; }

NelsonSiegel::NelsonSiegel(double b0, double b1, double b2, double tau): b0_(b0), b1_(b1), b2_(b2), tau_(tau){};
Svensson::Svensson(double b0, double b1, double b2, double b3, double tau1, double tau2): 
b0_(b0), b1_(b1), b2_(b2), b3_(b3), tau1_(tau1), tau2_(tau2){};

double NelsonSiegel::getBeta1() const { return b1_; }
double NelsonSiegel::getBeta2() const  { return b2_; }
double NelsonSiegel::getBeta0() const  { return b0_; }
void NelsonSiegel::setBeta0(double b0) { b0_=b0; }
void NelsonSiegel::setBeta1(double b1) { b1_=b1; }
void NelsonSiegel::setBeta2(double b2) { b2_=b2; }

double Svensson::getBeta1() const { return b1_; }
double Svensson::getBeta2() const  { return b2_; }
double Svensson::getBeta0() const  { return b0_; }
double Svensson::getBeta3() const { return b3_; }
void Svensson::setBeta0(double b0) { b0_=b0; }
void Svensson::setBeta1(double b1) { b1_=b1; }
void Svensson::setBeta2(double b2) { b2_=b2; }
void Svensson::setBeta3(double b3) { b3_ = b3; }

double NelsonSiegel::getTau() const  { return tau_; }
void NelsonSiegel::setTau(double tau) { tau_=tau; }

double Svensson::getTau1() const { return tau1_; }
double Svensson::getTau2() const { return tau2_; }
void Svensson::setTau1(double tau1) { tau1_ = tau1; }
void Svensson::setTau2(double tau2) { tau2_ = tau2; }

double NelsonSiegel::getRate(double t) const {
    return b0_ + b1_ * rateFuntion1(t,tau_) + b2_ * rateFuntion2(t,tau_);
}

double NelsonSiegel::getInstantaneousForwardRate(double t) const {
    return b0_ + b1_ * forwardRateFuntion1(t,tau_) + b2_ * forwardRateFuntion2(t,tau_);
}

double NelsonSiegel::getDerivativeInstantaneousForwardRate(double t) const {
    return forwardRateFuntion1(t, tau_)*(b2_*(1-t/tau_) - b1_);
}

double Svensson::getRate(double t) const {
    return b0_ + b1_ * rateFuntion1(t,tau1_) + b2_ * rateFuntion2(t,tau1_) + b3_ * rateFuntion2(t,tau2_);
}

double Svensson::getInstantaneousForwardRate(double t) const {
    return b0_ + b1_ * forwardRateFuntion1(t,tau1_) + b2_ * forwardRateFuntion2(t,tau1_) + b3_ * forwardRateFuntion2(t,tau2_);
}

double Svensson::getDerivativeInstantaneousForwardRate(double t) const {

    return forwardRateFuntion1(t, tau1_)*(b2_*(1-t/tau1_) - b1_) + b3_*forwardRateFuntion1(t, tau2_)*(1-t/tau2_)/tau2_;
}
