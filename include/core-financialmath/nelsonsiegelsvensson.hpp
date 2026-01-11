#pragma once 
#include <iostream>
#include <map>
#include "core-math/statistics/regression.hpp"
#include "core-math/interpolation.hpp"
#include "core-math/optim/neldermead.hpp"
#include "core-math/statistics/residuals.hpp"

// Model refrences 
// Parsimonious Modelling of Yield Curve - Nelson and Siegel (1987) : https://www.jstor.org/stable/2352957
// Estimating Forward Interest Rates with the Extended Nelson & Siegel Method - Svensson (1995) : https://larseosvensson.se/papers/95QRabs/

class NelsonSiegelFamily
{
    public:
        NelsonSiegelFamily(){}; 
        virtual ~NelsonSiegelFamily() = default;

        virtual double getRate(double t) const = 0; 
        virtual double getInstantaneousForwardRate(double t) const = 0;
        virtual double getDerivativeInstantaneousForwardRate(double t) const = 0;
        virtual double getBeta0() const = 0;
        virtual double getBeta1() const = 0;

        static double rateFuntion1(double t, double tau);
        static double rateFuntion2(double t, double tau);
        static double forwardRateFuntion1(double t, double tau);
        static double forwardRateFuntion2(double t, double tau);

};

class NelsonSiegel final: public NelsonSiegelFamily
{
    public: 
        NelsonSiegel(double b0, double b1, double b2, double tau); 
        ~NelsonSiegel() = default; 

        double getRate(double t) const override; 
        double getInstantaneousForwardRate(double t) const override;
        double getDerivativeInstantaneousForwardRate(double t) const override;

        double getBeta1() const override; 
        double getBeta2() const; 
        double getBeta0() const override; 
        double getTau() const; 

        void setBeta0(double b0);
        void setBeta1(double b1);
        void setBeta2(double b2);
        void setTau(double tau);
    
    private:
        double b0_;  
        double b1_;  
        double b2_;  
        double tau_;
        
}; 

class Svensson final: public NelsonSiegelFamily
{
    public: 
        Svensson(double b0, double b1, double b2, double b3, double tau1, double tau2); 
        ~Svensson() = default;

        double getRate(double t) const override; 
        double getInstantaneousForwardRate(double t) const override;
        double getDerivativeInstantaneousForwardRate(double t) const override;

        double getBeta1() const override; 
        double getBeta2() const; 
        double getBeta0() const override;  
        double getBeta3() const;
        double getTau1() const;
        double getTau2() const;

        void setBeta0(double b0);
        void setBeta1(double b1);
        void setBeta2(double b2);
        void setTau1(double tau1);
        void setBeta3(double b3);
        void setTau2(double tau2);

    private:
        double b0_;  
        double b1_;  
        double b2_;  
        double b3_; 
        double tau1_;
        double tau2_;
};

