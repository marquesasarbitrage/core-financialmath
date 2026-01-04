#include <iostream>
#include <cassert>
#include "../include/core-financialmath/lewis.hpp"

bool isClose(double a, double b, double tol = 1e-3) {
    return std::abs(a - b) <= tol;
}

void testBlackScholes()
{
    // Input parameters
    double S = 100.0;    // Stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.02;     // Dividend yield
    double sigma = 0.20; // Volatility
    double T = 2.0;      // Time to maturity

    // Option prices and Greeks
    double callPrice = 13.5218;
    double putPrice = 7.9266;

    double x = log(S*exp((r-q)*T)/K); 
    double normalizedSigma = sigma*sqrt(T);

    LewisEuropeanVanillaPrice::BlackScholes bs(sigma,T,64);
    double lewisPutPrice = bs.getPrice(S*exp((r-q)*T), K, false);
    double lewisCallPrice = bs.getPrice(S*exp((r-q)*T), K, true);

    assert(isClose(lewisCallPrice*std::exp(-r*T), callPrice, 1e-2));
    assert(isClose(lewisPutPrice*std::exp(-r*T), putPrice, 1e-2));

    std::cout << "Tests for the Lewis Black-Scholes object for european vanilla pricing have been passed! "<<std::endl;

}

void testHeston() {

    double S = 100.0;    // Stock price
    double K = 100.0;    // Strike price
    double r = 0.05;     // Risk-free rate
    double q = 0.01;     // Dividend yield
    double T = 1.5;      // Time to maturity
    double F = S * std::exp((r - q) * T); // Forward price
    
    LewisEuropeanVanillaPrice::Heston heston(2.0,0.05,0.3,0.45,0.05,T,64);
    double expectedHestonPrice = 13.2561; 
    double modelHestonPrice = heston.getPrice(F,K,true)*exp(-r*T);
    assert(isClose(modelHestonPrice, expectedHestonPrice, 1e-4));

    std::cout << "Tests for the Lewis Heston object for european vanilla pricing have been passed! \n";

}

int main() {

    testBlackScholes();
    testHeston();
    return 0; 
}