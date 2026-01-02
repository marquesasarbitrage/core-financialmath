#include <iostream>
#include <cassert>
#include <iomanip> 
#include <filesystem>
#include <fstream>
#include "../include/core-financialmath/letsberational.hpp"

bool isClose(double a, double b, double tol = 1e-3) {

    return std::abs(a - b) <= tol;
}

void testOptionPrice()
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

    double modelCallPrice = sqrt(S*exp((r-q)*T)*K)*LetsBeRational::getCallPrice(x,T,sigma)*exp(-r*T);
    double modelPutPrice = sqrt(S*exp((r-q)*T)*K)*LetsBeRational::getPutPrice(x,T,sigma)*exp(-r*T);
    LetsBeRational::ImpliedVolatilityResult impliedVolPut = LetsBeRational::getPutImpliedVolatility(x,T,LetsBeRational::getPutPrice(x,T,sigma));
    LetsBeRational::ImpliedVolatilityResult impliedVolCall = LetsBeRational::getCallImpliedVolatility(x,T,LetsBeRational::getCallPrice(x,T,sigma));
    assert(isClose(modelCallPrice, callPrice, 1e-4));
    assert(isClose(modelPutPrice, putPrice, 1e-4));
    assert(isClose(impliedVolPut.value_, .2, 1e-3));
    assert(isClose(impliedVolCall.value_, .2, 1e-3));

}

void writeNumericalResultLetsBeRational(double sigmaMin, double sigmaMax, double xMin, double xMax, std::string name) {

    std::filesystem::create_directory("LetsBeRationalResult");
    std::ofstream file("LetsBeRationalResult/" + name +".csv");
    if (!file.is_open()) {
        std::cerr << "Error: could not open file." << std::endl;
    }

    file << "logMoneyness,beta,normalizedSigma,estimatedSigma,relativeDiffNormalizedSigma,iterations,errorMessage \n";

    double sigmaDelta = (sigmaMax - sigmaMin)/100.0; 
    double xDelta = (xMax - xMin)/100.0; 
    double sigma = sigmaMin, x = xMin;
    std::vector<double> sigmas, moneyness; 
    double price, estSigma;
    for (int i = 0; i<100; i++) {
        sigmas.push_back(sigma);
        moneyness.push_back(x);
        x += xDelta; 
        sigma += sigmaDelta;
    }
    
    for (double s : sigmas) {

        for (double xx: moneyness) {

            price = LetsBeRational::_getPrice(xx,s,true);
            LetsBeRational::ImpliedVolatilityResult result = LetsBeRational::getNewtonNormalizedVolatility(price,xx,true);
            estSigma = result.value_; 
            std::string errorMessage;
            if (result.error_ != nullptr) {
                try {std::rethrow_exception(result.error_);}
                catch (const std::exception& e) {errorMessage = e.what();}
            } else { errorMessage = "None"; }
            file << xx << "," << price << "," << s << "," << estSigma << "," << std::abs((estSigma-s)/s)<< "," << result.iterations_ <<"," << errorMessage <<"," << "\n";
        }
    }

    file.close();
}

void testReproduceNumericalResultFromPaper() {

    // Figure 7, 8 and 9 from the original paper Let's Be Rational
    writeNumericalResultLetsBeRational(10e-7,1.22,0,3, "Figure7");
    writeNumericalResultLetsBeRational(10e-5,7.07,0,16, "Figure9");
    writeNumericalResultLetsBeRational(10e-5,0.18,0,10e-5, "Figure8");

}

int main() {

    testOptionPrice(); 
    testReproduceNumericalResultFromPaper();
    std::cout << "All tests for the Lets Be Rational object has been passed !"<<std::endl;
    return 0; 
}