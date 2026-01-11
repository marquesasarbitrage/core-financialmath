#include "../../include/core-financialmath/calibration/nelsonsiegelsvensson.hpp"

NelsonSiegelCalibration::NelsonSiegelCalibration(const std::map<double, double>& data, bool isSpotRate):data_(data), isSpotRate_(isSpotRate), gridSize_(.5){};

std::shared_ptr<NelsonSiegelFamily> NelsonSiegelCalibration::fitOLS(double tau1, double tau2, bool useSvensson) const {
    
    std::vector<double> y; 
    std::vector<std::vector<double>> x;
    for (const auto& d : data_) {

        y.push_back(d.second);
        std::vector<double> X = {};
        if (isSpotRate_){
            X.push_back(NelsonSiegel::rateFuntion1(d.first,tau1));
            X.push_back(NelsonSiegel::rateFuntion2(d.first,tau1));
            if (useSvensson) X.push_back(NelsonSiegel::rateFuntion2(d.first,tau2));
        }else {
            X.push_back(NelsonSiegel::forwardRateFuntion1(d.first,tau1));
            X.push_back(NelsonSiegel::forwardRateFuntion2(d.first,tau1));
            if (useSvensson)X.push_back(NelsonSiegel::forwardRateFuntion2(d.first,tau2));
        }
        x.push_back(X);
        
    }
    OrdinaryLeastSquare ols(y,x,true);
    double beta0 = ols.getIntercept(); 
    std::vector<double> betas = ols.getCoefficients();
    std::shared_ptr<NelsonSiegelFamily> nss = nullptr;
    if (useSvensson) nss = std::make_shared<Svensson>(beta0,betas[0],betas[1],betas[2],tau1,tau2);
    else nss = std::make_shared<NelsonSiegel>(beta0,betas[0],betas[1],tau1);
    return nss;
}

Residuals NelsonSiegelCalibration::getLoss(double tau1, double tau2, bool useSvensson) const {

    std::shared_ptr<NelsonSiegelFamily> nss = fitOLS(tau1,tau2, useSvensson);
    std::vector<double> estimates; 
    std::vector<double> trueValues;
    for (const auto& d: data_)
    {
        estimates.push_back(isSpotRate_ ? nss->getRate(d.first) : nss->getInstantaneousForwardRate(d.first));
        trueValues.push_back(d.second);
    }
    return Residuals(estimates,trueValues);
}

double NelsonSiegelCalibration::getNelsonSiegelIniatialTau() const {

    double tau, meanSquaredError, tauWinner, minMSE;
    double tMax = std::prev(data_.end())->first;
    while (tau<=tMax)
    {
        tau += gridSize_; 
        meanSquaredError = getLoss(tau,10.0,false).getMSE();
        if (tau == gridSize_){
            tauWinner = tau;
            minMSE = meanSquaredError;
        }else{
            if (meanSquaredError<minMSE){
                tauWinner = tau; 
                minMSE = meanSquaredError;
            }
        }
    }
    return tauWinner;
}

std::vector<double> NelsonSiegelCalibration::getSvenssonIniatialTau() const {

    double tau1, tau2, meanSquaredError, tauWinner1, tauWinner2, minMSE;
    double tMax = std::prev(data_.end())->first;
    while (tau1<=tMax)
    {
        tau1 += gridSize_;
        while (tau2<=tMax)
        {
            tau2 += gridSize_; 
            meanSquaredError = getLoss(tau1,tau2,true).getMSE();
            if (tau1 == gridSize_ and tau2 ==gridSize_){
                tauWinner1 = tau1;
                tauWinner2 = tau2;
                minMSE = meanSquaredError;
            }else{
                if (meanSquaredError<minMSE){
                    tauWinner1 = tau1;
                    tauWinner2 = tau2; 
                    minMSE = meanSquaredError;
                }
            }
        }
            
    }
    return {tauWinner1,tauWinner2};
}

void NelsonSiegelCalibration::setGridSize(double value) { gridSize_ = value; }

std::shared_ptr<NelsonSiegelFamily> NelsonSiegelCalibration::fitNelsonSiegel() const {

    std::function<double(std::vector<double>)> targetFunction = [*this](std::vector<double> params)
    { 
        if (params[0]==0.0) return 1e10;
        Residuals loss = getLoss(params[0],10.0, false);
        return loss.getMSE();
    };

    std::vector<double> initialTau = {getNelsonSiegelIniatialTau()};
    NelderMead nm = NelderMead(initialTau,targetFunction);
    nm.setInitSimplexMethod(NelderMead::InitSimplexMethod::SYMMETRIC); 
    nm.setPerturbationParam(gridSize_);
    nm.optimize();
    if (!nm.getError()) return fitOLS(nm.getResult()[0], 10.0, false);
    else return fitOLS(initialTau[0], 10.0, false);
}

std::shared_ptr<NelsonSiegelFamily>  NelsonSiegelCalibration::fitSvensson() const {
    
    std::function<double(std::vector<double>)> targetFunction = [*this](std::vector<double> params) {

        if (params[0]==0.0 or params[1]==0.0) return 1e10;
        Residuals loss = getLoss(params[0], params[1], true);
        return loss.getMSE();
    };
    std::vector<double> params0 = getSvenssonIniatialTau(); 
    NelderMead nm = NelderMead(params0,targetFunction);
    nm.optimize();
    nm.setInitSimplexMethod(NelderMead::InitSimplexMethod::SYMMETRIC); 
    nm.setPerturbationParam(gridSize_);
    if (!nm.getError()) return fitOLS(nm.getResult()[0], nm.getResult()[1], true);
    else return fitOLS(params0[0], params0[1], true);
}