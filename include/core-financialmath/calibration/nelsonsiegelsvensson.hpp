#pragma once 
#include <iostream>
#include "core-math/statistics/regression.hpp"
#include "core-math/interpolation.hpp"
#include "core-math/optim/neldermead.hpp"
#include "core-math/statistics/residuals.hpp"
#include "../../../include/core-financialmath/nelsonsiegelsvensson.hpp"

class NelsonSiegelCalibration {

    public:
        NelsonSiegelCalibration(const std::map<double, double>& data, bool isSpotRate);
        ~NelsonSiegelCalibration() = default;

        
        std::shared_ptr<NelsonSiegelFamily> fitOLS(double tau1, double tau2, bool useSvensson) const; 
        std::shared_ptr<NelsonSiegelFamily> fitNelsonSiegel() const; 
        std::shared_ptr<NelsonSiegelFamily> fitSvensson() const; 
        void setGridSize(double value); 

    private: 
        std::map<double, double> data_;
        bool isSpotRate_;
        double gridSize_;

        Residuals getLoss(double tau1, double tau2, bool useSvensson) const;
        double getNelsonSiegelIniatialTau() const;
        std::vector<double>  getSvenssonIniatialTau() const;

};



