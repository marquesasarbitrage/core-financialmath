#pragma once 
#include <iostream>
#include "errors/main.hpp"

namespace FinmathToolbox {

    void checkYearFraction(double value); 
    void checkStrikePrice(double value); 
    void checkSpotPrice(double value); 
    void checkImpliedVolatility(double value); 

}