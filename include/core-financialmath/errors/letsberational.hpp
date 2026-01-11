#pragma once 
#include <iostream>
#include "base.hpp"

namespace FinMathErrorRegistry {

    namespace LetsBeRational {

        class InvalidNormalizedPriceRangeError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };

        class OptionHasNoTimeValueError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };

        class NegativeOptionPriceError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };
    }
}
