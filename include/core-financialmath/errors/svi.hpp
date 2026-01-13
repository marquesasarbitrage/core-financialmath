#pragma once 
#include <iostream>
#include "base.hpp"

namespace FinMathErrorRegistry {

    namespace SVI {

        class InvalidSVIParamError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };

        class ButterflyArbitrageError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };


    }
}
