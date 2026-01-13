#include "../../include/core-financialmath/errors/svi.hpp"

namespace FinMathErrorRegistry
{

    namespace SVI {

        std::string InvalidSVIParamError::getErrorMessage() const {return "The SVI slice cannot be constructed because some parameters are incorrect.";}
        std::string ButterflyArbitrageError::getErrorMessage() const {return "The SVI slice is not free of butterfly arbitrage.";}
        


    }

};