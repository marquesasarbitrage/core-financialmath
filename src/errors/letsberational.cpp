#include "../../include/core-financialmath/errors/letsberational.hpp"

namespace FinMathErrorRegistry
{

    namespace LetsBeRational {

        std::string InvalidNormalizedPriceRangeError::getErrorMessage() const {return "The normalized price must be strictly superior to 0 and inferior to the b_max threshold.";}
        std::string OptionHasNoTimeValueError::getErrorMessage() const {return "The option has no time value remaining.";}
        std::string NegativeOptionPriceError::getErrorMessage() const {return "Option price cannot be negative.";}


    }

};