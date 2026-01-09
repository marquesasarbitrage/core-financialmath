#include "../include/core-financialmath/errors.hpp"

const char* FinMathLibraryError::what() const noexcept
{
    if (cachedMessage_.empty()) {
            cachedMessage_ = getErrorMessage();  
        }
    return cachedMessage_.c_str();
}

namespace FinMathErrorRegistry
{
    namespace LetsBeRational {

        std::string InvalidNormalizedPriceRangeError::getErrorMessage() const {return "The normalized price must be strictly superior to 0 and inferior to the b_max threshold.";}

    }

    namespace TimeSerie {

        std::string InvalidSizeError::getErrorMessage() const {return "A time serie must be represented by at least 2 points.";}
        std::string StartDateAfterEndDateError::getErrorMessage() const {return "When getting a sub time serie, the given start date must come before the end date.";}
        std::string InvalidReferenceDateError::getErrorMessage() const {return "The reference date must be within the min/max range observed within the time serie.";}
    }
    
};