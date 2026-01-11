#include "../../include/core-financialmath/errors/base.hpp"

namespace FinMathErrorRegistry
{
    const char* FinMathLibraryError::what() const noexcept {

        if (cachedMessage_.empty()) {
                cachedMessage_ = getErrorMessage();  
            }
        return cachedMessage_.c_str();
    }
    
    std::string NonPositiveYearFractionError::getErrorMessage() const {return "Year fraction must be non-positve to be valid.";}
    std::string NonPositiveStrikePriceError::getErrorMessage() const {return "Strike price must be non-positve to be valid.";}
    std::string NonPositiveSpotPriceError::getErrorMessage() const {return "Spot price must be non-positve to be valid.";}
    std::string NonPositiveImpliedVolatilityError::getErrorMessage() const {return "Implied volatility must be non-positve to be valid.";}

};