#include "../../include/core-financialmath/errors/base.hpp"

namespace FinMathErrorRegistry
{
    const char* FinMathLibraryError::what() const noexcept {

        if (cachedMessage_.empty()) {
                cachedMessage_ = getErrorMessage();  
            }
        return cachedMessage_.c_str();
    }
    
    namespace LetsBeRational {

        std::string InvalidNormalizedPriceRangeError::getErrorMessage() const {return "The normalized price must be strictly superior to 0 and inferior to the b_max threshold.";}

    }

};