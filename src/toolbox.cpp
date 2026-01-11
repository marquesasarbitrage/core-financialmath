#include "../include/core-financialmath/toolbox.hpp"

namespace FinmathToolbox {

    void checkYearFraction(double value) {

        if (value <= 0.0) throw FinMathErrorRegistry::NonPositiveYearFractionError();
    }

    void checkStrikePrice(double value) {

        if (value <= 0.0) throw FinMathErrorRegistry::NonPositiveStrikePriceError();
    }

    void checkSpotPrice(double value) {

        if (value <= 0.0) throw FinMathErrorRegistry::NonPositiveStrikePriceError();
    }

    void checkImpliedVolatility(double value) {

        if (value <= 0.0) throw FinMathErrorRegistry::NonPositiveImpliedVolatilityError();
    }

}