#pragma once 
#include <iostream>
#include <exception>

namespace FinMathErrorRegistry {

    class FinMathLibraryError: public std::exception {

        public:
            const char* what() const noexcept override;
            explicit FinMathLibraryError(){};
            virtual ~FinMathLibraryError() = default;
        protected: 
            virtual std::string getErrorMessage() const = 0; 
        private:
            mutable std::string cachedMessage_;  // must be mutable to modify in const what()
    };


    class NonPositiveYearFractionError final: public FinMathLibraryError
    {
        protected: 
            std::string getErrorMessage() const override; 
    };

    class NonPositiveStrikePriceError final: public FinMathLibraryError
    {
        protected: 
            std::string getErrorMessage() const override; 
    };

    class NonPositiveSpotPriceError final: public FinMathLibraryError
    {
        protected: 
            std::string getErrorMessage() const override; 
    };

    class NonPositiveImpliedVolatilityError final: public FinMathLibraryError
    {
        protected: 
            std::string getErrorMessage() const override; 
    };


}