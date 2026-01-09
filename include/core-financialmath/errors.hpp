#pragma once 
#include <iostream>
#include <exception>

class FinMathLibraryError: public std::exception 
{
    public:
        const char* what() const noexcept override;
        explicit FinMathLibraryError(){};
        virtual ~FinMathLibraryError() = default;
    protected: 
        virtual std::string getErrorMessage() const = 0; 
    private:
        mutable std::string cachedMessage_;  // must be mutable to modify in const what()
};

namespace FinMathErrorRegistry
{
    namespace LetsBeRational {

        class InvalidNormalizedPriceRangeError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };

    }

    namespace TimeSerie {

        class InvalidSizeError final: public FinMathLibraryError {
            protected: 
                std::string getErrorMessage() const override; 
        };

        class StartDateAfterEndDateError final: public FinMathLibraryError {
            protected: 
                std::string getErrorMessage() const override; 
        };

        class InvalidReferenceDateError final: public FinMathLibraryError {
            protected: 
                std::string getErrorMessage() const override; 
        };
    }
    
};