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

    namespace LetsBeRational {

        class InvalidNormalizedPriceRangeError final: public FinMathLibraryError
        {
            protected: 
                std::string getErrorMessage() const override; 
        };

    }

    

}