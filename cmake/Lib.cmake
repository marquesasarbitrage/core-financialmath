add_library(
    core-financialmath 
    STATIC 
        src/nelsonsiegelsvensson.cpp
        src/letsberational.cpp
        src/lewis.cpp
        src/errors.cpp)
target_link_libraries(core-financialmath PUBLIC core-datetime)
target_link_libraries(core-financialmath PUBLIC core-math)
target_include_directories(core-financialmath  PUBLIC include)
