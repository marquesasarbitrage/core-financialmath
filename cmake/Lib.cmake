add_library(
    core-financialmath 
    STATIC 
        src/errors/base.cpp
        src/nelsonsiegelsvensson.cpp
        src/letsberational.cpp
        src/baw.cpp
        src/lewis.cpp)
target_link_libraries(core-financialmath PUBLIC core-datetime)
target_link_libraries(core-financialmath PUBLIC core-math)
target_include_directories(core-financialmath  PUBLIC include)
