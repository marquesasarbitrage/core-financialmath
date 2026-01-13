add_library(
    core-financialmath 
    STATIC 
        src/errors/base.cpp
        src/errors/letsberational.cpp
        src/errors/svi.cpp

        src/toolbox.cpp

        src/calibration/nss.cpp

        src/letsberational.cpp
        src/baw.cpp
        src/svi.cpp
        src/lewis.cpp
        src/nss.cpp
)
target_link_libraries(core-financialmath PUBLIC core-datetime)
target_link_libraries(core-financialmath PUBLIC core-math)
target_include_directories(core-financialmath  PUBLIC include)
