add_executable(finmath-calibration-nss ${CMAKE_CURRENT_SOURCE_DIR}/tests/calibration/nelsonsiegelsvensson/nss.cpp)
target_link_libraries(finmath-calibration-nss PUBLIC core-financialmath)

add_executable(finmath-letsberational ${CMAKE_CURRENT_SOURCE_DIR}/tests/letsberational/lbr.cpp)
target_link_libraries(finmath-letsberational PUBLIC core-financialmath)

add_executable(finmath-lewis ${CMAKE_CURRENT_SOURCE_DIR}/tests/lewis.cpp)
target_link_libraries(finmath-lewis PUBLIC core-financialmath)

add_executable(finmath-baw ${CMAKE_CURRENT_SOURCE_DIR}/tests/baw.cpp)
target_link_libraries(finmath-baw PUBLIC core-financialmath)
