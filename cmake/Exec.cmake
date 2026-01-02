add_executable(finmath-nss ${CMAKE_CURRENT_SOURCE_DIR}/tests/nelsonsiegelsvensson/nss.cpp)
target_link_libraries(finmath-nss PUBLIC core-financialmath)

add_executable(finmath-letsberational ${CMAKE_CURRENT_SOURCE_DIR}/tests/letsberational/lbr.cpp)
target_link_libraries(finmath-letsberational PUBLIC core-financialmath)