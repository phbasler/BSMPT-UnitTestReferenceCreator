add_executable(C2HDM C2HDM.cpp)
target_link_libraries(C2HDM BSMPT::Minimizer BSMPT::Models )
target_compile_features(C2HDM PUBLIC cxx_std_14)
