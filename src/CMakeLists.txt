# SPDX-FileCopyrightText: 2021 Philipp Basler
#
# SPDX-License-Identifier: GPL-3.0-or-later

add_executable(C2HDM C2HDM.cpp)
target_link_libraries(C2HDM BSMPT::Minimizer BSMPT::Models )
target_compile_features(C2HDM PUBLIC cxx_std_14)


add_executable(R2HDM R2HDM.cpp)
target_link_libraries(R2HDM BSMPT::Minimizer BSMPT::Models )
target_compile_features(R2HDM PUBLIC cxx_std_14)

add_executable(RN2HDM RN2HDM.cpp)
target_link_libraries(RN2HDM BSMPT::Minimizer BSMPT::Models )
target_compile_features(RN2HDM PUBLIC cxx_std_14)
