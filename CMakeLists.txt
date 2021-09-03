# SPDX-FileCopyrightText: 2021 Philipp Basler
#
# SPDX-License-Identifier: GPL-3.0-or-later

cmake_minimum_required(VERSION 3.5)

project(
  BSMPT-UnitTestReferenceCreator
  VERSION 1.0.0
  LANGUAGES C CXX)
  
  # -------------------------------- setup cmake --------------------------------
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW) # use <packagename>_ROOT as search path
endif(POLICY CMP0074)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/bin")

find_package(PkgConfig)
find_package(BSMPT 2.2.1 REQUIRED)

add_subdirectory(src)
