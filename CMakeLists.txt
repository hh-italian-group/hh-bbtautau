project("hh-bbtautau")
cmake_minimum_required(VERSION 2.8)

set(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${AnalysisTools_DIR}/cmake/modules")
include("${AnalysisTools_DIR}/common.cmake")

add_library("hh-bbtautau" OBJECT ${HEADER_LIST})
