project(hh-bbtautau)
cmake_minimum_required(VERSION 2.8)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

add_library("hh-bbtautau" OBJECT ${HEADER_LIST})

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc" "${CMSSW_BASE_SRC}/HTT-utilities/*.h")
add_library(HTT-utilities STATIC "${HTT_UTILITIES_SRC}")

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" HTT-utilities)
endforeach()

set_target_properties(LimitConfigurationProducer bbetauAnalyzer bbtautauAnalyzer BjetSelectionStudy PROPERTIES EXCLUDE_FROM_ALL 1)
