project(hh-bbtautau)
cmake_minimum_required(VERSION 3.0)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE HH_KINFIT2_SRC "${CMSSW_BASE_SRC}/HHKinFit2/*.cpp")
add_library(HHKinFit2 STATIC ${HH_KINFIT2_SRC})
target_include_directories(HHKinFit2 PRIVATE "${CMSSW_BASE_SRC}/HHKinFit2/HHKinFit2/interface")
set_source_files_properties(${HH_KINFIT2_SRC} PROPERTIES COMPILE_FLAGS "-w")

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" HTT-utilities HHKinFit2)
endforeach()

set_target_properties(LimitConfigurationProducer bbetauAnalyzer bbmutauAnalyzer bbtautauAnalyzer BjetSelectionStudy BjetEffSF KinFitStudy PROPERTIES EXCLUDE_FROM_ALL 1)

add_library("hh-bbtautau" OBJECT ${HEADER_LIST} ${SOURCE_LIST} ${EXE_SOURCE_LIST} ${SCRIPT_LIST} ${CONFIG_LIST})
set_target_properties("hh-bbtautau" PROPERTIES EXCLUDE_FROM_ALL 1)
