project(hh-bbtautau)
cmake_minimum_required(VERSION 3.6)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE HH_KINFIT2_SRC "${CMSSW_BASE_SRC}/HHKinFit2/*.cpp")
add_library(HHKinFit2 STATIC ${HH_KINFIT2_SRC})
target_include_directories(HHKinFit2 PRIVATE "${CMSSW_BASE_SRC}/HHKinFit2/HHKinFit2/interface")
set_source_files_properties(${HH_KINFIT2_SRC} PROPERTIES COMPILE_FLAGS "-w")

#[[
set(PrecompHeader "${PROJECT_SOURCE_DIR}/Run/include/Precompiled.h")
set(PrecompHeaderOut "${PROJECT_SOURCE_DIR}/Run/include/Precompiled.h.pch")
string (REPLACE " " ";" PchFlags ${CMAKE_CXX_FLAGS})
get_property(PchIncludeRaw DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
set(PchInclude "")
foreach(dir ${PchIncludeRaw})
    set(PchInclude "${PchInclude};-I${dir}")
endforeach()
add_custom_command(OUTPUT "${PrecompHeaderOut}"
                   COMMAND ${CMAKE_CXX_COMPILER} -x c++-header -w ${PchFlags} ${PchInclude} "${PrecompHeader}"  -o "${PrecompHeaderOut}"
                   IMPLICIT_DEPENDS CXX "${PrecompHeader}"
                   VERBATIM)
add_custom_target(GeneratePch DEPENDS "${PrecompHeaderOut}")
]]
foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" HTT-utilities HHKinFit2)
#    add_dependencies("${exe_name}" GeneratePch)
#    target_compile_options("${exe_name}" BEFORE PRIVATE "-include;${PrecompHeader}")
endforeach()

set_target_properties(LimitConfigurationProducer BjetSelectionStudy BjetEffSF KinFitStudy SMWeight_t PROPERTIES EXCLUDE_FROM_ALL 1)

add_library("hh-bbtautau" OBJECT ${HEADER_LIST} ${SOURCE_LIST} ${EXE_SOURCE_LIST} ${SCRIPT_LIST} ${CONFIG_LIST})
set_target_properties("hh-bbtautau" PROPERTIES EXCLUDE_FROM_ALL 1)
