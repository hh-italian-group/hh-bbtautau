project(hh-bbtautau)
cmake_minimum_required(VERSION 3.6)

get_filename_component(AnalysisTools_DIR "${PROJECT_SOURCE_DIR}/../AnalysisTools" ABSOLUTE)
include("${AnalysisTools_DIR}/cmake/include/common.cmake")

file(GLOB_RECURSE ANA_TOOLS_SRC "${AnalysisTools_DIR}/*.cpp")
add_library(AnalysisTools STATIC ${ANA_TOOLS_SRC})
add_dependencies(AnalysisTools GenerateRootDict)

get_filename_component(h_tautau_DIR "${PROJECT_SOURCE_DIR}/../h-tautau" ABSOLUTE)
file(GLOB_RECURSE H_TAUTAU_SRC "${h_tautau_DIR}/*.cpp")
add_library(h-tautau STATIC ${H_TAUTAU_SRC})

file(GLOB_RECURSE HTT_UTILITIES_SRC "${CMSSW_BASE_SRC}/HTT-utilities/*.cc")
add_library(HTT-utilities STATIC ${HTT_UTILITIES_SRC})
set_source_files_properties(${HTT_UTILITIES_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE TauAnalysisTool_SRC "${CMSSW_BASE_SRC}/TauAnalysisTools/TauTriggerSFs/src/*.cc")
add_library(TauAnalysisTools STATIC ${TauAnalysisTool_SRC})

file(GLOB_RECURSE VBFTriggerSFs_SRC "${CMSSW_BASE_SRC}/VBFTriggerSFs/*.cc")
add_library(VBFTriggerSFs STATIC ${VBFTriggerSFs_SRC})
target_include_directories(VBFTriggerSFs PRIVATE "${CMSSW_BASE_SRC}/VBFTriggerSFs/VBFTriggerSFs/interface")
set_source_files_properties(${VBFTriggerSFs_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE TauPOG_SRC "${CMSSW_BASE_SRC}/TauPOG/TauIDSFs/src/*.cc")
add_library(TauPOG STATIC ${TauPOG_SRC})
set_source_files_properties(${TauPOG_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE HH_KINFIT2_SRC "${CMSSW_BASE_SRC}/HHKinFit2/*.cpp")
add_library(HHKinFit2 STATIC ${HH_KINFIT2_SRC})
target_include_directories(HHKinFit2 PRIVATE "${CMSSW_BASE_SRC}/HHKinFit2/HHKinFit2/interface")
set_source_files_properties(${HH_KINFIT2_SRC} PROPERTIES COMPILE_FLAGS "-w")

file(GLOB_RECURSE SVFIT_SRC "${CMSSW_BASE_SRC}/TauAnalysis/ClassicSVfit/src/*.cc")
add_library(SVfit STATIC ${SVFIT_SRC})
target_compile_definitions(SVfit PRIVATE USE_SVFITTF=1)
set_source_files_properties(${SVFIT_SRC} PROPERTIES COMPILE_FLAGS "-w")

#include("${AnalysisTools_DIR}/cmake/include/tensorflow.cmake")
if(${TF_AVAILABLE})
    file(GLOB_RECURSE MulticlassInference_SRC "${CMSSW_BASE_SRC}/MulticlassInference/MulticlassInference/src/*.cc")
    add_library(MulticlassInference STATIC ${MulticlassInference_SRC})
    target_link_libraries(MulticlassInference ${TF_LIBRARIES})
    target_include_directories(MulticlassInference PRIVATE ${TF_INCLUDES})
    set_source_files_properties(${MulticlassInference_SRC} PROPERTIES COMPILE_FLAGS "-w")

    file(GLOB_RECURSE cms_hh_tf_inference_SRC "${CMSSW_BASE_SRC}/cms_hh_tf_inference/inference/src/*.cc")
    add_library(cms_hh_tf_inference STATIC ${cms_hh_tf_inference_SRC})
    target_link_libraries(cms_hh_tf_inference ${TF_LIBRARIES})
    target_include_directories(cms_hh_tf_inference PRIVATE ${TF_INCLUDES})
    set_source_files_properties(${cms_hh_tf_inference_SRC} PROPERTIES COMPILE_FLAGS "-w")

    file(GLOB_RECURSE cms_hh_proc_interface_SRC "${CMSSW_BASE_SRC}/cms_hh_proc_interface/processing/src/*.cc")
    add_library(cms_hh_proc_interface STATIC ${cms_hh_proc_interface_SRC})
    set_source_files_properties(${cms_hh_proc_interface_SRC} PROPERTIES COMPILE_FLAGS "-w")

    target_link_libraries(CalcMulticlassDNN MulticlassInference ${TF_LIBRARIES})
    target_include_directories(CalcMulticlassDNN PRIVATE ${TF_INCLUDES})
    set_source_files_properties("${CMSSW_BASE_SRC}/hh-bbtautau/Analysis/source/CalcMulticlassDNN.cxx" PROPERTIES COMPILE_FLAGS "-w")

    target_link_libraries(CalcDNN cms_hh_tf_inference cms_hh_proc_interface ${TF_LIBRARIES})
    target_include_directories(CalcDNN PRIVATE ${TF_INCLUDES})
    set_source_files_properties("${CMSSW_BASE_SRC}/hh-bbtautau/Analysis/source/CalcDNN.cxx" PROPERTIES COMPILE_FLAGS "-Wno-implicit-float-conversion")
else()
    set_target_properties(CalcDNN CalcMulticlassDNN PROPERTIES EXCLUDE_FROM_ALL 1)
endif()

add_library("hh-bbtautau" STATIC ${SOURCE_LIST})

foreach(exe_name ${EXE_LIST})
    target_link_libraries("${exe_name}" hh-bbtautau h-tautau AnalysisTools HTT-utilities SVfit HHKinFit2 TauAnalysisTools VBFTriggerSFs TauPOG)
endforeach()

set_target_properties(LimitConfigurationProducer BjetSelectionStudy BjetEffSF KinFitStudy SMWeight_t PROPERTIES EXCLUDE_FROM_ALL 1)
