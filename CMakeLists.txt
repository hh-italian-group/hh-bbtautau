project("hh-bbtautau")
cmake_minimum_required(VERSION 2.8)

include("$ENV{CMSSW_BASE_SRC}/AnalysisTools/common.cmake")

add_library("hh-bbtautau" OBJECT ${HEADER_LIST})
