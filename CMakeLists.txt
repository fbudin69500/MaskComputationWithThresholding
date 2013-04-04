CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
CMAKE_POLICY(VERSION 2.6)
PROJECT(MaskComputationWithThresholding)


#####################     ITK    ############################################

FIND_PACKAGE(ITK REQUIRED)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
   MESSAGE(FATAL_ERROR "ITK not found. Please set ITK_DIR.") 
ENDIF(ITK_FOUND)

#  INCLUDE_DIRECTORIES( 
#    ${ITK_SOURCE_DIR}/Code
#    )
#  SET(ITK_LIBRARIES
#    ITKCommon
#    ITKIO
#    ITKNumerics
#    ITKBasicFilters
#    ITKAlgorithms
#    )
#  LINK_LIBRARIES(
#    ${ITK_LIBRARIES}
#    )


find_package(GenerateCLP REQUIRED)
if(GenerateCLP_FOUND)
  include(${GenerateCLP_USE_FILE})
else(GenerateCLP_FOUND)
   message(FATAL_ERROR "GenerateCLP Not Found")
endif(GenerateCLP_FOUND)

set(Filter MaskComputationWithThresholding)

set(${Filter}_SOURCE ${Filter}.cxx)
GENERATECLP(${Filter}_SOURCE ${Filter}.xml)
add_executable( ${Filter} ${${Filter}_SOURCE})
target_link_libraries( ${Filter}  ${ITK_LIBRARIES})

