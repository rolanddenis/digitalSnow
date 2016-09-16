PROJECT(Deformations)

#C++11
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

#Release
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -mtune=native")

#DGtal
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

FIND_PACKAGE(DGtal REQUIRED)
MESSAGE(STATUS "Found DGtal: ${DGTAL_LIBRARIES}")
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})


#Inclusion de Boost
include(FindBoost)
find_package(Boost 1.49.0 COMPONENTS program_options filesystem system timer REQUIRED )
message(STATUS "Found Boost: ${Boost_LIBRARIES} in ${Boost_LIBRARY_DIRS} and ${Boost_INCLUDE_DIRS} ")
#message(STATUS "More precisely: ${Boost_PROGRAM_OPTIONS_LIBRARY} ")
include_directories(${Boost_INCLUDE_DIRS})

#fftw
SET(FFTW fftw3 fftw3_threads m pthread)
INCLUDE_DIRECTORIES(/usr/include/)
LINK_DIRECTORIES(/usr/lib/)

##########################################
ADD_EXECUTABLE(Snow3D Snow3D)
TARGET_LINK_LIBRARIES(Snow3D ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW} )
