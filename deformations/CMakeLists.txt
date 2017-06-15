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
ADD_EXECUTABLE(3dVolViewer 3dVolViewer)
TARGET_LINK_LIBRARIES(3dVolViewer ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(3dVolProjector 3dVolProjector)
TARGET_LINK_LIBRARIES(3dVolProjector ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(generate3d generate3d)
TARGET_LINK_LIBRARIES(generate3d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

#ADD_EXECUTABLE(testLocalDeformation3d testLocalDeformation3d)
#TARGET_LINK_LIBRARIES(testLocalDeformation3d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

#ADD_EXECUTABLE(testLocalDeformation2d testLocalDeformation2d)
#TARGET_LINK_LIBRARIES(testLocalDeformation2d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(testDisk testDisk)
TARGET_LINK_LIBRARIES(testDisk ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(deformation2d deformation2d)
TARGET_LINK_LIBRARIES(deformation2d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(deformation3d deformation3d)
TARGET_LINK_LIBRARIES(deformation3d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(drip2d drip2d)
TARGET_LINK_LIBRARIES(drip2d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(drip3d drip3d)
TARGET_LINK_LIBRARIES(drip3d ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(imageBlurring imageBlurring)
TARGET_LINK_LIBRARIES(imageBlurring ${DGTAL_LIBRARIES} ${Boost_LIBRARIES} ${FFTW})

ADD_EXECUTABLE(BENCH_ApproximatedMultiImage BENCH_ApproximatedMultiImage)
TARGET_LINK_LIBRARIES(BENCH_ApproximatedMultiImage ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})

ADD_EXECUTABLE(extractCells extractCells)
TARGET_LINK_LIBRARIES(extractCells ${DGTAL_LIBRARIES} ${Boost_LIBRARIES})
