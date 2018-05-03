#Specify the version being used aswell as the language
cmake_minimum_required(VERSION 2.8)
#Name your project here
project(bimsim)



#set the module directory
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}")
set(CMAKE_AUTOMOC ON)

#add flags for gropf
#add_definitions(-pg)
#set(CMAKE_NVCC_FLAGS "-w -pg -g")
#set(CMAKE_CXX_FLAGS "-w -pg -g -fpermissive")
#set(CMAKE_C_FLAGS "-w -pg -g -fpermissive")

# As moc files are generated in the binary dir, tell CMake
  # to always look for includes there:
set(CMAKE_INCLUDE_CURRENT_DIR ON)

#set up CUDA
#set(CUDA_HOST_COMPILER /usr/bin/gcc-4.8)
find_package(CUDA REQUIRED)

#find OpenCV
find_package(OpenCV REQUIRED )

add_definitions(-DUSING_OPENCV)


#find the STIM library
find_package(STIM)

#if(CUDA_FOUND)
#	set(CUDA_NVCC_FLAGS -g; -G;-lineinfo)# -maxrregcount=16)
#endif()

#build the executable in the binary directory on MS Visual Studio
if ( MSVC )
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${OUTPUT_DIRECTORY}")
	SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${OUTPUT_DIRECTORY}")
endif ( MSVC )

if( CMAKE_COMPILER_IS_GNUCC )
	SET( CMAKE_CXX_FLAGS "--std=c++11")
	#SET( CMAKE_NVCC_FLAGS "--std=c++11")
endif( CMAKE_COMPILER_IS_GNUCC )

#set the include directories
include_directories(
	${CMAKE_CURRENT_BINARY_DIR}
        ${RTS_INCLUDE_DIR}
   	${Boost_INCLUDE_DIR}
	${CMAKE_CURRENT_SOURCE_DIR}
	${STIM_INCLUDE_DIRS}
)

#Assign source files to the appropriate variables
file(GLOB SRC_CPP "*.cpp")
file(GLOB SRC_H "*.h")
file(GLOB SRC_UI "*.ui")
file(GLOB SRC_CU "*.cu")


#create an executable
cuda_add_executable(bimsim 
                    ${SRC_CPP}
                    ${SRC_H}
                    ${SRC_CU}
		    ${SRC_RTS}
)


#set the link libraries
target_link_libraries(bimsim
			  ${CUDA_cufft_LIBRARY}
			  ${CUDA_cublas_LIBRARY}
			  ${OpenCV_LIBS}
		)



