# CMakeList.txt : CMake project for QM_RBM, include source and define
# project specific logic here.
#
cmake_minimum_required (VERSION 3.8)

project ("QM_RBM")

# TODO: Add tests and install targets if needed.

# Set CXX standard to C++11
set(CMAKE_CXX_STANDARD 11)

# Check if the user has specified a build type via the command line, e.g., (in the
# build-directory)
#
#   $cmake .. -DCMAKE_BUILD_TYPE=Release
#
# if not, set Debug as the default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()



# Add all compiler warnings. Remember, a warning in C++ can very often lead to a bug
# later on so _fix your warnings!_
############################################################################################set(CMAKE_CXX_FLAGS "-Wall -Wextra")
# Add debugging flag for tools such as gdb and valgrind.
set(CMAKE_CXX_FLAGS_DEBUG "-g")
# Add full compiler optimization when in the release build type. Other options include
# -O0, -O1, -O2.
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# Add directories containing header-files
include_directories(
    .
    Hamiltonians
    InitialStates
    Math
    WaveFunctions
)

# Add sources from all directories
file(
    GLOB SOURCES "*.cpp"
    GLOB SOURCES "Hamiltonians/*.cpp"
    GLOB SOURCES "InitialStates/*.cpp"
    GLOB SOURCES "Math/*.cpp"
    GLOB SOURCES "WaveFunctions/*.cpp"
)


# Create executable called "vmc"
add_executable(QM_RBM ${SOURCES} "Benchmark.h" "runtype.h" "sampler_gd.h" "sampler_gd.cpp" "Hamiltonians/ellipticOscillator.h" "Hamiltonians/ellipticOscillator.cpp" "WaveFunctions/correlatedGaussian.h" "WaveFunctions/correlatedGaussian.cpp")
