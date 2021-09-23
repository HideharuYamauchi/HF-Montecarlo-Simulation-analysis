# CMakeLists.txt for the "event" package. It creates a library and a main program.
# If ROOT is not installed in a default system location you need to tell CMake where to find it.
# Sourcing `thisroot.sh` already sets the required environment variables.
# Otherwise, you must tell the build system where to look for ROOT,
# for example by passing `-DROOT_DIR="/path/to/root/installation` at CMake configuration time.
   cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
   project(RF)
# Locate the ROOT package and define a number of useful targets and variables.
   find_package(ROOT REQUIRED)
   #include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
   include_directories(${PROJECT_SOURCE_DIR}/include)
   file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
   file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
# Create a shared library.
# Passing cmake targets such as `ROOT::RIO` as dependencies (rather than plain
# library names for example via ${ROOT_LIBRARIES}) ensures that properties such as required
# include directories and C++ standard are propagated to our libraries or executables.
# Note: To ensure compatibility with Cling, targets *must* be compiled using the
# same C++ standard as ROOT was compiled with.
   #add_library(RF SHARED Event.cxx)
   #target_link_libraries(RF ${ROOT_LIBRARIES})

# Create the main program using the library.
   add_executable(RF RF.cc ${sources} ${headers})
   #target_link_libraries(RF)
