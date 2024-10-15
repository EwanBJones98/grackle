
# overwrite the global variable that controls where the binaries compiled in
# this directory are placed
# -> because this file is only read by CMake via add_subdirectory and we are
#    are not making a CACHE variable, this won't affect anything outside of
#    this directory
# -> for a build directory, <build>, this places the binaries at 
#    <build>/examples
# -> assuming that the build-directory is at the root level of the grackle
#    repository, then executing the examples from within the examples directory
#    should ensure that the paths to the datafiles are correct
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../test)

add_executable(contiguous_looping contiguous_looping.c)
target_link_libraries(contiguous_looping Grackle::Grackle)