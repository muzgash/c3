cmake_minimum_required(VERSION 3.0)

include_directories(../inc)

add_library(c3 SHARED
../inc/c3.h
../src/c3.c)

add_executable(fft fft/fft_example.c)
target_link_libraries(fft c3 fftw3f fftw3 m)

