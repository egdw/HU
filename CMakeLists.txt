cmake_minimum_required(VERSION 3.20)
project(HU)
find_package(OpenCV 3 REQUIRED)
#find_package( PCL 1.7 REQUIRED )
#find_package(Eigen3 3.1.0 REQUIRED)
set(CMAKE_CXX_STANDARD 14)

add_executable(HU main.cpp main2.cpp main2.h)
target_link_libraries(HU ${OpenCV_LIBS})