cmake_minimum_required(VERSION 2.8)
project(depixelize)

include(build/CMakeLists.txt)

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${Boost_INCLUDE_DIRS}
  ${OpenCV_INCLUDE_DIRS}
)

add_subdirectory(depixelize)
