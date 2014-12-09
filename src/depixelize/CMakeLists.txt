add_executable(depixelize main.cpp pixel_grid.cpp spline_optimizer.cpp)
target_link_libraries(depixelize render ${Boost_LIBRARIES} ${OpenCV_LIBRARIES})

install(TARGETS depixelize DESTINATION ${PROJECT_SOURCE_DIR}/..)
