add_executable(depixelize main.cpp pixel_grid.cpp)
target_link_libraries(depixelize ${Boost_LIBRARIES} ${OpenCV_LIBRARIES})

install(TARGETS depixelize DESTINATION ${PROJECT_SOURCE_DIR}/..)
