cmake_minimum_required(VERSION 3.22)
project(double_pendulum LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_compile_options(-Wall -Wshadow -Wextra -O2)

set(SOURCES src/main.cpp src/RungeKutta.h src/RungeKutta.cpp src/const.h)
add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} sfml-graphics sfml-window sfml-system)

# add_custom_command(
#     TARGET ${PROJECT_NAME}
#     COMMENT "Copy assets directory"
#     PRE_BUILD COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/assets $<TARGET_FILE_DIR:${PROJECT_NAME}>/assets
# )
