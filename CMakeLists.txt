cmake_minimum_required(VERSION 3.18.4)
project(sparse-multivectors C)

add_executable(${PROJECT_NAME} "")
add_subdirectory(src)

target_link_libraries(${PROJECT_NAME} PUBLIC m)
