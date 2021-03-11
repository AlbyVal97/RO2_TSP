cmake_minimum_required(VERSION 2.8)

project(tsp)

include_directories(
    include
)

add_executable(${PROJECT_NAME}
    src/main.c
    src/chrono.c
    src/tsp.c
    src/tsp_utilities.c
    include/tsp.h
    include/tsp_utilities.h
)