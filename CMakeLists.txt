cmake_minimum_required(VERSION 2.8)

project(RO2_TSP)

include_directories(
    include
)

add_executable(${PROJECT_NAME}
    src/main.c
    src/chrono.c
    include/tsp.h
)