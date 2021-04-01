cmake_minimum_required(VERSION 3.4)

project(tsp)

#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#find_package(CPLEX)
#if (CPLEX_FOUND)
#    add_definitions(-DIL_STD)
#else()
#    message (FATAL_ERROR "CPLEX not found!")
#endif()

include_directories(
    include
)

add_executable(${PROJECT_NAME}
    src/main.c
    src/tsp.c
    src/tsp_utilities.c
    include/tsp.h
    include/tsp_utilities.h
    include/instance.h
)

#target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIRS})

#target_link_libraries(${PROJECT_NAME} PUBLIC ${CPLEX_LIBRARIES})

#set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-lconcert -lilocplex -lcplex -lpthread -ldl")