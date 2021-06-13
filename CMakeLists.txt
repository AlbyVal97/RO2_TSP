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
	src/convex_hull.c
    include/tsp.h
    include/tsp_utilities.h
    include/instance.h
	include/convex_hull.h
)

add_library(${PROJECT_NAME}_concorde_lib
    src/mincut.c
    src/connect.c
    src/util.c
    src/shrink.c
    src/cut_st.c
    src/allocrus.c
    src/sortrus.c
    src/urandom.c
    include/util.h
    include/machdefs.h
    include/macrorus.h
    include/cut.h
    include/config.h
)

target_link_libraries(${PROJECT_NAME}
    ${PROJECT_NAME}_concorde_lib
)

#target_include_directories(${PROJECT_NAME} PUBLIC ${CPLEX_INCLUDE_DIRS})

#target_link_libraries(${PROJECT_NAME} PUBLIC ${CPLEX_LIBRARIES})

#set_target_properties(${PROJECT_NAME} PROPERTIES LINK_FLAGS "-lconcert -lilocplex -lcplex -lpthread -ldl")

# Create a library called "Hello" which includes the source file "hello.cxx".
# The extension is already found. Any number of sources could be listed here.
#add_library (Hello hello.cxx)

# Make sure the compiler can find include files for our Hello library
# when other libraries or executables link to Hello
#target_include_directories (Hello PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add executable called "helloDemo" that is built from the source files
# "demo.cxx" and "demo_b.cxx". The extensions are automatically found.
#add_executable (helloDemo demo.cxx demo_b.cxx)

# Link the executable to the Hello library. Since the Hello library has
# public include directories we will use those link directories when building
# helloDemo
#target_link_libraries (helloDemo LINK_PUBLIC Hello)