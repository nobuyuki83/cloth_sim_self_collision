project(internal_cloth_eigen)

cmake_minimum_required(VERSION 2.8)
set( CMAKE_CXX_FLAGS "-Wall -Wno-deprecated-declarations -g" )

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories(
  ${OPENGL_INCLUDE_DIR}
  ${GLUT_INCLUDE_DIR}
  ../eigen-git-mirror
)

add_executable(${PROJECT_NAME}
  main.cpp
  ../cloth_internal_physics.h  
  ../solve_internal_eigen.h
  ../utility.h
  ../vector3d.h
)

target_link_libraries(${PROJECT_NAME} 
  ${GLUT_LIBRARY} 
  ${OPENGL_LIBRARY}
)
