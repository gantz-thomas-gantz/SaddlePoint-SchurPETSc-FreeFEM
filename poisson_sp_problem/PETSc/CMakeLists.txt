cmake_minimum_required(VERSION 3.10)

# Set the project name
project(poisson_sp_problem)

# Use custom PETSc installation
set(PETSC_DIR "~/petsc")
set(PETSC_ARCH "arch-debug")
set(PETSC_INCLUDE_DIRS "${PETSC_DIR}/include" "${PETSC_DIR}/${PETSC_ARCH}/include")
set(PETSC_LIBRARIES "${PETSC_DIR}/${PETSC_ARCH}/lib/libpetsc.so")

# Add the executable
add_executable(poisson_sp_problem main.cpp PETSCMatUtilities.cpp)

# Add the custom PETSC utilities modules
add_library(PETSCMatUtilities PETSCMatUtilities.cpp)
add_library(PETSCSolveUtilities PETSCSolveUtilities.cpp)

# Ensure PETSc and MPI include directories are added
include_directories(${PETSC_INCLUDE_DIRS})

# Link utilities to the main executable
target_link_libraries(poisson_sp_problem PETSCMatUtilities PETSCSolveUtilities)
# Link PETSC to the main executable
target_link_libraries(poisson_sp_problem ${PETSC_LIBRARIES})
