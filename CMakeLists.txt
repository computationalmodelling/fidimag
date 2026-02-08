cmake_minimum_required(VERSION 3.18...3.28)
project(fidimag LANGUAGES C)

# Policy settings for modern CMake
if(POLICY CMP0167)
    cmake_policy(SET CMP0167 NEW)  # FindBoost compatibility
endif()

# =============================================================================
# Build Configuration
# =============================================================================

# Default to Release build for optimization
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "RelWithDebInfo")
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

# C/C++ standards
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

# Position independent code (required for Python extensions)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# =============================================================================
# Add local library paths (for SUNDIALS/FFTW in ./local/)
# =============================================================================

list(APPEND CMAKE_PREFIX_PATH
    "${PROJECT_SOURCE_DIR}/local"
    "${PROJECT_SOURCE_DIR}/local/lib/cmake"
    "${PROJECT_SOURCE_DIR}/local/lib64/cmake"
)

# Also add to module path for custom Find modules
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")

message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

# =============================================================================
# Find Dependencies
# =============================================================================

# Python, NumPy, and Cython
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
message(STATUS "Python: ${Python_VERSION} (${Python_EXECUTABLE})")
message(STATUS "NumPy include: ${Python_NumPy_INCLUDE_DIRS}")

# Cython
find_program(CYTHON_EXECUTABLE NAMES cython cython3 REQUIRED)
message(STATUS "Cython: ${CYTHON_EXECUTABLE}")

# OpenMP (required for parallel computation)
find_package(OpenMP REQUIRED COMPONENTS C)
message(STATUS "OpenMP found: C ${OpenMP_C_FLAGS}")

# BLAS and LAPACK (linear algebra)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
message(STATUS "BLAS: ${BLAS_LIBRARIES}")
message(STATUS "LAPACK: ${LAPACK_LIBRARIES}")

# FFTW3 (Fast Fourier Transform)
# Look for both regular and OpenMP versions
find_library(FFTW3_LIBRARY NAMES fftw3 REQUIRED
    HINTS ${PROJECT_SOURCE_DIR}/local/lib ${PROJECT_SOURCE_DIR}/local/lib64
    DOC "FFTW3 library"
)
find_library(FFTW3_OMP_LIBRARY NAMES fftw3_omp REQUIRED
    HINTS ${PROJECT_SOURCE_DIR}/local/lib ${PROJECT_SOURCE_DIR}/local/lib64
    DOC "FFTW3 OpenMP library"
)
find_path(FFTW3_INCLUDE_DIR fftw3.h
    HINTS ${PROJECT_SOURCE_DIR}/local/include $ENV{FFTW_INC}
    DOC "FFTW3 include directory"
)
message(STATUS "FFTW3: ${FFTW3_LIBRARY}")
message(STATUS "FFTW3 OpenMP: ${FFTW3_OMP_LIBRARY}")
message(STATUS "FFTW3 include: ${FFTW3_INCLUDE_DIR}")

# SUNDIALS (ODE/DAE solver) - Custom find module
find_package(SUNDIALS 7.6 REQUIRED COMPONENTS cvodes nvecserial nvecopenmp)
message(STATUS "SUNDIALS ${SUNDIALS_VERSION} found")
message(STATUS "  Include: ${SUNDIALS_INCLUDE_DIRS}")
message(STATUS "  Libraries: ${SUNDIALS_LIBRARIES}")

# =============================================================================
# Common Compiler Flags
# =============================================================================

# C compiler flags
set(C_COMPILE_FLAGS
    -O3
    -Wall
    -Wno-cpp
    -Wno-unused-function
)

# Common libraries for all extensions
set(COMMON_LIBRARIES
    m  # Math library
    ${FFTW3_OMP_LIBRARY}
    ${FFTW3_LIBRARY}
    ${SUNDIALS_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
    OpenMP::OpenMP_C
)

# Common include directories
set(COMMON_INCLUDES
    ${Python_INCLUDE_DIRS}
    ${Python_NumPy_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/local/include
    ${FFTW3_INCLUDE_DIR}
    ${SUNDIALS_INCLUDE_DIRS}
)

# Add environment variable paths if set
if(DEFINED ENV{SUNDIALS_INC})
    list(APPEND COMMON_INCLUDES $ENV{SUNDIALS_INC})
    message(STATUS "Added SUNDIALS_INC: $ENV{SUNDIALS_INC}")
endif()

if(DEFINED ENV{FFTW_INC})
    list(APPEND COMMON_INCLUDES $ENV{FFTW_INC})
    message(STATUS "Added FFTW_INC: $ENV{FFTW_INC}")
endif()

# =============================================================================
# Helper Function: Create Cython Extension
# =============================================================================

function(add_fidimag_extension MODULE_NAME SOURCE_FILE)
    # Get the module name without "fidimag.extensions." prefix for target name
    string(REPLACE "." "_" TARGET_NAME ${MODULE_NAME})

    # Determine output paths
    get_filename_component(SOURCE_DIR ${SOURCE_FILE} DIRECTORY)
    get_filename_component(SOURCE_NAME ${SOURCE_FILE} NAME_WE)

    # Cython output file (always C)
    set(CYTHON_OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${SOURCE_DIR}/${SOURCE_NAME}.c")

    # Find additional C/C++ source files in the same directory
    file(GLOB EXTRA_SOURCES
        "${PROJECT_SOURCE_DIR}/fidimag/${SOURCE_DIR}/*.c"
        "${PROJECT_SOURCE_DIR}/fidimag/${SOURCE_DIR}/*.cpp"
    )

    # Remove the cythonized output if it's in the source dir
    list(FILTER EXTRA_SOURCES EXCLUDE REGEX ".*${SOURCE_NAME}\\.c$")
    list(FILTER EXTRA_SOURCES EXCLUDE REGEX ".*${SOURCE_NAME}\\.cpp$")

    message(STATUS "Extension ${MODULE_NAME}:")
    message(STATUS "  Source: ${SOURCE_FILE}")
    message(STATUS "  Cython output: ${CYTHON_OUTPUT}")
    message(STATUS "  Extra sources: ${EXTRA_SOURCES}")
    message(STATUS "  Language: C")

    # Ensure output directory exists
    get_filename_component(CYTHON_OUTPUT_DIR ${CYTHON_OUTPUT} DIRECTORY)
    file(MAKE_DIRECTORY ${CYTHON_OUTPUT_DIR})

    # Cythonize the .pyx file
    add_custom_command(
        OUTPUT ${CYTHON_OUTPUT}
        COMMAND ${CYTHON_EXECUTABLE}
            -3  # Python 3
            --directive linetrace=True
            --directive language_level=3
            -o ${CYTHON_OUTPUT}
            ${PROJECT_SOURCE_DIR}/fidimag/${SOURCE_FILE}
        DEPENDS ${PROJECT_SOURCE_DIR}/fidimag/${SOURCE_FILE}
        COMMENT "Cythonizing ${SOURCE_FILE}"
    )

    set_source_files_properties(${CYTHON_OUTPUT} PROPERTIES LANGUAGE C)
    set_source_files_properties(${EXTRA_SOURCES} PROPERTIES LANGUAGE C)

    # Create the Python extension module
    Python_add_library(${TARGET_NAME} MODULE
        ${CYTHON_OUTPUT}
        ${EXTRA_SOURCES}
        WITH_SOABI
    )

    # Set target properties
    set_target_properties(${TARGET_NAME} PROPERTIES
        OUTPUT_NAME ${SOURCE_NAME}
        PREFIX ""  # No lib prefix
        LINKER_LANGUAGE C
    )

    # Include directories
    target_include_directories(${TARGET_NAME} PRIVATE
        ${COMMON_INCLUDES}
        ${PROJECT_SOURCE_DIR}/fidimag/${SOURCE_DIR}
    )

    # Compiler flags
    target_compile_options(${TARGET_NAME} PRIVATE ${C_COMPILE_FLAGS})

    # Link libraries
    target_link_libraries(${TARGET_NAME} PRIVATE ${COMMON_LIBRARIES})

    # RPATH for finding local libraries at runtime
    set_target_properties(${TARGET_NAME} PROPERTIES
        INSTALL_RPATH "${PROJECT_SOURCE_DIR}/local/lib;${PROJECT_SOURCE_DIR}/local/lib64"
        BUILD_RPATH "${PROJECT_SOURCE_DIR}/local/lib;${PROJECT_SOURCE_DIR}/local/lib64"
    )

    # Install the extension to the proper location in the fidimag package
    # Extract the Python module path from MODULE_NAME
    string(REPLACE "." "/" MODULE_PATH ${MODULE_NAME})
    get_filename_component(INSTALL_DIR ${MODULE_PATH} DIRECTORY)

    install(TARGETS ${TARGET_NAME}
        LIBRARY DESTINATION ${INSTALL_DIR}
        COMPONENT python
    )
endfunction()

# =============================================================================
# Define All Cython Extensions
# =============================================================================

# Atomistic extensions
add_fidimag_extension("fidimag.extensions.clib" "atomistic/lib/clib.pyx")

# Micromagnetic extensions
add_fidimag_extension("fidimag.extensions.micro_clib" "micro/lib/micro_clib.pyx")
add_fidimag_extension("fidimag.extensions.baryakhtar_clib" "micro/lib/baryakhtar/baryakhtar_clib.pyx")

# Common extensions
add_fidimag_extension("fidimag.extensions.common_clib" "common/lib/common_clib.pyx")
add_fidimag_extension("fidimag.extensions.cvode" "common/sundials/cvode.pyx")
add_fidimag_extension("fidimag.extensions.dipolar" "common/dipolar/dipolar.pyx")
add_fidimag_extension("fidimag.extensions.nebm_clib" "common/neb_method/nebm_clib.pyx")

# =============================================================================
# Print Configuration Summary
# =============================================================================

message(STATUS "")
message(STATUS "========================================")
message(STATUS "fidimag Build Configuration Summary")
message(STATUS "========================================")
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C compiler: ${CMAKE_C_COMPILER}")
message(STATUS "Python: ${Python_VERSION}")
message(STATUS "NumPy: ${Python_NumPy_VERSION}")
message(STATUS "SUNDIALS: ${SUNDIALS_VERSION}")
message(STATUS "OpenMP: Enabled")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "========================================")
message(STATUS "")
