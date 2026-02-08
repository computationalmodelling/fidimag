# FindSUNDIALS.cmake
# Finds SUNDIALS library (v7.x)
#
# This module defines:
#  SUNDIALS_FOUND         - True if SUNDIALS is found
#  SUNDIALS_VERSION       - Version of SUNDIALS found
#  SUNDIALS_INCLUDE_DIRS  - Include directories for SUNDIALS
#  SUNDIALS_LIBRARIES     - Libraries to link against
#  SUNDIALS_<C>_LIBRARY   - Path to component C (e.g., SUNDIALS_CVODES_LIBRARY)
#
# Components supported:
#  cvodes, nvecserial, nvecopenmp, core
#
# Users can set SUNDIALS_ROOT to specify installation directory

# Use SUNDIALS_ROOT environment variable or CMake variable
if(NOT SUNDIALS_ROOT AND DEFINED ENV{SUNDIALS_ROOT})
    set(SUNDIALS_ROOT $ENV{SUNDIALS_ROOT})
endif()

# Find include directory
find_path(SUNDIALS_INCLUDE_DIR
    NAMES sundials/sundials_config.h
    HINTS
        ${SUNDIALS_ROOT}
        ${CMAKE_PREFIX_PATH}
    PATH_SUFFIXES include
    DOC "SUNDIALS include directory"
)

# Find sundials_config.h and extract version
if(SUNDIALS_INCLUDE_DIR)
    file(READ "${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h" SUNDIALS_CONFIG_H)

    # Try to extract version (SUNDIALS 7.x format)
    string(REGEX MATCH "#define SUNDIALS_VERSION \"([0-9]+)\\.([0-9]+)\\.([0-9]+)\"" _ ${SUNDIALS_CONFIG_H})
    if(CMAKE_MATCH_1)
        set(SUNDIALS_VERSION_MAJOR ${CMAKE_MATCH_1})
        set(SUNDIALS_VERSION_MINOR ${CMAKE_MATCH_2})
        set(SUNDIALS_VERSION_PATCH ${CMAKE_MATCH_3})
        set(SUNDIALS_VERSION "${SUNDIALS_VERSION_MAJOR}.${SUNDIALS_VERSION_MINOR}.${SUNDIALS_VERSION_PATCH}")
    endif()
endif()

# Components to find
set(SUNDIALS_COMPONENTS_AVAILABLE core cvodes nvecserial nvecopenmp)

# Determine which components to search for
if(SUNDIALS_FIND_COMPONENTS)
    set(SUNDIALS_COMPONENTS_TO_FIND ${SUNDIALS_FIND_COMPONENTS})
else()
    set(SUNDIALS_COMPONENTS_TO_FIND ${SUNDIALS_COMPONENTS_AVAILABLE})
endif()

# Find each component library
foreach(COMPONENT ${SUNDIALS_COMPONENTS_TO_FIND})
    string(TOUPPER ${COMPONENT} COMPONENT_UPPER)

    find_library(SUNDIALS_${COMPONENT_UPPER}_LIBRARY
        NAMES sundials_${COMPONENT}
        HINTS
            ${SUNDIALS_ROOT}
            ${CMAKE_PREFIX_PATH}
        PATH_SUFFIXES lib lib64
        DOC "SUNDIALS ${COMPONENT} library"
    )

    # Mark as advanced
    mark_as_advanced(SUNDIALS_${COMPONENT_UPPER}_LIBRARY)

    # Add to list if found
    if(SUNDIALS_${COMPONENT_UPPER}_LIBRARY)
        list(APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${COMPONENT_UPPER}_LIBRARY})
        set(SUNDIALS_${COMPONENT}_FOUND TRUE)
    else()
        set(SUNDIALS_${COMPONENT}_FOUND FALSE)
        if(SUNDIALS_FIND_REQUIRED_${COMPONENT})
            message(FATAL_ERROR "Required SUNDIALS component ${COMPONENT} not found")
        endif()
    endif()
endforeach()

# Remove duplicates
if(SUNDIALS_LIBRARIES)
    list(REMOVE_DUPLICATES SUNDIALS_LIBRARIES)
endif()

# Set SUNDIALS_INCLUDE_DIRS (plural)
set(SUNDIALS_INCLUDE_DIRS ${SUNDIALS_INCLUDE_DIR})

# Handle standard arguments
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS
    REQUIRED_VARS SUNDIALS_INCLUDE_DIR SUNDIALS_LIBRARIES
    VERSION_VAR SUNDIALS_VERSION
    HANDLE_COMPONENTS
)

# Mark cache variables as advanced
mark_as_advanced(
    SUNDIALS_INCLUDE_DIR
    SUNDIALS_VERSION
)

# Create imported target if found
if(SUNDIALS_FOUND AND NOT TARGET SUNDIALS::SUNDIALS)
    add_library(SUNDIALS::SUNDIALS INTERFACE IMPORTED)
    set_target_properties(SUNDIALS::SUNDIALS PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${SUNDIALS_LIBRARIES}"
    )

    # Create component-specific targets
    foreach(COMPONENT ${SUNDIALS_COMPONENTS_TO_FIND})
        string(TOUPPER ${COMPONENT} COMPONENT_UPPER)
        if(SUNDIALS_${COMPONENT_UPPER}_LIBRARY AND NOT TARGET SUNDIALS::${COMPONENT})
            add_library(SUNDIALS::${COMPONENT} UNKNOWN IMPORTED)
            set_target_properties(SUNDIALS::${COMPONENT} PROPERTIES
                IMPORTED_LOCATION "${SUNDIALS_${COMPONENT_UPPER}_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_INCLUDE_DIRS}"
            )
        endif()
    endforeach()
endif()

# Debug output
if(SUNDIALS_FOUND)
    message(STATUS "Found SUNDIALS ${SUNDIALS_VERSION}")
    message(STATUS "  Include dir: ${SUNDIALS_INCLUDE_DIR}")
    message(STATUS "  Libraries: ${SUNDIALS_LIBRARIES}")
endif()
