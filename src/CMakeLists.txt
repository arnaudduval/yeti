# Grab Lapack and openMP
find_package(
    LAPACK
    REQUIRED
)

find_package(OpenMP REQUIRED)

# define a macro to build f2py module
macro(CreateF2PYModule)
    set(options)
    set(oneValueArgs NAME DESTINATION)
    set(multiValueArgs FILES LINK_LIBRARIES SKIP)

    CMAKE_PARSE_ARGUMENTS(MODULE "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})

    set(generated_module_file ${CMAKE_CURRENT_BINARY_DIR}/${MODULE_NAME}.${Python3_SOABI}.so)

    add_custom_target(${MODULE_NAME} ALL
        DEPENDS ${generated_module_file}
    )

    # manage dependencies and crete link command
    set(link_string )
    foreach(dep ${MODULE_LINK_LIBRARIES})
        set(link_string ${link_string} -l${dep})
    endforeach()

    # manage skipped functions
    set(skip_string)
    foreach(fun ${MODULE_SKIP})
        set(skip_string ${skip_string} ${fun})
    endforeach()


    if(LEGACY)
        add_custom_command(
            OUTPUT ${generated_module_file}
            COMMAND "${Python3_EXECUTABLE}" -m numpy.f2py
            -m ${MODULE_NAME} #--quiet
            -c --opt=-O2
            -I${CMAKE_Fortran_MODULE_DIRECTORY}
            -L${LIBRARY_OUTPUT_PATH}
            ${link_string}
            -llapack -lblas
            ${FCOMPILER}
            --f90flags='-fopenmp' -lgomp
            ${MODULE_FILES}
            skip: ${skip_string}
            DEPENDS ${MODULE_FILES}
            DEPENDS ${MODULE_LINK_LIBRARIES}
            COMMENT "[F2PY] Building Fortran to Python interface module ${MODULE_NAME}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            COMMAND  ${CMAKE_COMMAND} -E copy ${generated_module_file} ${MODULE_DESTINATION}/${MODULE_NAME}.${Python3_SOABI}.so
        )
    else()
        # scikit-build : no copy of *.so file
        add_custom_command(
            OUTPUT ${generated_module_file}
            COMMAND "${Python3_EXECUTABLE}" -m numpy.f2py
            -m ${MODULE_NAME} #--quiet
            -c --opt=-O2
            -I${CMAKE_Fortran_MODULE_DIRECTORY}
            -L${LIBRARY_OUTPUT_PATH}
            ${link_string}
            -llapack -lblas
            ${FCOMPILER}
            --f90flags='-fopenmp' -lgomp
            ${MODULE_FILES}
            skip: ${skip_string}
            DEPENDS ${MODULE_FILES}
            DEPENDS ${MODULE_LINK_LIBRARIES}
            COMMENT "[F2PY] Building Fortran to Python interface module ${MODULE_NAME}"
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            # COMMAND  ${CMAKE_COMMAND} -E copy ${generated_module_file} ${MODULE_DESTINATION}/${MODULE_NAME}.${Python3_SOABI}.so
        )
        install(FILES ${generated_module_file} DESTINATION ${MODULE_DESTINATION})

    endif()

endmacro()


# Fortran object library for use with scikit-build
add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python3::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)


subdirs(yeti_iga)