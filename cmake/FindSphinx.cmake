set(_SPHINX_SCRIPT_DIR ${CMAKE_CURRENT_LIST_DIR})

message(STATUS "LIST " ${CMAKE_CURRENT_LIST_DIR})

include(FindPackageHandleStandardArgs)

# Sphinx should be near the Python interpreter
# Project needs Python
set(Python3_FIND_VIRTUALENV FIRST)
find_package(Python3
             COMPONENTS Interpreter
             REQUIRED)

message(STATUS ${Python3_EXECUTABLE})
if(Python3_Interpreter_FOUND)
    get_filename_component(_PYTHON_DIR "${Python3_EXECUTABLE}" DIRECTORY)
    set(
        _PYTHON_PATHS
        "${_PYTHON_DIR}"
    )
endif()

find_program(
    SPHINX_EXECUTABLE
    NAMES sphinx-build sphinx-build.exe
    HINTS ${_PYTHON_PATHS}
)

mark_as_advanced(SPHINX_EXECUTABLE)

find_package_handle_standard_args(Sphinx DEFAULT_MSG SPHINX_EXECUTABLE)

# If Sphinx has not been found
if(NOT Sphinx_FOUND)
    return()
endif()


# add_sphinx_document(
#   <name>
#   CONF_FILE <conf-py-filename>
#   [C_API <c-api-header-file>]
#   [SKIP_HTML] [SKIP_PDF]
#   <rst-src-file...)
#
# Function for creating Sphinx documentation targets.
function(add_sphinx_document TARGET_NAME)

    cmake_parse_arguments(
        ${TARGET_NAME}
        "SKIP_HTML;SKIP_PDF"
        "CONF_FILE"
        ""
        ${ARGN}
    )

    message(STATUS "TARGET NAME : " ${TARGET_NAME})
    message(STATUS "ARGN : " ${ARGN})
    message(STATUS "CONF : " ${${TARGET_NAME}_CONF_FILE})

    # Set up variable pfor Sphinx source, intermediate directory
    # and output directory
    get_filename_component(SRCDIR "${${TARGET_NAME}_CONF_FILE}" DIRECTORY)
    set(INTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/source")
    set(OUTDIR "${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME}/build")

    # Configure the conf.py file to include information from CMake

    string(TIMESTAMP SPHINX_TARGET_YEAR "%Y" UTC)

    add_custom_command(
        OUTPUT "${INTDIR}/conf.py"
        COMMAND "${CMAKE_COMMAND}" -E make_directory "${INTDIR}"
        COMMAND
            "${CMAKE_COMMAND}"
            "-DCONFIGURE_FILE_IN=${${TARGET_NAME}_CONF_FILE}"
            "-DCONFIGURE_FILE_OUT=${INTDIR}/conf.py"
            "-DSPHINX_TARGET_NAME=${TARGET_NAME}"
            "-DSPHINX_TARGET_VERSION=${PROJECT_VERSION}"
            "-DSPHINX_TARGET_VERSION_MAJOR=${PROJECT_VERSION_MAJOR}"
            "-DSPHINX_TARGET_VERSION_MINOR=${PROJECT_VERSION_MINOR}"
            "-DSPHINX_TARGET_YEAR=${SPHINX_TARGET_YEAR}"
            -P "${_SPHINX_SCRIPT_DIR}/BuildTimeConfigureFile.cmake"
        DEPENDS "${${TARGET_NAME}_CONF_FILE}"
    )

    set(SPHINX_DEPENDS "${INTDIR}/conf.py")

    # Loop to copy files to intermediate directory
    foreach(DOCFILE ${${TARGET_NAME}_UNPARSED_ARGUMENTS})
        message(STATUS "DOCFILE : " ${DOCFILE})
        get_filename_component(DOCFILE_INTDIR "${DOCFILE}" DIRECTORY)
        string(
            REPLACE
            "${SRCDIR}" "${INTDIR}"
            DOCFILE_INTDIR "${DOCFILE_INTDIR}"
        )

        get_filename_component(DOCFILE_DEST "${DOCFILE}" NAME)
        set(DOCFILE_DEST "${DOCFILE_INTDIR}/${DOCFILE_DEST}")

        add_custom_command(
            OUTPUT "${DOCFILE_DEST}"
            COMMAND
                "${CMAKE_COMMAND}" -E make_directory "${DOCFILE_INTDIR}"
            COMMAND
                "${CMAKE_COMMAND}" -E copy_if_different "${DOCFILE}" "${DOCFILE_DEST}"
            DEPENDS "${DOCFILE}"
        )

        list(APPEND SPHINX_DEPENDS "${DOCFILE_DEST}")

    endforeach()

    # Build html
    set(TARGET_DEPENDS)

    if(NOT ${TARGET_NAME}_SKIP_HTML)
        add_custom_command(
            OUTPUT "${OUTDIR}/html.stamp"
            # Create the _static directory required bu Sphinx in case
            # it wasn't added as one of the source files
            COMMAND "${CMAKE_COMMAND}" -E make_directory "${INTDIR}/_static"
            COMMAND "${SPHINX_EXECUTABLE}" -M html "${INTDIR}" "${OUTDIR}"
            COMMAND "${CMAKE_COMMAND}" -E touch "${OUTDIR}/html.stamp"
            DEPENDS ${SPHINX_DEPENDS}
        )

        list(APPEND TARGET_DEPENDS "${OUTDIR}/html.stamp")

    endif()

    ######
    # Add PDF generattion here
    ######

    add_custom_target(
        ${TARGET_NAME}
        DEPENDS ${TARGET_DEPENDS}
    )

    if(NOT TARGET doc)
        add_custom_target(doc)
    endif()
    add_dependencies(doc ${TARGET_NAME})

endfunction()
