find_package(Sphinx)

if(Sphinx_FOUND)
    add_sphinx_document(
        my-user-manual
        CONF_FILE "${CMAKE_CURRENT_LIST_DIR}/source/conf.py"
        "${CMAKE_CURRENT_LIST_DIR}/source/YETI.bib"
        "${CMAKE_CURRENT_LIST_DIR}/source/index.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/installation.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/using_yeti.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/coupling_modelling.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/load.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/elements_library.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/connectivity_convention.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/build_cmake.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/nodes-numbering.png"

    )


else()
    message(WARNING "No Sphinx found. Documentation target not available")
endif()
