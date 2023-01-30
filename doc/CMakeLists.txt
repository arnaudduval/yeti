find_package(Sphinx)

if(Sphinx_FOUND)
    add_sphinx_document(
        yeti-user-manual
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
        "${CMAKE_CURRENT_LIST_DIR}/source/usage/install_wsl.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/validation/benchs.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/theory/embedded_sensitivities.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/api/preprocessing.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/api/iga_manipulation.rst"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/nodes-numbering.png"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/logo-yeti_with_text.svg"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/install_WSL/01_launch_powershell_admin.png"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/install_WSL/02_launch_Ubuntu.png"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/benchs/dload/dload_geometry.svg"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/benchs/dload/dload_triangular_increase.svg"
        "${CMAKE_CURRENT_LIST_DIR}/source/_static/benchs/dload/dload_triangular_decrease.svg"
    )


else()
    message(WARNING "No Sphinx found. Documentation target not available")
endif()
