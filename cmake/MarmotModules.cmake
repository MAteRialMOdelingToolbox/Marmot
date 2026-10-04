#/* ---------------------------------------------------------------------
# *
# * Marmot's module system.
# *
# * Every directory modules/<category>/<Name>/ containing a module.cmake is a module. Discovery is
# * automatic; the *_MODULES variables (CORE_MODULES, MATERIAL_MODULES, ...) filter it.
# *
# * A module.cmake declares its module with a single call:
# *
# *   marmot_add_module(<Name>
# *       [REQUIRES <module>...]    modules whose headers <Name> includes
# *       [SOURCES <file>...]       default: src/*.cpp
# *       [LINK <library>...])      additional libraries, e.g. a shared library the module needs
# *
# * marmot_add_module() only records the declaration. After all module.cmake files are read,
# * marmot_build_modules() resolves the dependencies and creates one OBJECT library Marmot_<Name>
# * per module:
# *
# *   - A module whose required module is missing (filtered out, or itself skipped) is skipped with
# *     a warning, transitively. If the user explicitly asked for it in a *_MODULES filter, configuring
# *     fails instead.
# *   - A dependency cycle fails configuring.
# *   - Marmot_<Name> sees the headers of its own include/ directory and those of its required modules
# *     (transitively), nothing else: an include of an undeclared module's header fails to compile.
# *
# * All module objects are assembled into the one shared library libMarmot, so materials and elements
# * registering themselves with the factories through static initializers keep doing so.
# *
# * Module.cmake files of the older style (appending to the variable `sources`, calling
# * include_directories() or marmot_module_requires()) still work: their sources are collected into
# * one object library Marmot_legacy that sees the headers of all modules. They are deprecated.
# *
# * ---------------------------------------------------------------------
# */

# ── Declaration ───────────────────────────────────────────────────────────────────────────────────

# @brief Declare the module of the module.cmake being read; see the top of this file.
# @param name The module's name; must equal its directory name, which is what *_MODULES filters and
#             other modules' REQUIRES refer to.
function(marmot_add_module name)
    cmake_parse_arguments(PARSE_ARGV 1 _arg "" "" "REQUIRES;SOURCES;LINK")
    if(_arg_UNPARSED_ARGUMENTS)
        message(FATAL_ERROR "marmot_add_module(${name}): unknown arguments: ${_arg_UNPARSED_ARGUMENTS}")
    endif()

    get_filename_component(_dir_name "${CMAKE_CURRENT_LIST_DIR}" NAME)
    if(NOT name STREQUAL _dir_name)
        message(FATAL_ERROR "marmot_add_module(${name}) in ${CMAKE_CURRENT_LIST_DIR}: "
                            "the module name must equal its directory name '${_dir_name}'.")
    endif()

    get_property(_declared GLOBAL PROPERTY MARMOT_DECLARED_MODULES)
    if(name IN_LIST _declared)
        message(FATAL_ERROR "marmot_add_module(${name}): module declared twice.")
    endif()

    if(_arg_SOURCES)
        set(_sources "")
        foreach(_source IN LISTS _arg_SOURCES)
            if(NOT IS_ABSOLUTE "${_source}")
                set(_source "${CMAKE_CURRENT_LIST_DIR}/${_source}")
            endif()
            list(APPEND _sources "${_source}")
        endforeach()
    else()
        file(GLOB _sources CONFIGURE_DEPENDS "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp")
    endif()

    set_property(GLOBAL APPEND PROPERTY MARMOT_DECLARED_MODULES "${name}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_DIR "${CMAKE_CURRENT_LIST_DIR}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_REQUIRES "${_arg_REQUIRES}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_SOURCES "${_sources}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_LINK "${_arg_LINK}")
endfunction()

# @brief Deprecated dependency check of older module.cmake files: sets <result_var> to TRUE if all
#        listed modules were discovered, warning about each missing one. Use marmot_add_module(REQUIRES).
# @param result_var Variable set in the calling scope.
# @param module_name The calling module, for the messages.
function(marmot_module_requires result_var module_name)
    set(_ok TRUE)
    foreach(_dep IN LISTS ARGN)
        if(NOT _dep IN_LIST INSTALLED_MODULES)
            message(WARNING "Module ${module_name} will NOT be built: required module '${_dep}' is not enabled.")
            set(_ok FALSE)
        endif()
    endforeach()
    set(${result_var} ${_ok} PARENT_SCOPE)
endfunction()

# ── Resolution ────────────────────────────────────────────────────────────────────────────────────

# @brief Depth-first topological sort of the modules in <modules>, following their REQUIRES within
#        <modules>; fails on a cycle.
# @param out_var Variable receiving the modules, each after the modules it requires.
function(_marmot_sort_modules out_var modules)
    set(_sorted "")
    set(_on_path "")
    # CMake has no recursion-friendly local state, so the DFS keeps an explicit stack of
    # "<module>|<state>" entries; state 0 = enter, 1 = leave.
    foreach(_root IN LISTS modules)
        if(_root IN_LIST _sorted)
            continue()
        endif()
        set(_stack "${_root}|0")
        while(_stack)
            list(POP_BACK _stack _entry)
            string(REPLACE "|" ";" _entry "${_entry}")
            list(GET _entry 0 _module)
            list(GET _entry 1 _state)
            if(_state STREQUAL "1")
                list(REMOVE_ITEM _on_path "${_module}")
                list(APPEND _sorted "${_module}")
                continue()
            endif()
            if(_module IN_LIST _sorted)
                continue()
            endif()
            if(_module IN_LIST _on_path)
                message(FATAL_ERROR "Marmot modules: dependency cycle through ${_module} (path: ${_on_path}).")
            endif()
            list(APPEND _on_path "${_module}")
            list(APPEND _stack "${_module}|1")
            get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
            foreach(_dep IN LISTS _requires)
                if(_dep IN_LIST modules AND NOT _dep IN_LIST _sorted)
                    if(_dep IN_LIST _on_path)
                        message(FATAL_ERROR "Marmot modules: dependency cycle: ${_on_path};${_dep}")
                    endif()
                    list(APPEND _stack "${_dep}|0")
                endif()
            endforeach()
        endwhile()
    endforeach()
    set(${out_var} "${_sorted}" PARENT_SCOPE)
endfunction()

# @brief Resolve the declared modules and create their targets; see the top of this file.
#        Sets in the calling scope:
#          MARMOT_BUILT_MODULES     the built modules, each after the modules it requires
#          MARMOT_MODULE_TARGETS    their targets (plus Marmot_legacy, if any older-style module exists)
#          MARMOT_INCLUDE_DIRS      the include directories of all built modules
#          MARMOT_MODULE_LINK_LIBRARIES  the LINK libraries of all built modules
# @param legacy_modules The discovered modules whose module.cmake is of the older style.
# @param legacy_sources Their sources.
# @param legacy_include_dirs Their include directories.
# @param explicit_modules The modules the user named in a *_MODULES filter.
function(marmot_build_modules legacy_modules legacy_sources legacy_include_dirs explicit_modules)
    get_property(_declared GLOBAL PROPERTY MARMOT_DECLARED_MODULES)

    # Drop modules with a missing requirement until nothing changes (handles chains).
    set(_available ${_declared} ${legacy_modules})
    set(_dropped "")
    set(_changed TRUE)
    while(_changed)
        set(_changed FALSE)
        foreach(_module IN LISTS _declared)
            if(_module IN_LIST _dropped)
                continue()
            endif()
            get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
            foreach(_dep IN LISTS _requires)
                if(NOT _dep IN_LIST _available OR _dep IN_LIST _dropped)
                    set(_reason "required module '${_dep}' is not available")
                    if(_module IN_LIST explicit_modules)
                        message(FATAL_ERROR "Module ${_module} was requested explicitly, but cannot be built: ${_reason}.")
                    endif()
                    message(WARNING "Module ${_module} will NOT be built: ${_reason}.")
                    list(APPEND _dropped "${_module}")
                    set(_changed TRUE)
                    break()
                endif()
            endforeach()
        endforeach()
    endwhile()
    set(_modules ${_declared})
    if(_dropped)
        list(REMOVE_ITEM _modules ${_dropped})
    endif()
    _marmot_sort_modules(_modules "${_modules}")

    # Settings shared by all module objects. PUBLIC through the module targets, so that the tests
    # linking them compile alike (see add_marmot_test); the shared library does not link the module
    # targets, so consumers never see MARMOT_BUILDING_LIBRARY or the coverage flags.
    add_library(MarmotModuleSettings INTERFACE)
    target_compile_definitions(MarmotModuleSettings INTERFACE MARMOT_BUILDING_LIBRARY)
    target_link_libraries(MarmotModuleSettings INTERFACE Eigen3::Eigen autodiff::autodiff Fastor::Fastor)
    if(MARMOT_ENABLE_COVERAGE)
        target_compile_options(MarmotModuleSettings INTERFACE --coverage -O0)
    endif()

    set(_targets "")
    set(_include_dirs "")
    set(_link_libraries "")
    foreach(_module IN LISTS _modules)
        get_property(_dir GLOBAL PROPERTY MARMOT_MODULE_${_module}_DIR)
        get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
        get_property(_sources GLOBAL PROPERTY MARMOT_MODULE_${_module}_SOURCES)
        get_property(_link GLOBAL PROPERTY MARMOT_MODULE_${_module}_LINK)

        set(_target Marmot_${_module})
        if(_sources)
            add_library(${_target} OBJECT ${_sources})
            set(_scope PUBLIC)
            if(MARMOT_EXPORT_API_ONLY)
                set_target_properties(${_target} PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
            endif()
        else()
            # header-only module
            add_library(${_target} INTERFACE)
            set(_scope INTERFACE)
        endif()
        target_include_directories(${_target} ${_scope} $<BUILD_INTERFACE:${_dir}/include>)
        target_link_libraries(${_target} ${_scope} MarmotModuleSettings ${_link})
        foreach(_dep IN LISTS _requires)
            if(_dep IN_LIST legacy_modules)
                # an older-style module has no target of its own; its headers are on the legacy path
                target_include_directories(${_target} ${_scope} ${legacy_include_dirs})
            else()
                target_link_libraries(${_target} ${_scope} Marmot_${_dep})
            endif()
        endforeach()

        list(APPEND _targets ${_target})
        list(APPEND _include_dirs "${_dir}/include")
        list(APPEND _link_libraries ${_link})
    endforeach()

    if(legacy_modules)
        message(DEPRECATION "These modules use the older module.cmake style; declare them with "
                            "marmot_add_module() (see cmake/MarmotModules.cmake): ${legacy_modules}")
        if(legacy_sources)
            add_library(Marmot_legacy OBJECT ${legacy_sources})
            set(_scope PUBLIC)
            if(MARMOT_EXPORT_API_ONLY)
                set_target_properties(Marmot_legacy PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
            endif()
        else()
            add_library(Marmot_legacy INTERFACE)
            set(_scope INTERFACE)
        endif()
        target_include_directories(Marmot_legacy ${_scope} ${legacy_include_dirs})
        target_link_libraries(Marmot_legacy ${_scope} ${_targets} MarmotModuleSettings)
        list(APPEND _targets Marmot_legacy)
        list(APPEND _include_dirs ${legacy_include_dirs})
    endif()

    set(MARMOT_BUILT_MODULES "${_modules}" PARENT_SCOPE)
    set(MARMOT_MODULE_TARGETS "${_targets}" PARENT_SCOPE)
    set(MARMOT_INCLUDE_DIRS "${_include_dirs}" PARENT_SCOPE)
    set(MARMOT_MODULE_LINK_LIBRARIES "${_link_libraries}" PARENT_SCOPE)
endfunction()
