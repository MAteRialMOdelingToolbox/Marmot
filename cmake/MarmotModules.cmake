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
# *       [REQUIRES <module>...]    modules whose headers <Name> includes (see below)
# *       [SOURCES <file>...]       default: src/*.cpp
# *       [LINK <library>...])      additional libraries, e.g. a shared library the module needs; an imported
# *                                 target given here must also be found by consumers of the installed Marmot
# *
# * REQUIRES lists every module whose headers the module includes, except modules already required (directly or
# * transitively) by a listed module. A module requiring a module that does not exist fails configuring, unless a
# * *_MODULES filter is set: then it is skipped like a module whose requirement is filtered out (e.g. a module whose
# * dependency's repository is not checked out, on a machine building a filtered selection).
# *
# * marmot_add_module() only records the declaration. After all module.cmake files are read,
# * marmot_build_modules() resolves the dependencies and creates one OBJECT library Marmot_<Name> per module:
# *
# *   - A module whose required module is not built is skipped with a warning, transitively. Configuring fails
# *     instead if the module was requested explicitly in a *_MODULES filter, or if no filter is set at all.
# *   - A dependency cycle fails configuring.
# *   - Marmot_<Name> sees the headers of its own include/ directory and those of its required modules
# *     (transitively), nothing else: an include of an undeclared module's header fails to compile.
# *
# * All module objects are assembled into the one shared library libMarmot, so materials and elements
# * registering themselves with the factories through static initializers keep doing so.
# *
# * Module.cmake files of the older style (appending to `sources` and `INSTALLED_MODULE_INCLUDE_DIRS`, calling
# * include_directories() or marmot_module_requires()) still work, with a deprecation warning: such a module is
# * built if its module.cmake added sources or include directories, its marmot_module_requires() dependencies take
# * part in the resolution above, and the sources of all of them are compiled into one object library
# * Marmot_legacy, which sees the headers of all modules. A module of the new style requiring one of the older
# * style sees the headers of all modules, too, so its REQUIRES are not checked by compiling. Include directories
# * outside modules/ added by an older module.cmake are used to compile Marmot_legacy only, never installed.
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

    get_property(_declared GLOBAL PROPERTY MARMOT_MODULE_${name}_DECLARED)
    if(_declared)
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

    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_DECLARED TRUE)
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_REQUIRES "${_arg_REQUIRES}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_SOURCES "${_sources}")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/include")
    set_property(GLOBAL PROPERTY MARMOT_MODULE_${name}_LINK "${_arg_LINK}")
endfunction()

# @brief Deprecated dependency check of older module.cmake files: sets <result_var> to TRUE if all listed modules
#        were discovered, warning about each missing one. The dependencies are recorded for the resolution in
#        marmot_build_modules(). Use marmot_add_module(REQUIRES) instead.
# @param result_var Variable set in the calling scope.
# @param module_name The calling module, for the messages.
function(marmot_module_requires result_var module_name)
    get_filename_component(_module "${CMAKE_CURRENT_LIST_DIR}" NAME)
    set_property(GLOBAL APPEND PROPERTY MARMOT_MODULE_${_module}_REQUIRES ${ARGN})
    set(_ok TRUE)
    foreach(_dep IN LISTS ARGN)
        if(NOT _dep IN_LIST INSTALLED_MODULES)
            message(WARNING "Module ${module_name} will NOT be built: required module '${_dep}' is not enabled.")
            set(_ok FALSE)
        endif()
    endforeach()
    set(${result_var} ${_ok} PARENT_SCOPE)
endfunction()

# ── Reading the module.cmake files ────────────────────────────────────────────────────────────────

# @brief Read the module.cmake of the module in <module_dir>. In a function scope, so that what a module.cmake
#        of the older style appends to `sources`, `INSTALLED_MODULE_INCLUDE_DIRS` and `SHARED_LIBRARIES` is
#        captured per module; its include_directories() are captured from the directory property.
# @param module_dir The module's directory.
function(_marmot_read_module module_dir)
    get_filename_component(_module "${module_dir}" NAME)
    get_directory_property(_dir_includes_before INCLUDE_DIRECTORIES)
    set(sources "")
    set(INSTALLED_MODULE_INCLUDE_DIRS "")
    set(SHARED_LIBRARIES "")

    include("${module_dir}/module.cmake")

    get_directory_property(_dir_includes_after INCLUDE_DIRECTORIES)
    set(_added_dirs ${INSTALLED_MODULE_INCLUDE_DIRS})
    foreach(_dir IN LISTS _dir_includes_after)
        if(NOT _dir IN_LIST _dir_includes_before)
            list(APPEND _added_dirs "${_dir}")
        endif()
    endforeach()
    # Directories within modules/ hold module headers (installed, and part of libMarmot's interface); any other
    # directory is used only to compile Marmot_legacy.
    set(_include_dirs "")
    set(_extra_include_dirs "")
    foreach(_dir IN LISTS _added_dirs)
        string(FIND "${_dir}" "${MODULES_DIR}/" _pos)
        if(_pos EQUAL 0)
            list(APPEND _include_dirs "${_dir}")
        else()
            list(APPEND _extra_include_dirs "${_dir}")
        endif()
    endforeach()
    # include_directories() of an older module.cmake would apply to every target of this directory; the module's
    # directories are applied to Marmot_legacy only.
    set_directory_properties(PROPERTIES INCLUDE_DIRECTORIES "${_dir_includes_before}")

    get_property(_declared GLOBAL PROPERTY MARMOT_MODULE_${_module}_DECLARED)
    if(_declared)
        if(sources OR _added_dirs OR SHARED_LIBRARIES)
            message(FATAL_ERROR "${module_dir}/module.cmake calls marmot_add_module() and also uses the older "
                                "style (sources, include directories or SHARED_LIBRARIES); use marmot_add_module() only.")
        endif()
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY FALSE)
        return()
    endif()

    set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY TRUE)
    if(sources OR _added_dirs)
        list(REMOVE_DUPLICATES _include_dirs)
        list(REMOVE_DUPLICATES _extra_include_dirs)
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_DECLARED TRUE)
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_SOURCES "${sources}")
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_INCLUDE_DIRS "${_include_dirs}")
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_EXTRA_INCLUDE_DIRS "${_extra_include_dirs}")
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_LINK "${SHARED_LIBRARIES}")
    else()
        # it gated itself off (and said so)
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_DECLARED FALSE)
    endif()
endfunction()

# @brief Read the module.cmake files of all modules in <module_dirs>.
# @param module_dirs The directories of the discovered modules.
function(marmot_read_modules module_dirs)
    foreach(_module_dir IN LISTS module_dirs)
        get_filename_component(_module "${_module_dir}" NAME)
        get_property(_known_dir GLOBAL PROPERTY MARMOT_MODULE_${_module}_DIR)
        if(_known_dir AND NOT _known_dir STREQUAL _module_dir)
            message(FATAL_ERROR "Two modules are named ${_module}: ${_known_dir} and ${_module_dir}.")
        endif()
        set_property(GLOBAL PROPERTY MARMOT_MODULE_${_module}_DIR "${_module_dir}")
        _marmot_read_module("${_module_dir}")
    endforeach()
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

# @brief The target providing <module>: Marmot_<module>, or Marmot_legacy for a module of the older style.
# @param out_var Variable receiving the target name.
# @param module A built module.
function(marmot_module_target out_var module)
    get_property(_legacy GLOBAL PROPERTY MARMOT_MODULE_${module}_LEGACY)
    if(_legacy)
        set(${out_var} Marmot_legacy PARENT_SCOPE)
    else()
        set(${out_var} Marmot_${module} PARENT_SCOPE)
    endif()
endfunction()

# @brief Resolve the read modules and create their targets; see the top of this file.
#        Sets in the calling scope:
#          MARMOT_BUILT_MODULES          the built modules, each after the modules it requires
#          MARMOT_LEGACY_MODULES         those of them of the older style
#          MARMOT_MODULE_TARGETS         the targets of the built modules (Marmot_legacy once)
#          MARMOT_INCLUDE_DIRS           the include directories of the built modules
#          MARMOT_MODULE_LINK_LIBRARIES  the LINK libraries of the built modules
# @param discovered_modules The modules whose module.cmake was read.
# @param existing_modules All modules present in modules/, including those filtered out.
# @param explicit_modules The modules the user named in a *_MODULES filter.
# @param filtered TRUE if any *_MODULES filter is set.
function(marmot_build_modules discovered_modules existing_modules explicit_modules filtered)
    set(_candidates "")
    foreach(_module IN LISTS discovered_modules)
        get_property(_declared GLOBAL PROPERTY MARMOT_MODULE_${_module}_DECLARED)
        get_property(_legacy GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY)
        get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
        # Without a filter, a required module that does not exist is an error (typically a misspelled REQUIRES);
        # with one, it is treated like a filtered-out module below.
        if(NOT _legacy AND NOT filtered)
            foreach(_dep IN LISTS _requires)
                if(NOT _dep IN_LIST existing_modules)
                    message(FATAL_ERROR "Module ${_module} requires '${_dep}', which does not exist in modules/. "
                                        "Fix the name if it is misspelled; if its repository is not checked out, "
                                        "check it out or exclude ${_module} with a *_MODULES filter.")
                endif()
            endforeach()
        endif()
        if(_declared)
            list(APPEND _candidates "${_module}")
        elseif(_module IN_LIST explicit_modules)
            message(FATAL_ERROR "Module ${_module} was requested explicitly, but cannot be built (see above).")
        endif()
    endforeach()

    # Drop modules with a requirement that is not built until nothing changes (handles chains).
    set(_dropped "")
    set(_changed TRUE)
    while(_changed)
        set(_changed FALSE)
        foreach(_module IN LISTS _candidates)
            if(_module IN_LIST _dropped)
                continue()
            endif()
            get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
            get_property(_legacy GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY)
            foreach(_dep IN LISTS _requires)
                if(NOT _dep IN_LIST _candidates OR _dep IN_LIST _dropped)
                    if(_dep IN_LIST existing_modules)
                        set(_reason "required module '${_dep}' is not built")
                    else()
                        set(_reason "required module '${_dep}' does not exist in modules/")
                    endif()
                    if(_module IN_LIST explicit_modules)
                        message(FATAL_ERROR "Module ${_module} was requested explicitly, but cannot be built: ${_reason}.")
                    endif()
                    if(NOT filtered AND NOT _legacy)
                        message(FATAL_ERROR "Module ${_module} cannot be built: ${_reason}. Without a *_MODULES "
                                            "filter, every module must be buildable; exclude ${_module} with a filter.")
                    endif()
                    message(WARNING "Module ${_module} will NOT be built: ${_reason}.")
                    list(APPEND _dropped "${_module}")
                    set(_changed TRUE)
                    break()
                endif()
            endforeach()
        endforeach()
    endwhile()
    set(_modules ${_candidates})
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

    set(_legacy_modules "")
    set(_include_dirs "")
    set(_link_libraries "")
    foreach(_module IN LISTS _modules)
        get_property(_module_include_dirs GLOBAL PROPERTY MARMOT_MODULE_${_module}_INCLUDE_DIRS)
        get_property(_link GLOBAL PROPERTY MARMOT_MODULE_${_module}_LINK)
        get_property(_legacy GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY)
        list(APPEND _include_dirs ${_module_include_dirs})
        list(APPEND _link_libraries ${_link})
        if(_legacy)
            list(APPEND _legacy_modules "${_module}")
        endif()
    endforeach()

    set(_targets "")
    foreach(_module IN LISTS _modules)
        get_property(_legacy GLOBAL PROPERTY MARMOT_MODULE_${_module}_LEGACY)
        if(_legacy)
            continue()
        endif()
        get_property(_requires GLOBAL PROPERTY MARMOT_MODULE_${_module}_REQUIRES)
        get_property(_sources GLOBAL PROPERTY MARMOT_MODULE_${_module}_SOURCES)
        get_property(_module_include_dirs GLOBAL PROPERTY MARMOT_MODULE_${_module}_INCLUDE_DIRS)
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
        target_include_directories(${_target} ${_scope} ${_module_include_dirs})
        target_link_libraries(${_target} ${_scope} MarmotModuleSettings ${_link})
        foreach(_dep IN LISTS _requires)
            if(_dep IN_LIST _legacy_modules)
                # a module of the older style declares no dependencies of its own; see the top of this file
                target_include_directories(${_target} ${_scope} ${_include_dirs})
            else()
                target_link_libraries(${_target} ${_scope} Marmot_${_dep})
            endif()
        endforeach()
        list(APPEND _targets ${_target})
    endforeach()

    if(_legacy_modules)
        message(DEPRECATION "These modules use the older module.cmake style; declare them with "
                            "marmot_add_module() (see cmake/MarmotModules.cmake): ${_legacy_modules}")
        set(_legacy_sources "")
        set(_legacy_include_dirs "")
        set(_legacy_link "")
        foreach(_module IN LISTS _legacy_modules)
            get_property(_sources GLOBAL PROPERTY MARMOT_MODULE_${_module}_SOURCES)
            get_property(_module_include_dirs GLOBAL PROPERTY MARMOT_MODULE_${_module}_INCLUDE_DIRS)
            get_property(_extra_include_dirs GLOBAL PROPERTY MARMOT_MODULE_${_module}_EXTRA_INCLUDE_DIRS)
            get_property(_link GLOBAL PROPERTY MARMOT_MODULE_${_module}_LINK)
            list(APPEND _legacy_sources ${_sources})
            list(APPEND _legacy_include_dirs ${_module_include_dirs} ${_extra_include_dirs})
            list(APPEND _legacy_link ${_link})
        endforeach()
        if(_legacy_sources)
            add_library(Marmot_legacy OBJECT ${_legacy_sources})
            set(_scope PUBLIC)
            if(MARMOT_EXPORT_API_ONLY)
                set_target_properties(Marmot_legacy PROPERTIES CXX_VISIBILITY_PRESET hidden VISIBILITY_INLINES_HIDDEN ON)
            endif()
        else()
            add_library(Marmot_legacy INTERFACE)
            set(_scope INTERFACE)
        endif()
        target_include_directories(Marmot_legacy ${_scope} ${_legacy_include_dirs})
        target_link_libraries(Marmot_legacy ${_scope} ${_targets} MarmotModuleSettings ${_legacy_link})
        list(APPEND _targets Marmot_legacy)
    endif()

    set(MARMOT_BUILT_MODULES "${_modules}" PARENT_SCOPE)
    set(MARMOT_LEGACY_MODULES "${_legacy_modules}" PARENT_SCOPE)
    set(MARMOT_MODULE_TARGETS "${_targets}" PARENT_SCOPE)
    set(MARMOT_INCLUDE_DIRS "${_include_dirs}" PARENT_SCOPE)
    set(MARMOT_MODULE_LINK_LIBRARIES "${_link_libraries}" PARENT_SCOPE)
endfunction()
