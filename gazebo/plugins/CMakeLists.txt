cmake_minimum_required(VERSION 3.14)
project(MagneticMapPlugin VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Find Gazebo Harmonic packages ─────────────────────────────────────────────
find_package(gz-cmake3    REQUIRED)
find_package(gz-plugin2   REQUIRED COMPONENTS register)
find_package(gz-sim8      REQUIRED)
find_package(gz-math7     REQUIRED)
find_package(gz-transport13 REQUIRED)
find_package(gz-msgs10    REQUIRED)

# ── Build the plugin shared library ──────────────────────────────────────────
add_library(MagneticMapPlugin SHARED
    MagneticMapPlugin.cc
)

target_compile_options(MagneticMapPlugin PRIVATE
    -Wall -Wextra -O2
)

target_link_libraries(MagneticMapPlugin
    gz-plugin2::register
    gz-sim8::gz-sim8
    gz-math7::gz-math7
    gz-transport13::gz-transport13
    gz-msgs10::gz-msgs10
)

# ── Install destination ───────────────────────────────────────────────────────
# Default: install next to source for easy loading via GZ_SIM_SYSTEM_PLUGIN_PATH
install(TARGETS MagneticMapPlugin
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)

# ── Print post-build instructions ────────────────────────────────────────────
add_custom_command(TARGET MagneticMapPlugin POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "=== Build complete ==="
    COMMAND ${CMAKE_COMMAND} -E echo "Plugin: $<TARGET_FILE:MagneticMapPlugin>"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    COMMAND ${CMAKE_COMMAND} -E echo "To use the plugin:"
    COMMAND ${CMAKE_COMMAND} -E echo "  export GZ_SIM_SYSTEM_PLUGIN_PATH=$<TARGET_FILE_DIR:MagneticMapPlugin>"
    COMMAND ${CMAKE_COMMAND} -E echo "  (or copy the .so to a directory on GZ_SIM_SYSTEM_PLUGIN_PATH)"
    COMMAND ${CMAKE_COMMAND} -E echo ""
    VERBATIM
)

# ── Optional: a quick smoke-test target ──────────────────────────────────────
# Verifies the plugin registers correctly via gz-plugin loader
find_program(GZ_PLUGIN_CLI gz)
if(GZ_PLUGIN_CLI)
    add_custom_target(test_plugin
        COMMAND ${GZ_PLUGIN_CLI} plugin --info --plugin $<TARGET_FILE:MagneticMapPlugin>
        DEPENDS MagneticMapPlugin
        COMMENT "Checking plugin metadata..."
    )
endif()
