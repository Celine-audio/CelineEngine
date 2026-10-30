# CelineEngine -- the circuit simulator shared by Céline and Céline Player.
#
# A plugin brings it in with three lines, after its SharedCode target exists:
#
#     set(CELINE_ENGINE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/modules/CelineEngine" CACHE PATH "...")
#     include("${CELINE_ENGINE_DIR}/CelineEngine.cmake")
#     target_link_libraries(SharedCode INTERFACE CelineEngine)
#
# An INTERFACE library, like CelineUI: its sources compile inside each plugin,
# against that plugin's JUCE and with its flags, rather than into a library of
# their own. The engine needs nothing from the plugin, but building it the same
# way as everything else in SharedCode keeps one set of warnings and one
# optimisation level over the code that matters most for CPU.

include_guard(GLOBAL)

add_library(CelineEngine INTERFACE)

file(GLOB_RECURSE CelineEngineSources CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_LIST_DIR}/CelineEngine/*.cpp"
    "${CMAKE_CURRENT_LIST_DIR}/CelineEngine/*.h")

target_sources(CelineEngine INTERFACE ${CelineEngineSources})

# The repository root, not the source folder, so a plugin writes
# <CelineEngine/Engine.h> and every include says which repository the file
# lives in.
target_include_directories(CelineEngine INTERFACE "${CMAKE_CURRENT_LIST_DIR}")
