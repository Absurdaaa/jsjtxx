enable_language(C CXX)

file(GLOB_RECURSE PIC2D_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE PIC2D_HEADER_FILES "./include/*.h" "./include/*.hpp")

source_group("Header Files" FILES ${PIC2D_HEADER_FILES})

add_library(pic2d STATIC "${PIC2D_SOURCE_FILES}" "${PIC2D_HEADER_FILES}")
# Need access to both PIC headers and common grid definitions
target_include_directories(pic2d PRIVATE "./include" "${PROJECT_SOURCE_DIR}/common/include")

target_link_libraries(pic2d common)
target_link_libraries(pic2d "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")
