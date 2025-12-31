enable_language(C CXX)

file(GLOB_RECURSE PIC3D_SOURCE_FILES "./src/*.cpp")
file(GLOB_RECURSE PIC3D_HEADER_FILES "./include/*.h ./include/*.hpp")

source_group("Header Files" FILES ${PIC3D_HEADER_FILES})

add_library(pic3d STATIC "${PIC3D_SOURCE_FILES}" "${PIC3D_HEADER_FILES}")
target_include_directories(pic3d PRIVATE "./include")

include_directories("./include")

# common
target_link_libraries(pic3d common)

# eulerian3d (for MACGrid3d base class)
target_link_libraries(pic3d eulerian3d)

# glfw
target_link_libraries(pic3d "${PROJECT_SOURCE_DIR}/third_party/glfw/lib/glfw3.lib")
