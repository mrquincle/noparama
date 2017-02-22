
file(GLOB_RECURSE HEADER_FILES *.h *.hpp)

set(INCLUDE_DIRS_FLAGS "")
get_directory_property(INCLUDE_DIRS INCLUDE_DIRECTORIES)
foreach(D ${INCLUDE_DIRS})
    set(INCLUDE_DIRS_FLAGS "${INCLUDE_DIRS_FLAGS} -I${D}")
endforeach()

set(ALL_INCLUDES ${INCLUDE_DIRS_FLAGS} ${FLAGS})

set(header_commands_file ${PROJECT_BINARY_DIR}/header/compile_commands.json)

file(WRITE ${header_commands_file} "[\n")
foreach(F ${HEADER_FILES})
    file(APPEND ${header_commands_file} "{\n")
    file(APPEND ${header_commands_file} "  \"directory\": \"${CMAKE_BINARY_DIR}/header\",\n")
    file(APPEND ${header_commands_file} "  \"command\": \"/usr/bin/c++ ${ALL_INCLUDES}\",\n")
    file(APPEND ${header_commands_file} "  \"file\": \"${F}\"\n")
    file(APPEND ${header_commands_file} "},\n")
endforeach()
file(APPEND ${header_commands_file} "]")

