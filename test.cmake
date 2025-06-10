add_executable("CompressibleNeoHooke_executable" "${CMAKE_CURRENT_LIST_DIR}/test/test.cpp")
target_link_libraries("CompressibleNeoHooke_executable" Marmot)
add_test(NAME "CompressibleNeoHooke" COMMAND "CompressibleNeoHooke_executable")
