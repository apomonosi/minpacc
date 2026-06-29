add_executable(tester tester.c)
target_include_directories(tester PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(tester PRIVATE minpacc::minpacc)

add_test(NAME minpack_suite COMMAND tester)
