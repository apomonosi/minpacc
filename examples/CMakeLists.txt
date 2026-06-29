foreach(ex hybrd hybrd1 lmder1 lmdif1)
    add_executable(example_${ex} example_${ex}.c)
    target_link_libraries(example_${ex} PRIVATE minpacc::minpacc)
endforeach()
