add_library(minpacc minpack.c)
target_include_directories(minpacc
    PUBLIC
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_link_libraries(minpacc PUBLIC m)

add_library(minpacc::minpacc ALIAS minpacc)

install(TARGETS minpacc EXPORT minpaccTargets
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(FILES ${PROJECT_SOURCE_DIR}/include/minpack.h DESTINATION include)
