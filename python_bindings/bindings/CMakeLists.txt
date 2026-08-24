nanobind_add_module(_core raster_bindings.cpp)

target_link_libraries(_core PRIVATE lsdtt3::runtime)
target_compile_features(_core PRIVATE cxx_std_17)
set_target_properties(_core PROPERTIES
    BUILD_RPATH "$ORIGIN"
    INSTALL_RPATH "$ORIGIN"
)

install(TARGETS _core LIBRARY DESTINATION lsdtt3)
