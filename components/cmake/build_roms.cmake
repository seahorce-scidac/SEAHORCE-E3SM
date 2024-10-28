function(build_roms_model)

  message(STATUS "Found ROMS component")

  include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

  # Strip leading space, or else add_compile_definitions will add
  # an empty definition '-D `, which will cause compiler errors
  # string (REGEX REPLACE "^ " "" CPPDEFS "${CPPDEFS}")
  # add_compile_definitions expects a list of definitinos, each
  # one provided without -D. So 1) replace ' ' with ';', and
  # 2) strip '-D'/
  string (REPLACE " " ";" CPPDEFS "${CPPDEFS}")
  string (REPLACE "-D" "" CPPDEFS "${CPPDEFS}")
  add_compile_definitions ("${CPPDEFS}")

  set(ENABLE_ROMS true PARENT_SCOPE)
  set(ROMS_BUILD_DIR ${CMAKE_CURRENT_BUILD_DIR}/roms PARENT_SCOPE)
  add_subdirectory("roms")

endfunction(build_roms_model)
