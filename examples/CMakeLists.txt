# every example is built and run as a test: it returns non-zero when its result is wrong
set(aerobus_examples
  compensated_horner
  continued_fractions
  custom_taylor
  hermite
  make_polynomial
  modular_arithmetic
  polynomials_over_finite_field
  quotient_ring
  reed_solomon
)
foreach(example IN LISTS aerobus_examples)
  add_executable(example_${example} ${example}.cpp)
  target_link_libraries(example_${example} PRIVATE aerobus aerobus_dev_flags)
  add_test(NAME example.${example} COMMAND example_${example} WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
endforeach()
