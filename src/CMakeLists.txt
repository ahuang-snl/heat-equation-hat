INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
SET(main_SOURCES
  main.cpp
  Step01_LinearFunction.cpp
  Step01_SinXSinYFunction.cpp
  )

INCLUDE_DIRECTORIES ( ./ ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

LINK_DIRECTORIES (${Trilinos_LIBRARY_DIRS} ${Trilinos_TPL_LIBRARY_DIRS})

ADD_EXECUTABLE(
  step01.exe
  ${main_SOURCES}
  )

TARGET_LINK_LIBRARIES(step01.exe ${Trilinos_LIBRARIES}
${Trilinos_TPL_LIBRARIES})


FILE(COPY ./input.xml DESTINATION .)
