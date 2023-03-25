set(HEADERS_CXX_FILES main.h cayley.h sparse.h blades.h dense.h parser.h grade_einsum.h)
set(SOURCES_CXX_FILES main.c cayley.c sparse.c blades.c dense.c parser.c grade_einsum.c)

target_sources(${PROJECT_NAME} PUBLIC
               ${HEADERS_CXX_FILES}
               ${SOURCES_CXX_FILES})

target_include_directories(${PROJECT_NAME} PUBLIC
                           ${CMAKE_CURRENT_SOURCE_DIR})
