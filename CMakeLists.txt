project(FFTW)

# Bring the headers directory
include_directories(inc)

# The static library is in /usr/local/lib
set(PROJECT_LINK_LIBS fftw3.a)

file(GLOB SOURCES "src/*.cpp")

add_executable(fftw  ${SOURCES})
target_link_libraries(fftw ${PROJECT_LINK_LIBS})
