cmake_minimum_required(VERSION 3.5)

project(task11 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(Boost)
add_executable(${PROJECT_NAME}
    main.cpp)

target_link_libraries(${PROJECT_NAME} Boost::boost Boost::headers)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)
