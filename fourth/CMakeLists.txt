cmake_minimum_required(VERSION 3.5)

project(fourth VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets)
find_package(Qt6Charts)


add_executable(${PROJECT_NAME}
        main.cpp
        integrals.h
        integrals.cpp
        mainwindow.h mainwindow.cpp mainwindow.ui
)

# set(Boost_USE_STATIC_RUNTIME ON)
# if(Boost_FOUND)
#     target_include_directories(${PROJECT_NAME} PUBLIC ${Boost_INCLUDE_DIR})
#     target_link_libraries(${PROJECT_NAME} PUBLIC ${Boost_LIBRARIES})
#     message(STATUS "Found BOOST ${Boost_VERSION_STRING}")
# else()
#     message(STATUS "WEBSOCKET BOOST Not Found")
# endif()

target_link_libraries(${PROJECT_NAME} PUBLIC Qt${QT_VERSION_MAJOR}::Widgets Qt6::Charts)

include(GNUInstallDirs)
install(TARGETS ${PROJECT_NAME}
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

