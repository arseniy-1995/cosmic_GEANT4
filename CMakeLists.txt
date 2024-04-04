#----------------------------------------------------------------------------
# Setup the project
#
cmake_minimum_required(VERSION 3.16...3.21)
project(cosmic)

cmake_host_system_information(RESULT OS_NAME QUERY OS_NAME)
#message(STATUS "${OS_NAME}")

cmake_host_system_information(RESULT OS_PLATFORM QUERY OS_PLATFORM)
#message(STATUS "${OS_PLATFORM}")

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if (WITH_GEANT4_UIVIS)
    find_package(Geant4 REQUIRED ui_all vis_all)
else ()
    find_package(Geant4 REQUIRED)
    find_package(Geant4 10.4.0  QUIET OPTIONAL_COMPONENTS gdml qt)
endif ()

#ROOT

find_package(ROOT REQUIRED COMPONENTS RIO Net)
find_package(GSL REQUIRED)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#find_package(CERNLIB REQUIRED)

# ROOT version 6 required
if(ROOT_FOUND)
    STRING(REGEX MATCH "6.*" VERSION6MATCH ${ROOT_VERSION})
    if(NOT VERSION6MATCH)
        message(FATAL_ERROR "P02 requires ROOT 6")
    endif()
endif()

# P01 requires shared libraries
#
if(NOT Geant4_shared_FOUND)
    message(FATAL_ERROR "P02 must use shared libraries")
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
# Setup include directory for this project
#
include(${Geant4_USE_FILE})
#include_directories(${PROJECT_SOURCE_DIR}/include)

include(${ROOT_USE_FILE})
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
        ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
#file(GLOB lib ${PROJECT_SOURCE_DIR}/lib/*.a)

set(CMAKE_CXX_STANDARD 14)

# FOR GENBOS
#####find_package(cernlib REQUIRED)

#include(${CERNLIB_INCLUDE_DIR})
#include(${CERNLIB_LIBRARIES})

#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/linux_arm64/libgenbos.a -lkernlib -lphtools -lmathlib -lgfortran -lm)
#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.a -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib -lkernlib -lphtools -lmathlib -lgfortran -lm )

#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.a -lkernlib -lphtools -lmathlib /opt/local/lib/gcc12/libgfortran.a -lm -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib)

##set(GENBOSLIBS $(PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.a -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib -lkernlib -lphtools -|mathlib -lgfortran -lm)
##set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.a -lkernlib -lphtools -lmathlib /opt/local/lib/gcc12/libgfortran.a -lm -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib)
##set(GENBOSLIBS $(PROJECT_SOURCE_DIR}/lib/linux_arm64/libgenbos.a -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib)

#####LINUX x86_64 OK

if (${OS_NAME} MATCHES "Linux")
    if (${OS_PLATFORM} MATCHES "x86_64")
        set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/linux_x64/libgenbos.a -L/usr/local/lib/ -lphtools -lmathlib -lgfortran -lkernlib -lnsl -lcrypt -ldl -lphtools -lmathlib -lgeant321 -lpacklib -lgfortran -lm -lnsl -lphtools)
    endif ()
endif ()

#####

#####LINUX ARM OK

if (${OS_NAME} MATCHES "Linux")
    if (${OS_PLATFORM} MATCHES "aarch64")
        set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/linux_arm64/libgenbos.a -lphtools -lmathlib -lgfortran -lkernlib -lnsl -lcrypt -ldl -lphtools -lmathllb -lgeant321 -lpacklib -lgfortran -lm -lns -lphtools)
    endif ()
endif ()


#####

##### MAC ARM OK

if (${OS_NAME} MATCHES "macOS")
    if (${OS_PLATFORM} MATCHES "arm64")
        set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.a -L/usr/local/lib/ -lphtools -lmathlib -lgfortran -lkernlib -lm -L/usr/local/lib/cernlib/2023/lib)
    endif ()
endif ()


#####

##### MAC x86_64 OK

if (${OS_NAME} MATCHES "macOS")
    if (${OS_PLATFORM} MATCHES "x86_64")
        set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_x64/libgenbos.a -L/usr/local/lib/ -lphtools -lmathlib -lgfortran -lkernlib -lm -L/usr/local/lib/cernlib/2023/lib)
    endif ()
endif ()


#####

#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/linux_arm64/libgenbos.a -L/opt/homebrew/Cellar/cernlib/2023.08.14.0/lib/cernlib/2023/lib)

#  /opt/homebrew/Cellar/gcc/13.2.0/lib/gcc/current/libgfortran.a

#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/linux_x64/libgenbos.a -lkernlib -lphtools -lmathlib -lgfortran -lm )
#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/libgenbos.a -lgfortran -lm)
#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/libgenbos.a -lgfortran -lm)

#set(GENBOSLIBS ${PROJECT_SOURCE_DIR}/lib/mac_arm64/libgenbos.so ${CERNLIB_LIBRARIES} )


#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
add_executable(cosmic cosmic.cc ${sources} ${headers})
#target_link_libraries(cosmic PUBLIC ROOT::RIO ROOT::Net ${ROOT_LIBRARIES} ${Geant4_LIBRARIES})

# FOR GENBOS
target_link_libraries(cosmic PUBLIC ROOT::RIO ROOT::Net ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} ${GENBOSLIBS})

#target_link_libraries(cosmic PUBLIC ${Geant4_LIBRARIES} ${GENBOSLIBS})

#target_link_libraries(cosmic ${Geant4_LIBRARIES})
#----------------------------------------------------------------------------
# Copy all scripts to the build directory, i.e. the directory in which we
# build Cosmic. This is so that we can run the executable directly because it
# relies on these scripts being in the current working directory.
#
set(COSMIC_SCRIPTS
        cosmic.out
        cosmic.in
        gui.mac
        init_vis.mac
        plotHisto.C
        plotNtuple.C
        run1.mac
        run2.mac
        GENBOS.mac
        Proton.mac
        Neutron.mac
        ProtonNeutron.mac
        Gamma.mac
        Cosmic.mac
        LQ.mac
        vis.mac
        )

foreach (_script ${COSMIC_SCRIPTS})
    configure_file(
            ${PROJECT_SOURCE_DIR}/${_script}
            ${PROJECT_BINARY_DIR}/${_script}
            COPYONLY
    )
endforeach ()


# Add program to the project targets
# (this avoids the need of typing the program name after make)
#
add_custom_target(cosmic_ DEPENDS cosmic)
#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS cosmic DESTINATION bin)
