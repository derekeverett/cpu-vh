set (SOURCES
    CommandLineArguments.cpp
    DynamicalVariables.cpp
    EnergyMomentumTensor.cpp
    EquationOfState.cpp
    FileIO.cpp
    FiniteDifference.cpp
    FluxFunctions.cpp
    FluxLimiter.cpp
    FullyDiscreteKurganovTadmorScheme.cpp
    GlauberModel.cpp
    HalfSiteExtrapolation.cpp
    HydroParameters.cpp
    HydroPlugin.cpp
    HydroWrapper.cpp
    InitialConditionParameters.cpp
    InitialConditions.cpp
    LatticeParameters.cpp
    LocalPropagationSpeed.cpp
    MonteCarloGlauberModel.cpp
    Properties.cpp
    RunHydroWrapper.cpp
    SemiDiscreteKurganovTadmorScheme.cpp
    SourceTerms.cpp
    SpectralRadius.cpp
    )

add_library(cpu-vh_lib SHARED ${SOURCES})
add_executable(cpu-vh ${SOURCES})

target_link_libraries (cpu-vh_lib ${GSL_LIBRARIES} /home/everett.165/libconfig-1.5/lib/.libs)
target_link_libraries (cpu-vh ${GSL_LIBRARIES} /home/everett.165/libconfig-1.5/lib/.libs)

#libconfig not working yet
#target_link_libraries (cpu-vh_lib ${GSL_LIBRARIES} /home/everett.165/libconfig-1.5/lib/.libs)
#target_link_libraries (cpu-vh ${GSL_LIBRARIES} /home/everett.165/libconfig-1.5/lib/.libs)

install(TARGETS cpu-vh DESTINATION ${CMAKE_HOME_DIRECTORY})
