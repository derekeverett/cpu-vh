set (SOURCES
    CommandLineArguments.cpp
    DynamicalSources.cpp
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
    Vorticity.cpp
    )

add_library(cpu-vh_lib SHARED ${SOURCES})
add_executable(cpu-vh ${SOURCES})


target_link_libraries (cpu-vh_lib ${GSL_LIBRARIES})
target_link_libraries (cpu-vh ${GSL_LIBRARIES})

#target_link_libraries (cpu-vh_lib ${GSL_LIBRARIES} ${Config_LIBRARIES})
#target_link_libraries (cpu-vh ${GSL_LIBRARIES} ${Config_LIBRARIES})

install(TARGETS cpu-vh DESTINATION ${CMAKE_HOME_DIRECTORY})
