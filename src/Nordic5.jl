module Nordic5

    using BlockSystems
    using ModelingToolkit
    using Plots
    using DifferentialEquations
    using Distributions
    using ForwardDiff    

    # load component library
    # The components for the Nordic5 system were defined following the publication [1] and the corresponding documentation in the GitHub Repo [2], unless otherwise stated.
    # [1] Dynamic Virtual Power Plant Design for Fast Frequency Reserves: Coordinating Hydro and Wind, J. Björk et. al., 2021, 
    # [2] https://github.com/joakimbjork/Nordic5/blob/main/Nordic5_doc.pdf 
    include(joinpath("components", "auxiliary.jl"))
    include(joinpath("components", "build_nodes.jl"))
    include(joinpath("components", "hydro.jl"))
    include(joinpath("components", "thermal.jl"))
    include(joinpath("components", "wind_turbine.jl"))
    include(joinpath("components", "parameters.jl"))
    include(joinpath("components", "swing.jl"))

    # load utils
    include("utils.jl")
    # load code for ND networks
    include("nd_utils.jl")

    # Utils for ProBeTune: https://github.com/PIK-ICoNe/ProBeTune.jl
    include("probetune_utils.jl")

    # load samplers and loss functions
    include("input_samplers.jl")
    include("loss_functions.jl")

    export NDProblem, make_p_dict, solve_sys_spec_n, pbt_loss, pbt_tuning
    export frequency_metric, single_node_perturbation
    export get_easy_inverter
    export get_hydro_wind_simple, get_thermal_wind_simple
    export get_controlled_swing
    export get_hydro, get_thermal
    export get_hydro_1_p, get_hydro_1_p, get_hydro_2_p, get_hydro_3_p, get_thermal_5_p, get_thermal_4_p, get_wt_simple_p
    export voltage_mag_series, voltage_angle_series, get_u0, frequency_series, has_frequency_series
    export StaticLine
    export get_current, total_current, power_series
end # module
