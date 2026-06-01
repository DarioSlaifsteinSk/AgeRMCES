cd(@__DIR__)

# Import all the necessary packages
using JuMP, InfiniteOpt, KNITRO, EMSmodule, LiiBRA
using LinearAlgebra, Distributions, Statistics, Parameters, Random, Revise, Test
using DataFrames, LaTeXStrings, Printf, JSON3, Makie, CairoMakie, GLMakie
import EMSmodule.@unpack_JinAgingParams
Random.seed!(1234);
# nEV=1:2;
# includet("../functions/makeEMSobjs_perp.jl") # data structures
includet("../functions/build_data.jl") # data structures
# includet("../functions/makeEMSplots.jl") # data structures
# includet("../functions/EMSfns_perp.jl") # General functions
# includet("../functions/ESSfns_perp.jl") # BESS functions
# includet("../functions/testEMS_perp.jl") # Test functions
# includet("../functions/simEMS.jl") # Plotting functions

## Define all the necessary functions
function build_model(optimizer, # model optimizer
    sets::modelSettings, # number of EVs, discrete time supports, etc.
    data::Dict, # information for the model (mainly parameters), these are the devices (EV, BESS, PV, etc.), the costs (interests, capex, etc) and the revenues
    preRes::Dict, # previous results to warm start the model,
    )
    # Sets
    tend=sets.dTime[end]
    t0=sets.dTime[1]; # initial time, can´t divide by 0
    model = InfiniteModel(optimizer) # create model
    tunerpath = "../src/tunerfile-explore.opt" # path to the tuner file
    optpath = "../src/optfile.opt" # path to the tuner file


    # Optimizer attributes
    # KNITRO attributes
    # set_optimizer_attributes(model, "outdir"=> "..\\data\\output\\Run21",
    #                         "outmode"=>2,
    #                         "outname"=>"run_14_0_0.log",
    #                         "outlev"=>2,
    #                         )
    set_attributes(model,
                # "tuner"=>1,
                "scale"=>1,
                "presolve"=>1,
                "outlev"=>3,
                # "ms_enable"=>1,
                # "ms_maxsolves"=>20,
                "nlp_algorithm"=>3,
                "numthreads"=>4,
                "ms_numthreads"=>4,                
                )
    set_time_limit_sec(model, 3*60*60);

    # define cont t
    @infinite_parameter(model, t ∈ [t0, tend], supports = collect(sets.dTime),
                        derivative_method = FiniteDifference(Backward()))
                        # derivative_method = OrthogonalCollocation(3)) # nope
                        # derivative_method = FiniteDifference(Forward(),true))

    # Add devices
    # Electrical
    spv!(model, data);
    bess!(model, sets, data);
    ev!(model, sets, data);
    # Thermal
    st!(model, data);
    heatpump!(model, data);
    tess!(model, data);
    gridThermal!(model, sets, data); # thermal balance

    # grids and balances
    gridConn!(model, data); # grid connection
    pei!(model, sets, data); # electrical balance
    # Add objective function
    costFunction!(model, sets, data)
    # check if the previous solution is not empty
    !isempty(preRes) ? set_warm_start(model, preRes) : nothing;
    set_silent(model)
    return model
end;

function runDayAhead(sDA::modelSettings,
    sCT::modelSettings,
    inputW_DA::Dict, # for opt
    inputW_CT::Dict, # for sim
    preSol::Dict,
    manager,
    )
    modelDA = build_model(() -> KNITRO.Optimizer(; license_manager = manager), sDA, inputW_DA, preSol)
    set_attributes(modelDA,
        "ms_enable"=>1,
        "ms_maxsolves"=>5,
        "nlp_algorithm"=>1, # interior point for day-ahead, active sets for MPC
        )
    set_silent(modelDA)
    optimize!(modelDA)
    # save results, primal
    results = getResultsFull(modelDA)
    # free the license
    KNITRO.KN_release_license(manager)
    # Upsample the day-ahead results to match the MPC time resolution
    for key ∈ keys(results)
        try
            results[key] = repeat(results[key], inner = Int(sDA.Δt/sCT.Δt))
        catch e
        end
    end
    DtDA = supports(modelDA[:t])
    results["t"] = collect(DtDA[1]:(sCT.Δt*3600):(DtDA[end]+3600-sCT.Δt*3600))
    results, ctrlRes = handleInfeasible!(results, preSol, 1, sDA, CT_MPC()) # CHECK
    simTransitionFun!(results, inputW_CT, sCT, CT_MPC())
    sync_ess!(inputW_CT, inputW_DA[:"BESS"].AgingParameters)
    sync_S0!(inputW_DA, inputW_CT, sDA, inputW_DA[:"BESS"].AgingParameters)
    return results, ctrlRes, inputW_DA, inputW_CT, manager
end

function solvePolicies(optimizer, # model optimizer
    sets::modelSettings, # number of EVs, discrete time supports, etc.
    data::Dict, # information for the model (mainly parameters), these are the devices (EV, BESS, PV, etc.), the costs (interests, capex, etc) and the revenues
    preRes::Dict, # previous results to warm start the model
    )
    # Sets
    tend=sets.dTime[end]
    t0=sets.dTime[1]; # initial time, can´t divide by 0
    model = InfiniteModel(optimizer) # create model
    tunerpath = "../src/tunerfile-explore.opt" # path to the tuner file
    optpath = "../src/optfile.opt" # path to the tuner file


    # Optimizer attributes
    # KNITRO attributes
    # set_optimizer_attributes(model, "outdir"=> "..\\data\\output\\Run21",
    #                         "outmode"=>2,
    #                         "outname"=>"run_14_0_0.log",
    #                         "outlev"=>2,
    #                         )
    set_optimizer_attributes(model,
                            "tuner"=>1,
                            "scale"=>1,
                            "presolve"=>1,
                            "outlev"=>3,
                            # "ms_enable"=>1,
                            # "ms_maxsolves"=>20,
                            "algorithm"=>3,
                            "numthreads"=>8,
                            # "tuner_file"=>tunerpath,
                            # "mip_maxnodes" => 3500,
                            # "convex"=>1,
                            # "opttol"=>1e-3,
                            # "feastol"=>1e-3,
                            # "mip_opt_gap_abs"=>1e-2,
                            # "mip_multistart"=>1,
                            # # options
                            # "hessopt"=>1,
                            # "hessian_no_f"=>1,
                            # "mip_method" =>	1,
                            # "mip_nodealg" => 1,
                            # "mip_selectrule" =>	2,
                            # "mip_branchrule" =>	2,
                            # "mip_heuristic_strategy" =>	2,
                            # "mip_heuristic_feaspump" =>1,
                            # "mip_heuristic_localsearch" =>	1,
                            # "mip_heuristic_mpec" => 1,
                            )
    set_time_limit_sec(model, 6*60*60);

    # define cont t
    @infinite_parameter(model, t ∈ [t0, tend], supports = collect(sets.dTime),
                        derivative_method = FiniteDifference(Backward()))
                        # derivative_method = OrthogonalCollocation(3)) # nope
                        # derivative_method = FiniteDifference(Forward(),true))

    # Add devices
    # Electrical
    spv!(model, data);
    bess!(model, sets, data);
    ev!(model, sets, data);
    # Thermal
    st!(model, data);
    heatpump!(model, data);
    tess!(model, data);
    gridThermal!(model, sets, data); # thermal balance

    # grids and balances
    gridConn!(model, data); # grid connection
    pei!(model, sets, data); # electrical balance
    # Add objective function
    costFunction!(model, sets, data)
    # check if the previous solution is not empty
    !isempty(preRes) ? set_warm_start(model, preRes) : nothing;
    # Solve model
    optimize!(model)
    results = getResults_perp(model)
    ddict, sdpdict = getDuals(model)
    MOI.empty!(InfiniteOpt.backend(model))
    return results, ddict, sdpdict;
end;

function handleInfeasible!(results::Dict, preSol::Dict, ts, s, typeOpt::CT_MPC)
    # this function is used to handle infeasible solutions
    # it is called in the solvePolicies() function
    # If an MPC solution is infeasible implement the previous solution 
    shift = 1

    println("\r Checking the status of the last set of results")
    if results[:"status"] == LOCALLY_INFEASIBLE
        # discard the last set of results and implement the second step of the previous solution at ts-1
        println("\r Infeasible solution at step $ts, implementing the second step of the previous solution at step $(ts-1)")
        results=copy(preSol);
        # the first setpoints from results[ts-1] were already implemented in the previous step, 
        # so we need to shift the setpoints one step, notice that the end of the time window won't move! 
        for k in keys(results)
            # Check if the key requires special handling and skip the "status" key
            if k == :"status"
                results[k] = LOCALLY_INFEASIBLE;
            elseif k == :"γf"
                if s.nEV !== 1
                    results[k] = [results[k][n][2:end] for n ∈ 1:s.nEV];
                else
                    results[k] = results[k][2:end];
                end
            else
                length(results[k]) == 1 ? continue : results[k] = results[k][2:end];
                # results[k] = results[k][2:end];
            end
        end
    else
        println("\r Locally feasible solution at step $ts")
    end
    # update the controller setpoints so that for next iteration this is stored somewhere
    preSol = copy(results);
    return results, preSol;
end

function build_pbrom_sim_cell!(data::Dict, s::modelSettings)
    nEV = s.nEV; cellID = s.cellID;
    ## BESS MODEL
    PowerLim = [-12.5, 12.5]; P0 = 0; # Power limits and initial condition [kW]. 
    SoCLim = [0.15, 1.]; SoC0 = 0.4; # SoC limits and initial condition [p.u.].
    ηC = 0.95; # Charger efficiency [p.u.].
    bessOCV=JSON3.read(open("../data/input/cell_models/$cellID-modelocv.json", "r"),
                Dict{String, Vector{Float64}});
    # Dynamic ECM model, with np=1
    bessDYN=JSON3.read(open("../data/input/cell_models/$cellID-modeldyn-no-hys.json", "r"),
                Dict{String, Union{Vector{Float64}, Vector}});
    # Bucket model info
    Q0 = bessOCV["OCVQ"]; Q0 = mean(Q0); # [Ah/cell]
    η = bessOCV["OCVeta"]; η = mean(η); # [p.u.]
    # ECM info
    T=bessDYN["temps"]; T=T[1,:]; # [C]
    R0Param=bessDYN["R0Param"]; # [Ohm]
    RParam=bessDYN["RParam"]; # [Ohm]
    RCParam=bessDYN["RCParam"]; # τ [s]
    Tind25=findall(x->x==25, T)
    # State of Health
    SoHQ = 1; SoHR0 = R0Param[Tind25][1];
    if cellID == "SYNSANYO"
        # lets try a smaller BESS
        Np = 7; Ns = 100; # Branches in parallel and series cells per branch.
        ocv_params=OCVlinearPerfParams()
        vLim = [2.8, 4.2]
        # Ageing submodel
        aging_params=JinAgingParams();
        perf_params = CIDRAPBROMPerfParams()
        perf_params.Cell.Neg.RFilm = R0Param[Tind25][1];
    elseif cellID == "A123"
        Np = 25; Ns = 110; # Branches in parallel and series cells per branch.
        # first option, 
        # ocv_params.ocvLine = [2.5, 1.1]
        # second option from 20% to 95% SoC
        ocv_params=OCVlinearPerfParams(ocvLine = [3.2, 0.2105])
        vLim = [2.0, 3.6]
        # Ageing submodel
        # aging_params=empAgingParams();
        Cell = Construct(cellID);
        perf_params = CIDRAPBROMPerfParams(Cell=Cell)
        aging_params=JinAgingParams(Rs = Cell.Neg.Rs,
                            An = Cell.Const.CC_A,
                            Ln = Cell.Neg.L,
                            z100p = Cell.Neg.θ_100,
                            z0p = Cell.Neg.θ_0,
                            # εₑ0 = Cell.Neg.ϵ_e,
                            t⁺₀ = Cell.Const.t_plus,
                            # DeRef= Cell.Neg.De,
                            ce_avg = Cell.Const.ce0,
                            # ce_max=Cell.Const.ce0,
                            σn = Cell.Neg.σ,
                            εₛ = Cell.Neg.ϵ_s, # check
                            );
    end
    initVal=494.246; # Cost info [USD/kWh]
    # General info definition
    gen_params = Generic(PowerLim, P0, SoCLim, SoC0, 6., Q0, SoHQ, SoHR0, Np, Ns, η, ocv_params, vLim, ηC, initVal);
    
    # wrap everything in a BESSData type
    battModel = BESSData(gen_params, perf_params, aging_params, cellID);

    ## EV MODEL
    # Battery pack definition
    PowerLim = [-12.5, 12.5] #check
    # The battery pack has to be around 400Vdc and 50kWh
    if cellID == "SYNSANYO"
        Ns = 100; Np = 25;
    elseif cellID == "A123"
        Ns = 110; Np = 61;
    end
    P0=[0, 0]; # initial
    Q0 = Q0.*ones(nEV); # [Ah/cell]
    SoC0 = [0.6, 0.8]; # initial
    vt0=[gen_params.OCVParam.ocvLine[1]+gen_params.OCVParam.ocvLine[2].*SoC0[n] for n ∈ 1:nEV]
    evModel = Vector{EVData}(undef, nEV);
    for n ∈ 1:nEV
        gen_params = Generic(PowerLim, P0[n], SoCLim, SoC0[n], 6., Q0[n], SoHQ, SoHR0, Np, Ns, η, ocv_params, vLim, ηC, initVal)
        batteryPack = BESSData(gen_params, perf_params, aging_params, cellID);
        # wrap everything in a EVData type
        evModel[n] = EVData(batteryPack, data["EV"][n].driveInfo)
    end

    return merge!(data, Dict(:"BESSpbrom" => battModel,
                        :"EVpbrom" => evModel));
end

function age_cell!(data::Dict, α::Float64)
    nEV = length(data[:"EV"]);
    β = (1-α)/2; # SoHR0 [p.u.]
    # Capacity fade
    data[:"BESS"].GenInfo.SoHQ = α;
    data[:"BESS"].GenInfo.initQ *= α;
    data[:"BESS"].AgingParameters.z100p *= α;
    # Power fade
    
    data[:"BESS"].GenInfo.SoHR0 *= (1+β);
    @unpack_JinAgingParams data[:"BESS"].AgingParameters;
    data[:"BESS"].AgingParameters.εₑ0 += data[:"BESS"].AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * β * data[:"BESS"].PerfParameters.R0Param[1] / An *1e-4;
    data[:"BESS"].AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"BESS"].PerfParameters.R0Param[1] / An *1e-4;
    data[:"BESS"].PerfParameters.R0Param *= (1+β);
    data[:"BESS"].AgingParameters.initT = 2 * 365 * 24 * 3600; # 2 years
    # Modify the EV model
    for n in 1:nEV
        # Capacity fade
        data[:"EV"][n].carBatteryPack.GenInfo.SoHQ = α
        data[:"EV"][n].carBatteryPack.GenInfo.initQ *= α
        data[:"EV"][n].carBatteryPack.AgingParameters.z100p *= α
        # Power fade
        data[:"EV"][n].carBatteryPack.GenInfo.SoHR0 *= (1+β)
        @unpack_JinAgingParams data[:"EV"][n].carBatteryPack.AgingParameters;
        data[:"EV"][n].carBatteryPack.AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"EV"][n].carBatteryPack.PerfParameters.R0Param[1] / An *1e-4;
        data[:"EV"][n].carBatteryPack.AgingParameters.εₑ0 += data[:"EV"][n].carBatteryPack.AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * 
                                                        β * data[:"EV"][n].carBatteryPack.PerfParameters.R0Param[1] / An *1e-4;
        data[:"EV"][n].carBatteryPack.PerfParameters.R0Param *= (1+β)
        data[:"EV"][n].carBatteryPack.AgingParameters.initT = 2 * 365 * 24 * 3600.; # 2 years
    end
    return data;
end

function age_cell_pbrom!(data::Dict, α::Float64)
    nEV = length(data[:"EV"]);
    β = (1-α)/2; # SoHR0 [p.u.]
    # Capacity fade
    data[:"BESSpbrom"].GenInfo.SoHQ = α;
    data[:"BESSpbrom"].GenInfo.initQ *= α;
    data[:"BESSpbrom"].AgingParameters.z100p *= α;
    # Power fade
    data[:"BESSpbrom"].GenInfo.SoHR0 *= (1+β);
    @unpack_JinAgingParams data[:"BESSpbrom"].AgingParameters;
    # data[:"BESSpbrom"].AgingParameters.εₑ0 += data[:"BESSpbrom"].AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * β * data[:"BESSpbrom"].PerfParameters.R0Param[1] / An *1e-4;
    # data[:"BESSpbrom"].AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"BESSpbrom"].PerfParameters.R0Param[1] / An *1e-4;
    data[:"BESSpbrom"].AgingParameters.εₑ0 += data[:"BESSpbrom"].AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * β * data[:"BESSpbrom"].GenInfo.SoHR0 / An *1e-4;
    data[:"BESSpbrom"].AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"BESSpbrom"].GenInfo.SoHR0 / An *1e-4;
    data[:"BESSpbrom"].AgingParameters.initT = 2 * 365 * 24 * 3600; # 2 years
    # Modify the EV model
    for n in 1:nEV
        # Capacity fade
        data[:"EVpbrom"][n].carBatteryPack.GenInfo.SoHQ = α
        data[:"EVpbrom"][n].carBatteryPack.GenInfo.initQ *= α
        data[:"EVpbrom"][n].carBatteryPack.AgingParameters.z100p *= α
        # Power fade
        data[:"EV"][n].carBatteryPack.GenInfo.SoHR0 *= (1+β)
        @unpack_JinAgingParams data[:"EV"][n].carBatteryPack.AgingParameters;
        # data[:"EVpbrom"][n].carBatteryPack.AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"EV"][n].carBatteryPack.PerfParameters.R0Param[1] / An *1e-4;
        # data[:"EVpbrom"][n].carBatteryPack.AgingParameters.εₑ0 += data[:"EV"][n].carBatteryPack.AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * 
        #                                                 β * data[:"EV"][n].carBatteryPack.PerfParameters.R0Param[1] / An *1e-4;
        data[:"EVpbrom"][n].carBatteryPack.AgingParameters.δSEI0 += (κref .* εₑ0^brug)  / εₛ * β * data[:"EV"][n].carBatteryPack.GenInfo.SoHR0 / An *1e-4;
        data[:"EVpbrom"][n].carBatteryPack.AgingParameters.εₑ0 += data[:"EV"][n].carBatteryPack.AgingParameters.as * (κref .* εₑ0^brug)  / εₛ * 
                                                        β * data[:"EV"][n].carBatteryPack.GenInfo.SoHR0 / An *1e-4;
        data[:"EVpbrom"][n].carBatteryPack.AgingParameters.initT = 2 * 365 * 24 * 3600.; # 2 years
    end
    return data;
end

function rollingHorizon(sDA, sCT, typeOpt::CT_MPC, α::Float64=1.) # MPC or day-ahead
    ## Rolling Horizon Simulation.
    # This function simulates the EMS for a given number of steps.
    # The EMS is initialized at time t0 and then it is solved for a time window of Tw hours.
    # Then, the EMS is solved again for the next Tw hours, but this time the initial conditions are the ones obtained from the previous solution.
    # This process is repeated until the number of steps is reached.
    # The output of this function are:
    # - results::Vector{Dict}, where each dictionary contains the results of the EMS for each time window.
    # - rhDict::Dict, which is a dictionary with the concatenated results of the EMS for each time window.
    manager = KNITRO.LMcontext()
    # allocate memory
    steps = sCT.steps
    fsCT = 1/sCT.Δt # sampling frequency, per hour
    steps_per_day = Int(24 * fsCT)
    total_days = sDA.steps
    res_DA = Vector{Dict}(undef, total_days)
    results = Vector{Dict}(undef, steps*total_days)
    ctrlRes = Vector{Dict}(undef, steps*total_days)
    # Build the data dictionary
    inputW_DA=build_data(;nEV = sDA.nEV,
                        season = sDA.season,
                        profType = sDA.profType,
                        loadType = sDA.loadType,
                        year = sDA.year,
                        fs = Int(1/sDA.Δt),
                        cellID = sDA.cellID,
                        );

    inputW_CT=build_data(;nEV = sCT.nEV,
                        season = sCT.season,
                        profType = sCT.profType,
                        loadType = sCT.loadType,
                        year = sCT.year,
                        fs = Int(1/sCT.Δt),
                        cellID = sCT.cellID,
                        );
    start_data!(inputW_CT, inputW_DA, sCT)
    α < 1. ? age_cell_pbrom!.([inputW_DA, inputW_CT], α) : nothing;
    # Initial solution
    for day in 1:total_days
        # Day-ahead schedule at midnight
        iDA = (day-1)*steps_per_day + 1
        sDA.t0 = (day - 1.) * 24. + sDA.Δt
        sDA.tend = sDA.t0 + sDA.Tw
        sDA.dTime = collect(sDA.t0:sDA.Δt:sDA.tend) * 3600
        # Solve the day-ahead problem
        iDA == 1 ? preSol = Dict() : preSol = ctrlRes[iDA-1] # CHECK THIS
        try
            results[iDA]["compTime"] = @elapsed begin
                results[iDA], ctrlRes[iDA], inputW_DA, inputW_CT, manager = runDayAhead(sDA, sCT, inputW_DA, inputW_CT, preSol, manager)
            end
        catch e
            println("Error at day $day")
            showerror(stdout, e, catch_backtrace())
            return results, res_DA, ctrlRes, inputW_DA, inputW_CT, sDA, sCT
        end
        
        # save the day-ahead results in another dictionary
        res_DA[day] = copy(results[iDA])
        # MPC for the rest of the day
        sCT.t0 = sDA.t0 + sCT.Δt
        sCT.tend = sCT.t0 + sCT.Tw
        sCT.dTime = collect(sCT.t0:sCT.Δt:sCT.tend) * 3600
        # Solve the MPC problem
        for ts in 2:steps
            global_index = (day-1)*steps_per_day + ts
            try
                results[global_index]["compTime"] = @elapsed begin
                    # modelCT = build_model(() -> KNITRO.Optimizer(; license_manager = manager),
                    #      sCT, inputW_CT, ctrlRes[global_index-1], res_DA[day])
                    modelCT = build_model(() -> KNITRO.Optimizer(; license_manager = manager),
                         sCT, inputW_CT, ctrlRes[global_index-1])
                    set_silent(modelCT)
                    optimize!(modelCT)
                    # save results, primal
                    results[global_index] = getResultsFull(modelCT)
                    # free the license
                    KNITRO.KN_release_license(manager)
                    results[global_index], ctrlRes[global_index] = handleInfeasible!(results[global_index],
                            ctrlRes[global_index-1], global_index, sCT, CT_MPC())
                    simTransitionFun!(results[global_index], inputW_CT, sCT, CT_MPC())
                    sync_ess!(inputW_CT, inputW_CT[:"BESS"].AgingParameters)
                    sCT.dTime .+= (sCT.Δt * 3600.0)
                end
                println("\r Step $global_index out of $(sCT.steps) done!")
            catch e
                println("Error at step $global_index")
                showerror(stdout, e, catch_backtrace())
                return results, res_DA, ctrlRes, inputW_DA, inputW_CT, sDA, sCT
            end
        end
        # Sync initial conditions
        sync_S0!(inputW_DA, inputW_CT, sDA, inputW_DA[:"BESS"].AgingParameters)
    end
    return results, res_DA, ctrlRes, inputW_DA, inputW_CT, sDA, sCT
end

function rollingHorizon(s;
    typeOpt::String="MPC",
    α::Float64=1. # MPC or day-ahead
    ) # MPC or day-ahead
    ## Rolling Horizon Simulation.
    # This function simulates the EMS for a given number of steps.
    # The EMS is initialized at time t0 and then it is solved for a time window of Tw hours.
    # Then, the EMS is solved again for the next Tw hours, but this time the initial conditions are the ones obtained from the previous solution.
    # This process is repeated until the number of steps is reached.
    # The output of this function are:
    # - results::Vector{Dict}, where each dictionary contains the results of the EMS for each time window.
    # - rhDict::Dict, which is a dictionary with the concatenated results of the EMS for each time window.
    @assert typeOpt ∈ ["MPC", "day-ahead"] "typeOpt must be either MPC or DA"
    if typeOpt == "day-ahead"
        @assert (s.Tw + s.Δt) % 24 == 0 "Tw must be a multiple of fs times 24"
    end

    steps = s.steps;
    Dt = s.dTime[1]:(s.dTime[2]-s.dTime[1]):s.dTime[end];
    # allocate memory
    results=Vector{Dict}(undef, steps); controllerRes=Vector{Dict}(undef, steps);
    ddict=Vector{Dict}(undef, steps); sdpdict=Vector{Dict}(undef, steps);
    # Initialize an array to store the times
    times = Vector{Float64}(undef, steps)

    # Build the data dictionary
    EMSData=build_data(;nEV = s.nEV,
                        season = s.season,
                        profType = s.profType,
                        loadType = s.loadType,
                        year = s.year,
                        cellID = s.cellID,
                        );
    EMSData[:"grid"].λ /= 1000; # convert to $/kWh
    EMSData[:"grid"].λ /= 3600; # convert to $/kWs
    # Modify the TESS init condition
    EMSData[:"TESS"].SoC0 = 0.4;
    # Modify the BESS perf. model
    EMSData[:"BESS"].GenInfo.SoC0 = 0.6;
    EMSData[:"BESS"].GenInfo.SoCLim[1] = 0.15;
    # EMSData[:"BESS"].PerfParameters = bucketPerfParams();
    # # Modify the EV model
    [EMSData[:"EV"][n].carBatteryPack.GenInfo.SoCLim[1] = 0.15 for n ∈ 1:s.nEV]
    # [EMSData[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams() for n ∈ 1:s.nEV];
    
    build_pbrom_sim_cell!(EMSData, s); # add simulation cell models

    α < 1. ? age_cell_pbrom!(EMSData, α) : nothing;

    # Initial solution
    times[1] = @elapsed begin
        results[1], ddict[1], sdpdict[1] = solvePolicies(KNITRO.Optimizer, s, EMSData, Dict());
        controllerRes[1] = copy(results[1]);
        # Check the status of the last set of results
        handleInfeasible!(results, controllerRes, 1, s; typeOpt=typeOpt);
        simTransitionFun_cs3!(results[1], EMSData, s; typeOpt = typeOpt)
        # move time window
        Dt = Dt .+ 24 * 3600.0
        s.dTime=collect(Dt);
    end

    for ts in 2:steps
        try
            times[ts] = @elapsed begin
                # build+solve model
                results[ts], ddict[ts], sdpdict[ts] =solvePolicies(KNITRO.Optimizer, s, EMSData, Dict());
                controllerRes[ts] = copy(results[ts]);
                # Check the status of the last set of results
                handleInfeasible!(results, controllerRes, ts, s; typeOpt=typeOpt);
                simTransitionFun_cs3!(results[ts], EMSData, s; typeOpt = typeOpt)
                # move time window
                Dt = Dt .+ 24 * 3600.0
                s.dTime=collect(Dt);
            end
            # Base.GC.gc(); # collect garbage
            println("\r Step $ts out of $steps done!")
        catch e
            println("Error at step $ts")
            showerror(stdout, e, catch_backtrace())
            return results, controllerRes, EMSData, s, ddict, sdpdict;
        end
    end
    [results[ts]["compTime"]=times[ts] for ts in 1:steps]
    return results, controllerRes, EMSData, s, ddict, sdpdict;
end

# define settings
resultsDict = Dict(:α=>[], :results=>[], :controller => []);
# Run the rolling horizon simulation
for α ∈ 1.0:-0.1:0.9
    Random.seed!(1234);
    # Wloss=0.1
    Wgrid = 1; Wloss = 0.; W=[Wgrid 1000 Wloss 1000];
    sDA=modelSettings(nEV=1, t0=0,Tw=48-1, Δt=1, steps=29, costWeights=W, season="winter",
                profType="monthly", loadType="GV", year=2023, cellID="SYNSANYO");
    sCT=modelSettings(nEV=1, t0=1/4,Tw=24-1/4, Δt=1/4, steps=96, costWeights=W, season="winter",
                profType="monthly",loadType="GV", year=2023, cellID="SYNSANYO");
    # Define the weights
    results, res_DA, ctrlRes, ~, ~, ~, ~  = rollingHorizon(sDA, sCT, CT_MPC(), α);
    # associate the results with the weights
    push!(resultsDict[:α], α);
    push!(resultsDict[:results], concatResultsFlex(results, CT_MPC()));
    push!(resultsDict[:controller], concatResultsFlex(ctrlRes, CT_MPC()));
    println("Okay α=$α is done!")
    # Base.GC.gc();
end

## Save results using JSON3
folder = "../data/output/Run24/"
open(folder * "cs3_R24_MPC_BNoDeg_NMC_winter.json", "w") do f
    JSON3.pretty(f, resultsDict)
end

# Plotting functions
EMSData = build_data(;nEV = 1,
                    season = "winter",
                    profType = "monthly",
                    loadType = "GV",
                    year = 2023,
                    cellID = "SYNSANYO",
                    );
build_pbrom_sim_cell!(EMSData, s); # add simulation cell models
makeEMSplots(resultsDict[:results][1], EMSData; backend="GLMakie")
makeEMSplots(resultsDict[:controller][1], EMSData; backend="GLMakie", plot_ageing=false)