cd(@__DIR__)

# Import all the necessary packages
using JuMP, InfiniteOpt, KNITRO, EMSmodule, LiiBRA
using LinearAlgebra, Distributions, Statistics, Parameters, Random, Revise, Test
using DataFrames, LaTeXStrings, Printf, JSON3, Makie, CairoMakie, GLMakie
Random.seed!(1234);
includet("../functions/build_data.jl") # data structures

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

function handleInfeasible!(results::Dict, preSol::Dict, ts, s, typeOpt::DayAhead)
    # this function is used to handle infeasible solutions
    # it is called in the solvePolicies() function
    # If an MPC solution is infeasible implement the previous solution 
    # take the length from results because its the one handled by handleInfeasible()
    shift = Int(ceil(length(results["t"])/2))-1

    println("\r Checking the status of the last set of results")
    if results[:"status"] == LOCALLY_INFEASIBLE
    # discard the last set of results and implement the second step of the previous solution at ts-1
        println("\r Infeasible solution at step $ts, implementing the second step of the previous solution at step $(ts-1)")
        results=copy(preSol);
        for k in keys(results)
            # Check if the key requires special handling and skip the "status" key
            if k == :"status"
                results[k] = LOCALLY_INFEASIBLE;
            elseif k == :"γf"
                if s.nEV !== 1
                    results[k] = [results[k][n][(shift+1):end] for n ∈ 1:s.nEV];
                else
                    results[k] = results[k][(shift+1):end];
                end
            else
                length(results[k]) == 1 ? continue : results[k] = results[k][(shift+1):end];
            end
        end
        println(results["t"])
    else
        println("\r Locally feasible solution at step $ts")
    end
    # update the controller setpoints so that for next iteration this stored somewhere
    preSol = copy(results);
    return results, preSol;
end

function rollingHorizon(sDA, sCT, typeOpt::CT_MPC) # MPC or day-ahead
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

# define settings
resultsDict = Dict("A123"=>Dict("results"=>Dict{String, Any}(), "controller"=>Dict{String, Any}()),
     "SYNSANYO"=>Dict("results"=>Dict{String, Any}(), "controller"=>Dict{String, Any}()));
# Run the rolling horizon simulation
for cellID ∈ ["A123", "SYNSANYO"]
    Random.seed!(1234);
    # Define the weights
    Wgrid = 1; Wloss = 0.1;
    W = [Wgrid 1000 Wloss 1000];
    sDA=modelSettings(nEV=1, t0=0,Tw=48-1, Δt=1, steps=29, costWeights=W, season="summer",
                profType="monthly", loadType="GV", year=2023, cellID=cellID);
    sCT=modelSettings(nEV=1, t0=1/4,Tw=24-1/4, Δt=1/4, steps=96, costWeights=W, season="summer",
                profType="monthly",loadType="GV", year=2023, cellID=cellID);
    results, res_DA, ctrlRes, inputW_DA, inputW_CT, sDA, sCT = rollingHorizon(sDA, sCT, CT_MPC());
    # associate the results with the weights
    dict = Dict("results"=>concatResultsFlex(results, CT_MPC()),
                "controller"=>concatResultsFlex(ctrlRes, CT_MPC()));
    resultsDict[cellID] = dict;
    println("Okay $(cellID) is done!")
    # Base.GC.gc();
end

## Save results using JSON3
folder = "../data/output/Run24/"
open(folder * "cs2_R24_MPC_CPBDeg_W11e-1_LFPvsNMC_summer.json", "w") do f
    JSON3.pretty(f, resultsDict)
end

# Plotting functions
EMSData = build_data(;nEV = 1, season = "winter", profType = "monthly", loadType="GV", year = 2023, cellID = "A123");
makeEMSplots(resultsDict["A123"][1]["results"], EMSData; backend="GLMakie")
makeEMSplots(resultsDict["A123"][1]["controller"], EMSData; backend="GLMakie", plot_ageing=false)
makeEMSplots(resultsDict["SYNSANYO"][1]["results"], EMSData; backend="GLMakie")
## Save results using JSON3
open("../data/output/Run17/cs2_R17_DA_CPBDeg_W11_summer.json","w") do f
    JSON3.pretty(f, results)
end