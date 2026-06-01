cd(@__DIR__)

# Import all the necessary packages
using JuMP, InfiniteOpt, KNITRO
using LinearAlgebra, Distributions, Statistics, Parameters, Random, Revise, Test
using DataFrames, LaTeXStrings, Printf, JSON3, Makie, CairoMakie, GLMakie
Random.seed!(1234);
# nEV=1:2;
# includet("../functions/makeEMSobjs.jl") # data structures
includet("../functions/makeEMSobjs_perp.jl") # data structures
includet("../functions/build_data.jl") # data structures
includet("../functions/makeEMSplots.jl") # data structures
# includet("../functions/EMSfns.jl") # General functions
# includet("../functions/ESSfns.jl") # BESS functions
# includet("../functions/testEMS.jl") # Test functions
includet("../functions/EMSfns_perp.jl") # General functions
includet("../functions/ESSfns_perp.jl") # BESS functions
includet("../functions/testEMS_perp.jl") # Test functions
includet("../functions/simEMS.jl") # Plotting functions
# includet("../functions/makeForecasts.jl") # EMS functions

## Define all the necessary functions
function solvePolicies(optimizer, # model optimizer
    sets::modelSettings, # number of EVs, discrete time supports, etc.
    data::Dict, # information for the model (mainly parameters), these are the devices (EV, BESS, PV, etc.), the costs (interests, capex, etc) and the revenues
    preRes::Dict, # previous results to warm start the model,
    manager
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
                            "tuner"=>1,
                            "scale"=>1,
                            "presolve"=>1,
                            "outlev"=>3,
                            # "ms_enable"=>1,
                            # "ms_maxsolves"=>20,
                            "nlp_algorithm"=>3,
                            "numthreads"=>4,
                            "ms_numthreads"=>4,
                            "mip_numthreads"=>4,
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
    # Solve model
    optimize!(model)
    results = getResults_perp(model)
    # ddict, sdpdict = getDuals(model)
    # MOI.empty!(InfiniteOpt.backend(model))# free the license
    KNITRO.KN_release_license(manager)
    # return results, ddict, sdpdict;
    return results, manager
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
            elseif k == :"γ_cont"
                if s.nEV !== 1
                    results[k] = [results[k][n][(shift+1):end] for n ∈ 1:s.nEV];
                else
                    results[k] = results[k][(shift+1):end];
                end
            else
                length(results[k]) == 1 ? continue : results[k] = results[k][(shift+1):end];
            end
        end
    else
        println("\r Locally feasible solution at step $ts")
    end
    # update the controller setpoints so that for next iteration this stored somewhere
    preSol = copy(results);
    return results, preSol;
end

function rollingHorizon(s::modelSettings, termCond::Float64;
    typeOpt::String="MPC") # MPC or day-ahead
    ## Rolling Horizon Simulation.
    # This function simulates the EMS for a given number of steps.
    # The EMS is initialized at time t0 and then it is solved for a time window of Tw hours.
    # Then, the EMS is solved again for the next Tw hours, but this time the initial conditions are the ones obtained from the previous solution.
    # This process is repeated until the number of steps is reached.
    # The output of this function are:
    # - results::Vector{Dict}, where each dictionary contains the results of the EMS for each time window.
    # - rhDict::Dict, which is a dictionary with the concatenated results of the EMS for each time window.
    manager = KNITRO.LMcontext()
    steps = s.steps; Tw = s.Tw; Δt = s.Δt;
    Dt = s.dTime[1]:(s.dTime[2]-s.dTime[1]):s.dTime[end];
    # allocate memory
    results=Vector{Dict}(undef, steps);
    # ddict=Vector{Dict}(undef, steps); sdpdict=Vector{Dict}(undef, steps);
    controllerRes=Vector{Dict}(undef, steps);
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
    EMSData[:"grid"].λ /= 1000 # €/MWh --> €/kWh
    EMSData[:"grid"].λ /= 3600 # €/kWh --> €/kWs
    # Modify the TESS init condition
    EMSData[:"TESS"].SoC0=0.4;
    # # Modify the BESS perf. model
    EMSData[:"BESS"].GenInfo.SoC0=0.6;
    EMSData[:"BESS"].GenInfo.termCond=termCond;
    if s.cellID == "A123"
        EMSData[:"BESS"].GenInfo.SoCLim = [0.1, 0.95];
        [EMSData[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.1, 0.95] for n ∈ 1:s.nEV];
        # EMSData[:"BESS"].GenInfo.PowerLim = [-10.0, 10.0];
    elseif s.cellID == "SYNSANYO"
        EMSData[:"BESS"].GenInfo.SoCLim = [0.15, 1.0];
        [EMSData[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.15, 1.0] for n ∈ 1:s.nEV];
    end
    EMSData[:"BESS"].PerfParameters = bucketPerfParams();
    # Modify the EV model
    [EMSData[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams() for n ∈ 1:s.nEV];
    # # Modify the BESS ageing model
    # EMSData[:"BESS"].AgingParameters = empAgingParams();
    # # Modify the EV model
    # [EMSData[:"EV"][n].carBatteryPack.AgingParameters = empAgingParams() for n ∈ 1:s.nEV];

    # Initial solution
    results[1]["compTime"] = @elapsed begin
        results[1], manager = solvePolicies(() -> KNITRO.Optimizer(; license_manager = manager), s, EMSData, Dict(), manager);
        # controllerRes[1] = copy(results[1]);
        # Check the status of the last set of results
        results[1], controllerRes[1] = handleInfeasible!(results[1], Dict(), 1, s, DayAhead());
        simTransitionFun!(results[1], EMSData, s, DayAhead())
        # move time window
        # typeOpt == "MPC" ? Dt = Dt .+ Δt*3600.0 : Dt = Dt.+ 24*3600.0;
        Dt = Dt.+ 24*3600.0;
        # update_forecasts(EMSData, Dt)
        s.dTime=collect(Dt);
    end

    for ts in 2:steps
        try
            results[ts]["compTime"] = @elapsed begin
                # build+solve models
                results[ts], manager =solvePolicies(() -> KNITRO.Optimizer(; license_manager = manager), s, EMSData, Dict(), manager);
                # Check the status of the last set of results
                results[ts], controllerRes[ts] = handleInfeasible!(results[ts], controllerRes[ts-1], ts, s, DayAhead());
                simTransitionFun!(results[ts], EMSData, s, DayAhead())
                # move time window
                Dt = Dt.+ 24*3600.0;
                # update_forecasts(EMSData, Dt)
                s.dTime=collect(Dt);
            end
            println("\r Step $ts out of $steps done!")
        catch e
            println("Error at step $ts")
            showerror(stdout, e, catch_backtrace())
            return results, EMSData,s, controllerRes;
        end
    end
    return results, EMSData,s, controllerRes;
end

function rollingHorizon(s,
    seed :: Vector{Dict}; # seed for the warm start
    typeOpt::String="MPC", # MPC or day-ahead
    )
    ## Rolling Horizon Simulation.
    # This function simulates the EMS for a given number of steps.
    # The EMS is initialized at time t0 and then it is solved for a time window of Tw hours.
    # Then, the EMS is solved again for the next Tw hours, but this time the initial conditions are the ones obtained from the previous solution.
    # This process is repeated until the number of steps is reached.
    # The output of this function are:
    # - results::Vector{Dict}, where each dictionary contains the results of the EMS for each time window.
    # - rhDict::Dict, which is a dictionary with the concatenated results of the EMS for each time window.
    @assert typeOpt ∈ ["MPC", "day-ahead"] "typeOpt must be either MPC or DA"
    @assert length(seed) == s.steps "The length of the seed must be equal to the number of steps"
    if typeOpt == "day-ahead"
        @assert (s.Tw + s.Δt) % 24 == 0 "Tw must be a multiple of fs times 24"
    end

    steps = s.steps; Tw = s.Tw; Δt = s.Δt;
    Dt = s.dTime[1]:(s.dTime[2]-s.dTime[1]):s.dTime[end];
    # allocate memory
    results=Vector{Dict}(undef, steps);
    ddict=Vector{Dict}(undef, steps);
    sdpdict=Vector{Dict}(undef, steps);
    controllerRes=Vector{Dict}(undef, steps);
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
    EMSData[:"grid"].λ /= 1000 # €/MWh --> €/kWh
    EMSData[:"grid"].λ /= 3600 # €/kWh --> €/kWs
    # Modify the TESS init condition
    EMSData[:"TESS"].SoC0=0.4;
    # # Modify the BESS perf. model
    EMSData[:"BESS"].GenInfo.SoC0=0.6;
    if s.cellID == "A123"
        EMSData[:"BESS"].GenInfo.SoCLim = [0.1, 0.95];
        [EMSData[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.1, 0.95] for n ∈ 1:s.nEV];
        # EMSData[:"BESS"].GenInfo.PowerLim = [-10.0, 10.0];
    elseif s.cellID == "SYNSANYO"
        EMSData[:"BESS"].GenInfo.SoCLim = [0.15, 1.0];
        [EMSData[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.15, 1.0] for n ∈ 1:s.nEV];
    end
    # EMSData[:"BESS"].PerfParameters = bucketPerfParams();
    # # Modify the EV model
    # [EMSData[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams() for n ∈ 1:s.nEV];
    # Modify the BESS ageing model
    # EMSData[:"BESS"].AgingParameters = empAgingParams();
    # # Modify the EV model
    # [EMSData[:"EV"][n].carBatteryPack.AgingParameters = empAgingParams() for n ∈ 1:s.nEV];

    # Initial solution
    times[1] = @elapsed begin
        results[1], ddict[1], sdpdict[1] = solvePolicies(KNITRO.Optimizer, s, EMSData, seed[1]);
        controllerRes[1] = copy(results[1]);
        # Check the status of the last set of results
        handleInfeasible!(results, controllerRes, 1, s; typeOpt=typeOpt);
        simTransitionFun!(results[1], EMSData, s; typeOpt = typeOpt)
        # move time window
        typeOpt == "MPC" ? Dt = Dt .+ Δt*3600.0 : Dt = Dt.+ 24*3600.0;
        # update_forecasts(EMSData, Dt)
        s.dTime=collect(Dt);
    end

    for ts in 2:steps
        try
            times[ts] = @elapsed begin
                # build+solve model
                # results[ts]=solvePolicies(KNITRO.Optimizer, s, EMSData, controllerRes[ts-1]);
                results[ts], ddict[ts], sdpdict[ts] =solvePolicies(KNITRO.Optimizer, s, EMSData, seed[ts]);
                controllerRes[ts] = copy(results[ts]);
                # Check the status of the last set of results
                handleInfeasible!(results, controllerRes, ts, s; typeOpt=typeOpt);
                # EMSData=update_measurements(results[ts], s, EMSData; typeOpt=typeOpt);
                simTransitionFun!(results[ts], EMSData, s; typeOpt = typeOpt)
                # move time window
                # typeOpt == "MPC" ? Dt = Dt .+ Δt*3600.0 : Dt = Dt.+ (Tw+Δt)/d*3600.0;
                typeOpt == "MPC" ? Dt = Dt .+ Δt*3600.0 : Dt = Dt.+ 24*3600.0;
                # update_forecasts(EMSData, Dt)
                s.dTime=collect(Dt);
            end
            # Base.GC.gc(); # collect garbage
            println("\r Step $ts out of $steps done!")
        catch e
            println("Error at step $ts")
            showerror(stdout, e, catch_backtrace())
            return results, EMSData,s, controllerRes;
        end
    end
    [results[ts]["compTime"]=times[ts] for ts in 1:steps]
    return results, EMSData,s, controllerRes;
end

# Define the weights
# opt weights
Wgrid = 1;
Wloss = 0;
W=[Wgrid 1000 Wloss 1000];
# Run the rolling horizon simulation
terminal_conditions = [0., 4., 6., 12., 16.]
rlts = Vector{NamedTuple}(undef, length(terminal_conditions))
costs = Vector{NamedTuple}(undef, length(terminal_conditions))
for (i, tc) ∈ enumerate(terminal_conditions)
    s=modelSettings(nEV=1, t0=1/4,Tw=48-1/4, Δt=1/4, steps=5, costWeights=W, season="summer",
                profType="weekly", loadType="GV", year=2023, cellID="SYNSANYO");
    results,EMSData,s, controllerRes = rollingHorizon(s, tc; typeOpt="day-ahead")
    rhDict=concatResultsRH(results, DayAhead()); # concatenate results
    objective = calcFullObj(rhDict, EMSData, s)
    println("terminal condition $(tc) is ok.")
    rlts[i] = (termCond = tc, rlt = rhDict)
    costs[i] = (termCond = tc, obj = objective)
end
function plotTerminalCond(rlts::Vector{NamedTuple}, costs::Vector{NamedTuple}; filename="cs1_BNoDeg_termCond.pdf")
    CairoMakie.activate!()
    set_theme!(theme_latexfonts(), fontsize = 20)
    f = Figure(size=(900, 400))
    ax_soc = Axis(f[1,1:2],
        xlabel = "t [hr]",
        ylabel = L"SoC_{\text{BESS},t}");
    ax_cost = Axis(f[1,3],
        xlabel = "t [hr]",
        ylabel = L"J");
    terminal_conditions = [rlts[i].termCond for i ∈ 1:length(rlts)]
    for (i, tc) ∈ enumerate(terminal_conditions)
        t = rlts[i].rlt["t"]/3600
        SoCbess = rlts[i].rlt["SoCbess"]
        stairs!(ax_soc, t, SoCbess, label = "$tc")
    end
    axislegend(ax_soc, title = L"Periodic Condition $t_1$", position = :lt)
    objValues = [cost.obj.objVal for cost ∈ costs]
    barplot!(ax_cost, terminal_conditions, objValues)
    save(filename,f)
    return f
end
folder = "../images/AppEnergy/cs1/"
plotTerminalCond(rlts, costs; filename = folder*"cs1_BNoDeg_summer_termCond.pdf")
s=modelSettings(nEV=1, t0=1/4,Tw=48-1/4, Δt=1/4, steps=5, costWeights=W, season="summer",
                profType="weekly", loadType="GV", year=2023, cellID="SYNSANYO");
folder = "../data/output/Run23/"
open(folder * "cs1_R23_rlts_DA_BNoDeg_W10_$(s.cellID)_y$(s.year)$(s.season)d$(s.steps)_termCond.json", "w") do f
    JSON3.pretty(f, rlts)
end
open(folder * "cs1_R23_costs_DA_BNoDeg_W10_$(s.cellID)_y$(s.year)$(s.season)d$(s.steps)_termCond.json", "w") do f
    JSON3.pretty(f, costs)
end