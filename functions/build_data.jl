################ EXOGENOUS INFO BUILDING ################
using Parameters, CSV, Dates, Serialization, Interpolations

function processPrices(df;
    type::String="raw",
    upSampRatio::Int=4, # samples per hour
    profType::String="daily",
    season::String="summer",
    market::String = "day-ahead",
    )
    # This function processes the prices data from EPEX in several formats.
    # type:
    # - "raw" raw data from EPEX
    # - "summary". Summary data coming from "makePriceForecast.ipynb" dashboard.
    # The price data is in [€/MWh]
    @assert type ∈ ["raw", "summary"] "Invalid type of data"
    @assert profType ∈ ["daily", "weekly", "biweekly","monthly","yearly"] "Invalid profile type"
    profType == "yearly" ? nothing : @assert season ∈ ["summer", "winter"] "Invalid season";
    
    if type=="raw"
        # Raw EPEX FTP server data
        # eliminate clearence hour B column
        if market == "day-ahead"
            select!(df, Not(:"Hour 3B"));
        else # intraday
            select!(df, Not([:"Hour 3B Q1", :"Hour 3B Q2", :"Hour 3B Q3", :"Hour 3B Q4"]));
        end
        df=coalesce.(df, 0) # replace missing values with 0
        
        # Now we need to reshape the DataFrame into a timeseries
        priceData = Vector{Float64}();
        for r in 1:nrow(df)
            row=Vector{Float64}(df[r,2:25])
            append!(priceData, row)
        end
    elseif type=="summary"
        priceData=Vector{Float64}(df[:,:mean])
    end
    # Change resolution of the prices. From 1h/sample to 15min/sample
    priceData = repeat(priceData, inner=upSampRatio);
    priceData=[0; priceData]
    if profType == "yearly"
        return priceData = [priceData priceData.*0.95]
    end
    # for the rest of the profiles you continue with the seasonal profiles
    priceData=getSeasonalProfiles(priceData; type=profType, n_samples_per_hour = upSampRatio)[season];
    # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
    if profType == "biweekly"
        priceData=repeat(priceData[1:(end-upSampRatio*24)], outer=2); append!(priceData, priceData[1:upSampRatio*24])
    end
    # Buy < sell prices 
    priceData = [priceData priceData.*0.95]
    # priceData=CSV.read("energy_prices.csv", DataFrame);  # old data from Wil
    return priceData;
end;

function build_data(; nEV::Int64=2, # number of EVs in the system
    season::String = "winter", # "summer" or "winter"
    profType::String, # "weekly" or "daily"
    loadType::String, # "GV" or "mffbas" or "base_models"
    year::Int64,
    fs::Int64 = 4, # samples per hour
    cellID::String = "SYNSANYO" # cell ID for the battery packs
    )
    # This function builds the data for the optimization problem.
    # The inputs are:
    # - season: "summer" or "winter"
    # - profType: "weekly" or "daily"
    # - loadType: "GV" or "mffbas" or "base_models"
    # - year: the year of the simulation
    # The output is a dictionary that contains the models for all the different devices in the Multi-carrier Energy System.
    # Check inputs
    # @assert typeof(nEV) == UnitRange{Int64} "nEV must be a range of integers"
    @assert typeof(nEV) == Int64 "max nEV must be a integer"
    @assert season ∈ ["summer", "winter"] "Invalid season"
    @assert profType ∈ ["daily", "weekly", "biweekly", "monthly", "yearly"] "Invalid profile type"
    @assert loadType ∈ ["GV", "mffbas", "base_models"] "Invalid load type"
  
    ## PV model
    if loadType == "base_models" # From Joel's Base models
        MPPT = CSV.read("../Base models/PV_15min.csv", DataFrame, header=false) # measurement of the max. power point tracking
    else
        pvData = CSV.read("../data/input/SPVMPPTData.csv", DataFrame, delim=',', header = 1)
        MPPT = pvData.MPPT; # measurement of the max. power point tracking
    end
    MPPT=Array(MPPT)[1:Int(4/fs):end];
    Npv=1; # number of pv panels
    if profType ≠ "yearly"
        # get the daily seasonal profile for the PV
        # MPPT=getSeasonalProfiles(MPPT, type = profType)[season];
        MPPT=getSeasonalProfiles(MPPT; n_samples_per_hour = fs, type = profType)[season];
        # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
        if profType == "biweekly" MPPT=repeat(MPPT[1:(end-fs*24)], outer=2); end
    end
    append!(MPPT, MPPT[1:fs*24*2]); # append the first day to the end
    # MPPT=MPPT[1:end-1];
    spvModel = SPVData(MPPTData = MPPT*Npv);
    
    ## BESS model
    # A battery pack model is composed of three parts:
    # - its general information, contained in the GenInfo type.
    # - its performance submodel, contained in the PerfParameters type. This describes the SoC and terminal voltage of a cell.
    # - its ageing submodel, contained in the AgingParameters type. This describes the evolution of the performance sub-model parameters.
    
    # Read data from the E2 of the ESCtoolbox from Plett.
        # OCV non linear data model
        bessOCV=JSON3.read(open("../data/input/cell_models/$cellID-modelocv.json", "r"),
                    Dict{String, Vector{Float64}});
        # Dynamic ECM model, with np=1
        bessDYN=JSON3.read(open("../data/input/cell_models/$cellID-modeldyn-no-hys.json", "r"),
                    Dict{String, Union{Vector{Float64}, Vector}});
        # bessOCVmat = matread("../data/input/cell_models/E2model-ocv.mat")
        # # Dynamic ECM model
        # bessDYNmat = matread("../data/input/cell_models/E2model.mat")
        # # Unpack
        # bessOCVmat=bessOCVmat["model"];
        # bessDYNmat=bessDYNmat["model"];

        # Bucket model info
        Q0 = bessOCV["OCVQ"]; Q0 = mean(Q0); # [Ah/cell]
        η = bessOCV["OCVeta"]; η = mean(η); # [p.u.]

        # ECM info
        T=bessDYN["temps"]; T=T[1,:]; # [C]
        R0Param=bessDYN["R0Param"]; # [Ohm]
        RParam=bessDYN["RParam"]; # [Ohm]
        RCParam=bessDYN["RCParam"]; # τ [s]
        Tind25=findall(x->x==25, T)
    
    # Heliox 43kWh battery pack seems a little bit too much lets use a 20kWh pack.
    PowerLim = [-12.5, 12.5]; P0 = 0; # Power limits and initial condition [kW]. 
    SoCLim = [0.2, 1.]; SoC0 = 0.4; # SoC limits and initial condition [p.u.].
    termCond = 6.; # Termination condition [hr].
    ηC = 0.95; # Charger efficiency [p.u.].
    # State of Health
    SoHQ = 1; SoHR0 = R0Param[Tind25][1];
    if cellID == "SYNSANYO"
        # Np = 10; Ns = 100; # Branches in parallel and series cells per branch.
        # lets try a smaller BESS
        Np = 7; Ns = 100; # Branches in parallel and series cells per branch.
        ocv_params=OCVlinearPerfParams()
        vLim = [2.8, 4.2]
        # Ageing submodel
        # aging_params=empAgingParams();
        aging_params=JinAgingParams();
        Cell = Construct("LG M50")
        SList = [1., 0.15]
    elseif cellID == "A123"
        Np = 25; Ns = 110; # Branches in parallel and series cells per branch.
        # first option, 
        # ocv_params.ocvLine = [2.5, 1.1]
        # second option from 20% to 95% SoC
        ocv_params=OCVlinearPerfParams(ocvLine = [3.2, 0.2105])
        vLim = [2.0, 3.6]
        # Ageing submodel
        # aging_params=empAgingParams();
        Cell = Construct("A123");
        SList = [0.95, 0.1]
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
    gen_params = Generic(PowerLim, P0, SoCLim, SoC0, termCond, Q0, SoHQ, SoHR0, Np, Ns, η, ocv_params, vLim, ηC, initVal);
    # Performance submodel
    perf_params = ECMPerfParams(R0Param=R0Param[Tind25], RParam=RParam[Tind25], RCParam=RCParam[Tind25])
    # Physics-based performance submodel
    Cell.RA.H1 = Cell.RA.H2 = [1:2500; 3000:3500; 4000:4500]
    Base.invokelatest(Spatial!, Cell, 4, 2)
    A, B, C, D = Base.invokelatest(Realise, Cell, SList);
    perf_PB_params = ROMPerfParams(SList=SList, Cell=Cell, A=A, B=B, C=C, D=D);
    perf_PB_params.Cell.Neg.RFilm = R0Param[Tind25][1];
    # wrap everything in a BESSData type
    battModel = BESSData(gen_params, perf_params, aging_params, cellID);
    battPBModel = BESSData(gen_params, perf_PB_params, aging_params, cellID);

    ## EV MODEL
    # An electric vehicle model is composed of two parts:
    # - its battery pack, contained in a BESSData type.
    # - its driving submodel, contained in the driveData type. This describes the availability, 
    # times of departure and arrival, and the reference SoC of the EV.
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
    gen_params = [Generic(PowerLim, P0[n], SoCLim, SoC0[n], termCond, Q0[n], SoHQ, SoHR0, Np, Ns, η, ocv_params, vLim, ηC, initVal) for n ∈ 1:nEV]
    perf_params=[ECMPerfParams(R0Param=R0Param[Tind25], RParam=RParam[Tind25],
                RCParam=RCParam[Tind25], vt0=vt0[n]) for n ∈ 1:nEV]
    batteryPack = [BESSData(gen_params[n], perf_params[n], aging_params, cellID) for n ∈ 1:nEV]
    PB_batteryPack = [BESSData(gen_params[n], perf_PB_params, aging_params, cellID) for n ∈ 1:nEV]

    # Driving information definition
    μD = 0.5; σD=1.; # Parameters for the Gaussian distributions
    refSoC=[0.85, 0.85]; # user requirement
    # availability
    avObj = [availabilityEV(length(MPPT), fs,  μD, σD, Ns, Np, vLim[2], Q0[n], PowerLim[2] * ηC) for n ∈ 1:nEV];
    γ = [avObj[n][1] for n ∈ 1:nEV];
    tDep = [avObj[n][2] for n ∈ 1:nEV];
    tArr = [avObj[n][3] for n ∈ 1:nEV];
    Pdrive = [avObj[n][4] for n ∈ 1:nEV];
    drive_info = [driveData(Pdrive[n], refSoC[n], γ[n], tDep[n], tArr[n]) for n in 1:nEV]

    # wrap everything in a EVData type
    evModel = [EVData(batteryPack[n], drive_info[n]) for n in 1:nEV]
    evPBModel = [EVData(PB_batteryPack[n], drive_info[n]) for n in 1:nEV]

    # Solar thermal model
    Pn=0.5; capex= 1500; ηST = 0.6;
    stModel = ElectroThermData(Pn, capex, ηST);

    # TESS model
    tessModel = TESSData();
    # Power Electronic Interface model
    peiModel = peiData();

    ## Grid Model
    # The grid is represented by:
    # - λ energy prices. [buy; sell] 
    # - loadE electrical load measurement
    # - loadTh thermal load measurement
    # for anual profiles the length is (365*24*fs)+1=35041.

    ## pick file paths, read and convert to array
    # for the prices
    # Raw EPEX FTP server data
    pricePath = "../data/input/EPEX/auction_spot_prices_netherlands_$year.csv"
    spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=2);

    #= for IECON 2023 =#
        # summary data
        # pricePath = "../data/input/EPEX/summaryPrice_Winter.csv";
        # spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=1);
        # pricePath = "../data/input/EPEX/summaryPrice_Summer.csv";
        # spotPriceDF = CSV.read(pricePath, DataFrame, delim=',', header=1);

    # for the electrical load
    if loadType == "GV"
        # processed synthezided load profile 1 year
        loadEPath = "../data/input/GV/Load_1.csv";
        loadE = CSV.read(loadEPath, DataFrame, header = false) # electric load
        # loadE = CSV.read("../data/input/GV/Load_1.csv", DataFrame, header=false) # electric load
        loadE= Vector(loadE[!,1]);
    elseif loadType == "mffbas"
        # From Market Facilitation Forum (MFF) and the Beheerder Afspraken Stelsel (BAS) i.e. mffbas
        loadEPath = "../data/input/mffbas/summaryE1_$year.csv"
        #= for IECON 2023 =#
        # loadEPath = "../data/input/mffbas/summaryE1_mean_Winter.csv"
        # loadEPath = "../data/input/mffbas/summaryE1_mean_Summer.csv"
        loadE = CSV.read(loadEPath, DataFrame) # electric load
        loadE= Vector(loadE[!,2]); 
        
    elseif loadType == "base_models" # Joel's base models
        loadEPath = "../data/input/Base models/Load_Profile_15min.csv";
        loadE = CSV.read(loadEPath, DataFrame, transpose=true, header=false)
        loadE= Vector(loadE[!,1]); 
    end
    # normalize loadE
    loadE=loadE./maximum(loadE);
    # peak of 2kW
    loadE=loadE*2.5;
    # downsample
    loadE = loadE[1:Int(4/fs):end]
    # for the thermal load
    loadThPath= "../Base models/Thermal_load_15min.csv";    
    loadTh = CSV.read(loadThPath, DataFrame, header = false) # thermal load
    loadTh= Vector(loadTh[!,1]); loadTh=loadTh[1:Int(4/fs):end]*1e-3;
    
    # check length of loadE, loadTh
    length(loadE) == (365*24*fs+1) ? nothing : loadE = [loadE[1]; loadE];

    # Data processing (upsampling, seasonal patterns, etc.)
    priceData = processPrices(spotPriceDF; upSampRatio = fs, type="raw", profType = profType, season = season);
    if profType ≠ "yearly"
        # priceData = processPrices(spotPriceDF; type="raw", profType = profType, season = season);
        # get the profile for the loads
        loadE=getSeasonalProfiles(loadE; n_samples_per_hour = fs, type = profType)[season];
        loadTh=getSeasonalProfiles(loadTh; n_samples_per_hour = fs, type = profType)[season];
        # for biweekly profiles we need to repeat the weekly profile twice and append the first day to the end
        if profType == "biweekly"
            loadE=repeat(loadE[1:(end-fs*24)], outer=2);
            loadTh=repeat(loadTh[1:(end-fs*24)], outer=2);
        end
    else
        # priceData = processPrices(spotPriceDF; type="raw", profType = profType);
        priceData = priceData[1:end-1,:];
        loadE = loadE[1:end-1];
        loadTh = loadTh[1:end-1];
    end
    priceData = vcat(priceData, priceData[1:fs*24*2,:])
    append!(loadE, loadE[1:fs*24*2])
    append!(loadTh, loadTh[1:fs*24*2])
    # Grid connection limits
    gridModel = gridData([-17, 17], 0.9, priceData, loadE, loadTh);
    
    # Heat pump model
    # only a uniderectional (heating) HP for now
    Pn=4; capex= 500; ηHP = 3; 
    # The initial condition follow the power balance
    hpModel = ElectroThermData(Pn, capex, ηHP);

    data=Dict("SPV"=>spvModel, 
        "BESS"=>battModel, "BESSpbrom"=>battPBModel,
        "EV"=>evModel, "EVpbrom"=>evPBModel,
        "ST"=>stModel, "HP"=>hpModel, "TESS"=>tessModel,
        "grid"=>gridModel, "PEI"=>peiModel);
    return data
end

function start_data!(data::Dict, s::modelSettings)
    data[:"grid"].λ /= (1000*3600) # €/MWh --> €/kWs
    # Modify the TESS init condition
    data[:"TESS"].SoC0=0.4;
    # # Modify the BESS perf. model
    data[:"BESS"].GenInfo.SoC0=0.6;
    if s.cellID == "A123"
        data[:"BESS"].GenInfo.SoCLim = [0.1, 0.95];
        [data[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.1, 0.95] for n ∈ 1:s.nEV];
        # data[:"BESS"].GenInfo.PowerLim = [-10.0, 10.0];
    elseif s.cellID == "SYNSANYO"
        data[:"BESS"].GenInfo.SoCLim = [0.15, 1.0];
        [data[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.15, 1.0] for n ∈ 1:s.nEV];
    end
    data[:"BESS"].PerfParameters = bucketPerfParams();
    # Modify the EV model
    [data[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams() for n ∈ 1:s.nEV];
    # # Modify the BESS ageing model
    # data[:"BESS"].AgingParameters = empAgingParams();
    # # Modify the EV model
    # [data[:"EV"][n].carBatteryPack.AgingParameters = empAgingParams() for n ∈ 1:s.nEV];
    return data
end

function start_data!(inputW_CT::Dict, inputW_DA::Dict, s)
    # Exogen Info
    for n ∈ 1:s.nEV
        inputW_CT["EV"][n].driveInfo.Pdrive = copy(inputW_DA["EV"][n].driveInfo.Pdrive)
        inputW_CT["EV"][n].driveInfo.tDep = copy(inputW_DA["EV"][n].driveInfo.tDep)
        inputW_CT["EV"][n].driveInfo.tArr = copy(inputW_DA["EV"][n].driveInfo.tArr)
        inputW_CT["EV"][n].driveInfo.γ = repeat(inputW_DA["EV"][n].driveInfo.γ,inner=4)
    end
    # Convert the cost of electricity from €/MWh to €/kWs
    inputW_DA["grid"].λ /= (1000*3600)
    inputW_CT["grid"].λ = hcat([repeat(inputW_DA["grid"].λ[:,r], inner = 4) for r ∈ 1:2])
    # Set initial conditions
    # Modify the TESS init condition
    inputW_DA[:"TESS"].SoC0=0.4; inputW_CT[:"TESS"].SoC0=0.4;
    # Modify the BESS perf. model
    inputW_DA[:"BESS"].GenInfo.SoC0=0.6; inputW_CT[:"BESS"].GenInfo.SoC0=0.6;
    inputW_CT["BESS"].GenInfo.termCond = 0.
    if s.cellID == "A123"
        inputW_CT[:"BESS"].GenInfo.SoCLim = [0.1, 0.95];
        inputW_DA[:"BESS"].GenInfo.SoCLim = [0.1, 0.95];
        for n ∈ 1:s.nEV
            inputW_CT[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.1, 0.95]
            inputW_DA[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.1, 0.95];
        end
        # data[:"BESS"].GenInfo.PowerLim = [-10.0, 10.0];
    elseif s.cellID == "SYNSANYO"
        inputW_CT[:"BESS"].GenInfo.SoCLim = [0.15, 1.0];
        inputW_DA[:"BESS"].GenInfo.SoCLim = [0.15, 1.0];
        for n ∈ 1:s.nEV
            inputW_DA[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.15, 1.0]
            inputW_CT[:"EV"][n].carBatteryPack.GenInfo.SoCLim = [0.15, 1.0]
        end
    end
    # inputW_CT[:"BESS"].PerfParameters = bucketPerfParams();
    # inputW_DA[:"BESS"].PerfParameters = bucketPerfParams();

    # # Modify the EV model
    # for n ∈ 1:s.nEV
    #     inputW_CT[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams();
    #     inputW_DA[:"EV"][n].carBatteryPack.PerfParameters = bucketPerfParams();
    # end
    # Modify the BESS ageing model
    # inputW_DA[:"BESS"].AgingParameters = empAgingParams();
    # inputW_CT[:"BESS"].AgingParameters = empAgingParams();
    # # Modify the EV model
    # [inputW_DA[:"EV"][n].carBatteryPack.AgingParameters = empAgingParams() for n ∈ 1:s.nEV];
    # [inputW_CT[:"EV"][n].carBatteryPack.AgingParameters = empAgingParams() for n ∈ 1:s.nEV];
    return inputW_CT, inputW_DA
end

function sync_S0!(inputW_CT::Dict, inputW_DA::Dict, sDA::modelSettings)
    inputW_CT["TESS"].SoC0 = copy(inputW_DA["TESS"].SoC0)
    inputW_CT["BESS"].GenInfo.SoC0 = copy(inputW_DA["BESS"].GenInfo.SoC0)
    inputW_CT["BESS"].GenInfo.SoHQ = copy(inputW_DA["BESS"].GenInfo.SoHQ)
    inputW_CT["BESS"].GenInfo.SoHR0 = copy(inputW_DA["BESS"].GenInfo.SoHR0)
    inputW_CT[:"BESS"].AgingParameters.z100p = copy(inputW_DA[:"BESS"].AgingParameters.z100p)
    inputW_CT[:"BESS"].AgingParameters.δSEI0 = copy(inputW_DA[:"BESS"].AgingParameters.δSEI0)
    inputW_CT[:"BESS"].AgingParameters.εₑ0 = copy(inputW_DA[:"BESS"].AgingParameters.εₑ0)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm)

    for n ∈ 1:sDA.nEV
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoC0)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHQ)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHR0)
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.z100p = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.z100p);
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.δSEI0 = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.δSEI0);
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.εₑ0 = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.εₑ0);
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100)
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm)
    end
    return inputW_CT, inputW_DA    
end

function sync_S0!(inputW_CT::Dict, inputW_DA::Dict, sDA::modelSettings, bessAging::JinAgingParams)
    inputW_CT["TESS"].SoC0 = copy(inputW_DA["TESS"].SoC0)
    inputW_CT["BESS"].GenInfo.SoC0 = copy(inputW_DA["BESS"].GenInfo.SoC0)
    inputW_CT["BESS"].GenInfo.SoHQ = copy(inputW_DA["BESS"].GenInfo.SoHQ)
    inputW_CT["BESS"].GenInfo.SoHR0 = copy(inputW_DA["BESS"].GenInfo.SoHR0)
    inputW_CT[:"BESS"].AgingParameters.z100p = copy(inputW_DA[:"BESS"].AgingParameters.z100p)
    inputW_CT[:"BESS"].AgingParameters.δSEI0 = copy(inputW_DA[:"BESS"].AgingParameters.δSEI0)
    inputW_CT[:"BESS"].AgingParameters.εₑ0 = copy(inputW_DA[:"BESS"].AgingParameters.εₑ0)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm)

    for n ∈ 1:sDA.nEV
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoC0)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHQ)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHR0)
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.z100p = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.z100p);
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.δSEI0 = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.δSEI0);
        inputW_CT[:"EV"][n].carBatteryPack.AgingParameters.εₑ0 = copy(inputW_DA[:"EV"][n].carBatteryPack.AgingParameters.εₑ0);
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100)
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm)
    end
    return inputW_CT, inputW_DA    
end
function sync_S0!(inputW_CT::Dict, inputW_DA::Dict, sDA::modelSettings, bessAging::empAgingParams)
    inputW_CT["TESS"].SoC0 = copy(inputW_DA["TESS"].SoC0)
    inputW_CT["BESS"].GenInfo.SoC0 = copy(inputW_DA["BESS"].GenInfo.SoC0)
    inputW_CT["BESS"].GenInfo.SoHQ = copy(inputW_DA["BESS"].GenInfo.SoHQ)
    inputW_CT["BESS"].GenInfo.SoHR0 = copy(inputW_DA["BESS"].GenInfo.SoHR0)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.θ_100)
    inputW_CT[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"BESSpbrom"].PerfParameters.Cell.Neg.RFilm)

    for n ∈ 1:sDA.nEV
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoC0)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHQ)
        inputW_CT["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(inputW_DA["EV"][n].carBatteryPack.GenInfo.SoHR0)
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100 = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.θ_100)
        inputW_CT[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm = copy(inputW_DA[:"EVpbrom"][n].carBatteryPack.PerfParameters.Cell.Neg.RFilm)
    end
    return inputW_CT, inputW_DA    
end

function sync_ess!(data::Dict)
    nEV = length(data["EVpbrom"]);
    data["BESS"].GenInfo.SoC0 = copy(data["BESSpbrom"].GenInfo.SoC0);
    data["BESS"].GenInfo.SoHQ = copy(data["BESSpbrom"].GenInfo.SoHQ);
    data["BESS"].GenInfo.SoHR0 = copy(data["BESSpbrom"].GenInfo.SoHR0);
    # data[:"BESS"].AgingParameters.z100p = copy(data["BESSpbrom"].AgingParameters.z100p);
    # data[:"BESS"].AgingParameters.δSEI0 = copy(data["BESSpbrom"].AgingParameters.δSEI0);
    # data[:"BESS"].AgingParameters.εₑ0 = copy(data["BESSpbrom"].AgingParameters.εₑ0);
    for n ∈ 1:nEV
        data["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoC0)
        data["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHQ)
        data["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHR0)
        # data[:"EV"][n].carBatteryPack.AgingParameters.z100p = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.z100p);
        # data[:"EV"][n].carBatteryPack.AgingParameters.δSEI0 = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.δSEI0);
        # data[:"EV"][n].carBatteryPack.AgingParameters.εₑ0 = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.εₑ0);
    end    
    return data
end

function sync_ess!(data::Dict, bessAging::JinAgingParams)
    nEV = length(data["EVpbrom"]);
    data["BESS"].GenInfo.SoC0 = copy(data["BESSpbrom"].GenInfo.SoC0);
    data["BESS"].GenInfo.SoHQ = copy(data["BESSpbrom"].GenInfo.SoHQ);
    data["BESS"].GenInfo.SoHR0 = copy(data["BESSpbrom"].GenInfo.SoHR0);
    data[:"BESS"].AgingParameters.z100p = copy(data["BESSpbrom"].AgingParameters.z100p);
    data[:"BESS"].AgingParameters.δSEI0 = copy(data["BESSpbrom"].AgingParameters.δSEI0);
    data[:"BESS"].AgingParameters.εₑ0 = copy(data["BESSpbrom"].AgingParameters.εₑ0);
    for n ∈ 1:nEV
        data["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoC0)
        data["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHQ)
        data["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHR0)
        data[:"EV"][n].carBatteryPack.AgingParameters.z100p = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.z100p);
        data[:"EV"][n].carBatteryPack.AgingParameters.δSEI0 = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.δSEI0);
        data[:"EV"][n].carBatteryPack.AgingParameters.εₑ0 = copy(data["EVpbrom"][n].carBatteryPack.AgingParameters.εₑ0);
    end    
    return data
end
function sync_ess!(data::Dict, bessAging::empAgingParams)
    nEV = length(data["EVpbrom"]);
    data["BESS"].GenInfo.SoC0 = copy(data["BESSpbrom"].GenInfo.SoC0);
    data["BESS"].GenInfo.SoHQ = copy(data["BESSpbrom"].GenInfo.SoHQ);
    data["BESS"].GenInfo.SoHR0 = copy(data["BESSpbrom"].GenInfo.SoHR0);    
    for n ∈ 1:nEV
        data["EV"][n].carBatteryPack.GenInfo.SoC0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoC0)
        data["EV"][n].carBatteryPack.GenInfo.SoHQ = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHQ)
        data["EV"][n].carBatteryPack.GenInfo.SoHR0 = copy(data["EVpbrom"][n].carBatteryPack.GenInfo.SoHR0)
    end    
    return data
end