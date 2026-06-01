# plotFigs_AppEnergyDA.jl
cd(@__DIR__)
using JuMP, InfiniteOpt, KNITRO, EMSmodule, LiiBRA
using LinearAlgebra, Distributions, Statistics, Parameters, Random, Revise, Test
using DataFrames, LaTeXStrings, Printf, JSON3, Makie, CairoMakie, GLMakie, ColorSchemes
Random.seed!(1234);
# nEV=1:2;
# includet("../functions/makeEMSobjs_perp.jl") # data structures
includet("../functions/build_data.jl") # data structures
# includet("../functions/makeEMSplots.jl") # data structures
# includet("../functions/EMSfns_perp.jl") # General functions
# includet("../functions/ESSfns_perp.jl") # BESS functions
# includet("../functions/testEMS_perp.jl") # Test functions
# const labels = ["BNoDeg_W10", "CEmpDeg_W11", "CEmpDeg_W11e-2", "CPBDeg_W11", "CPBDeg_W11e-1", "CPBDeg_W11e-2"]
const labels = ["BNoDeg_W10", "CEmpDeg_W11e-2", "CPBDeg_W11e-2"]
# const labels = ["CPBDeg_W11", "CPBDeg_W11e-1", "CPBDeg_W11e-2"]

# load all the data
# inputs = build_data(;nEV=1, season="winter", profType="yearly", loadType="GV", year=2023)
results_DA = Dict(:summer=>Dict(), :winter=>Dict())
inputs = Dict(:summer=>Dict(), :winter=>Dict())
folder = "../data/output/Run24/"
steps = 29; cellID = "SYNSANYO"; profType = "weekly"; loadType = "GV"; year = 2023;
for season ∈ [:summer, :winter]
    inputs[season] = build_data(; nEV=1, 
            season=String(season),
            profType=profType,
            loadType=loadType,
            year=year,
            cellID=cellID)
    inputs[season]["grid"].λ /= 1000 # convert to €/kWh
    inputs[season]["grid"].λ /= 3600 # convert to €/kWs
    for DAtype ∈ labels
        results_DA[season][Symbol(DAtype)] = JSON3.read(open(folder*"cs1_R24_MPC_$(DAtype)_$(cellID)_y$(year)$(String(season))d$(steps).json", "r"),
                                                Dict{String, Vector{Dict}});
    end
end

for DAtype ∈ labels
    results_DA[Symbol(DAtype)] = JSON3.read(open(folder*"cs1_R21_DA_$(DAtype)_y2023d365.json", "r"), Vector{Dict});
end
# add the R20 with the warm-start
folder = "../data/output/Run20/"
results_DA[:CPBDeg_W11_wst] = JSON3.read(open(folder*"cs1_R20_DA_CPBDeg_W11_y2023d90.json", "r"), Vector{Dict});

# Function to create an Axis with common properties
function create_axis(figure, position, title)
    return Axis(figure[position];
        title = title,
        xlabel=L"$t$ [hr]", ylabel=L"$P_{sa}$ [kW]",
        ytickformat = values -> [L"%$(value)" for value in values],
        xtickformat = values -> [L"%$(value)" for value in values],
    )
end

function plotCompT(compTime; backend::String="GLMakie", filename::String="fig.png")
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    # Process the data
    # category_labels = ["BNoDeg_W10", "CEmpDeg_W11", "CPBDeg_W11", "CPBDeg_W11e4", "CPBDeg_W11_wst","CPBDeg_W10"]
    category_labels = ["BNoDeg_W10", "CEmpDeg_W11", "CPBDeg_W11", "CPBDeg_W11e4"]
    # data_array = vcat([log.(compTime[Symbol(DAtype)]) for DAtype ∈ category_labels]...)
    data_array = vcat([compTime[Symbol(DAtype)] for DAtype ∈ category_labels]...)
    category_labels = repeat(category_labels, inner = 29)

    f=Figure(size=(850, 600)) # create figure
    # colors=ColorSchemes.tab10.colors; # get colors
    colors=Makie.wong_colors();
    set_theme!(theme_latexfonts(), fontsize = 15) # set LaTeXStrings and fontsize
    
    # main axis
    ax1=Axis(f[1,1];
            # title = L"\textrm{Comp. time}",
            # xlabel="Planners", 
            # xlabelsize = 15,
            ylabel = L"$t_{\textrm{comp}}$ [s]",
            yscale = log10,
            yminorticksvisible = true,
            yminorgridvisible = true,
            # yminorticks = IntervalsBetween(10),
            )
    rainclouds!(ax1,category_labels, data_array;
        plot_boxplots = true, 
        clouds = hist,
        hist_bins = 2000,
        cloud_width = 0.9,
        boxplot_width = 0.15,
        gap = 0.2,
        dodge_gap = 1.5,
        # violin_limits = extrema,
        color = colors[indexin(category_labels, unique(category_labels))],
        )
        
    # inset axis
    inset_ax = Axis(f[2:3, 1],
        # width=Relative(0.9),
        # height=Relative(0.6),
        # halign=0.8,
        # valign=0.9,
        # backgroundcolor = :white,
        yscale = log10,
        ylabel = L"$t_{\textrm{comp}}$ [s]",
        xticks = (collect(1:length(unique(String.(category_labels)))), unique(String.(category_labels))),
        yminorticksvisible = true,
        yminorgridvisible = true,
        yminorticks = IntervalsBetween(5),
        )
    # translate!(inset_ax.scene, 0, 0, 10)
    # # this needs separate translation as well, since it's drawn in the parent scene
    # translate!(inset_ax.elements[:background], 0, 0, 9)
    # translate!(inset_ax.elements[:xgridlines], 0, 0, 9.5)
    # translate!(inset_ax.elements[:ygridlines], 0, 0, 9.5)
    rainclouds!(inset_ax,category_labels, data_array;
        plot_boxplots = true, 
        clouds = hist,
        hist_bins = 2000,
        cloud_width = .9,
        boxplot_width = 0.15,
        gap = 0.2,
        dodge_gap = 1.5,
        color = colors[indexin(category_labels, unique(category_labels))])
    limits!(inset_ax, nothing, nothing, 10^(1.88), 10^2.4)
    hidexdecorations!(ax1, grid = false)
    # axislegend(ax1, position=:rb)
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function gridHist(results_DA, inputs; backend::String="GLMakie", filename::String="fig.png")
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end
    grid = []
    valueCap = []
    # Process the data
    # category_labels = ["BNoDeg_W10", "CEmpDeg_W11", "CPBDeg_W11", "CPBDeg_W11e4", "CPBDeg_W11_wst","CPBDeg_W10"]
    # for (i, DAtype) ∈ enumerate(category_labels)
    for (i, DAtype) ∈ enumerate(labels)
        concRes = concatResultsRH(results_DA[Symbol(DAtype)]; typeOpt="day-ahead")
        Pg = convert.(Float64, concRes["Pg"])
        λ = inputs["grid"].λ[1:length(Pg),:]*1e-3
        value = zeros(length(Pg))
        value[Pg .> 0] .= λ[Pg .> 0, 1] .* Pg[Pg .> 0] # imports
        value[Pg .< 0] .= λ[Pg .< 0, 2] .* Pg[Pg .< 0] # exports
        push!(grid, Pg)
        push!(valueCap, value)
    end
    grid = vcat(grid...)
    valueCap = vcat(valueCap...)
    # data_array = vcat.([results_DA[Symbol(DAtype)][:"Pg"] for DAtype ∈ category_labels])
    # category_labels = repeat(category_labels, inner = Int(length(grid)/length(category_labels)))
    category_labels = repeat(labels, inner = Int(length(grid)/length(labels)))

    f=Figure(size=(850, 600)) # create figure
    # colors=ColorSchemes.tab10.colors; # get colors
    colors=Makie.wong_colors();
    set_theme!(theme_latexfonts(), fontsize = 15) # set LaTeXStrings and fontsize
    
    # main axis
    ax1=Axis(f[1,1];
            # title = L"\textrm{Comp. time}",
            # xlabel="Planners", 
            # xlabelsize = 15,
            ylabel = L"$\lambda_{\textrm{buy/sell}} P_{\textrm{g}}$ [€/h]",
            # yscale = log10,
            yminorticksvisible = true,
            yminorgridvisible = true,
            xticks = (collect(1:length(unique(String.(category_labels)))),
                     unique(String.(category_labels))),
            # yminorticks = IntervalsBetween(10),
            )
    # rainclouds!(ax1,category_labels, grid;
    rainclouds!(ax1,category_labels, valueCap;
        plot_boxplots = true, 
        # clouds=hist,
        cloud_width=1,
        boxplot_width=0.15,
        gap = 0.02,
        dodge_gap = 1.5,
        # violin_limits = extrema,
        color = colors[indexin(category_labels, unique(category_labels))],
        )
    # limits!(inset_ax, nothing, nothing, 10^(1.7), 10^2.2)
    # hidexdecorations!(ax1, grid = false)
    # # axislegend(ax1, position=:rb)
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function create_Pfig(results_DA, P_key, label; backend::String="GLMakie", filename::String="fig.png")
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end
    
    set_theme!(theme_latexfonts(), fontsize = 20)
    f=Figure(size=(800, 800/1.618))
    # colors=ColorSchemes.tab10.colors;
    colors=Makie.wong_colors();
    
    ax1=Axis(f[1,1];
            ylabel=label,
            xticks = 0:24*3:24*29,
            xautolimitmargin = (0, 0), yautolimitmargin = (0.1, 0.15),
            )
    ax2=Axis(f[2,1];
            xlabel=L"$t$ [hr]", ylabel=label,
            xticks = 0:24*3:24*29,
            xautolimitmargin = (0, 0), yautolimitmargin = (0.1, 0.15),
            )
    ax1_inset = Axis(f[1, 1],
                width=Relative(0.45/1.618),
                height=Relative(0.45),
                halign=.4, valign = .25,
                # yticks = 0.2:0.2:1,
                yaxisposition = :left,
                xautolimitmargin = (0, 0), yautolimitmargin = (0.05, 0.25),
                )
    # ax2_inset = Axis(f[2, 1],
    #             width=Relative(0.45/1.618),
    #             height=Relative(0.45),
    #             halign=.85, valign = 0.15,
    #             # yticks = 0.1:0.03:16,
    #             yaxisposition = :right,
    #             xautolimitmargin = (0.1, 0.1), yautolimitmargin = (0.1, 0.1),
    #             )

    # axes = [ax1 ax2;
    #         ax1_inset ax2_inset];
    axes = [ax1 ax2;
            ax1_inset ax1_inset];
    #loop over seasons
    for (j,season) ∈ enumerate([:summer, :winter])
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            lbl = split(DAtype, "_")[1]
            concRes = concatResultsFlex(results_DA[season][Symbol(DAtype)]["sim"], CT_MPC())
            t = convert.(Float64, concRes["t"])/3600
            Psa = convert.(Float64, concRes[P_key])
            γ_cont = concRes[:"γf"][:]
            j == 1 ? day = 12 : day = 14;
            it0_zoom = day*4*24
            itend_zoom = it0_zoom + 24*4*2 - 1;
            axes[2,j].xticks = t[it0_zoom]:6:t[itend_zoom]
            # stairs!(axes[1,j], t, Psa, color=colors[i], label="$(DAtype[1:end-4])", linewidth=2, step=:post)
            stairs!(axes[1,j], t, Psa, color=colors[i], label= lbl, linewidth=2, step=:post)
            if j == 1
                stairs!(axes[2,j], t[it0_zoom:itend_zoom], Psa[it0_zoom:itend_zoom], color=colors[i], linewidth=2, step=:post)
                limits!(axes[2,j], t[it0_zoom], t[itend_zoom], nothing, nothing);
            end
            
            # iArr = findall(x -> x == 1, diff(γ_cont))
            # # departure when γ_cont changes from 1 to 0
            # iDep = findall(x -> x == -1, diff(γ_cont))
            # # check if they are the same length.
            # iDep = Int.(iDep); iArr = Int.(iArr);
            # vspan!(axes[1,j], t[iDep],t[iArr], ymax=100, color = (:grey, 0.15))
            # vspan!(axes[2,j], t[iDep],t[iArr], ymax=100, color = (:grey, 0.15))
            
        end
    end
    # translate!.([ax1_inset.blockscene, ax2_inset.blockscene], 0, 0, 150)
    translate!.([ax1_inset.blockscene], 0, 0, 150)

    linkxaxes!(ax1, ax2);
    # hidexdecorations!.([ax1, ax1_inset, ax2_inset], grid = false, ticks = false);
    hidexdecorations!.([ax1, ax1_inset], grid = false, ticks = false);
    axislegend(ax1, position=:rb, size = 12)
    #add a text annotation
    textlabel!(axes[1,1], 450, .1, text= "Summer", fontsize = 20)
    textlabel!(axes[1,2], 150, .35, text= "Winter", fontsize = 20)
    textlabel!(axes[2,1], 322, 0.94, text= "Days 12-13", fontsize = 18)
    # textlabel!(axes[2,2], 24*15, 0.22, text= "Days 15-16", fontsize = 18)
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function compareFEC(results_DA,
        i_key,
        Q_key;
        backend::String="GLMakie",
        filename::String="fig.png",
        step::Int=850,
        )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    set_theme!(theme_latexfonts(), fontsize = 20)
    f=Figure(size=(900, 700))
    # colors=ColorSchemes.tab10.colors;
    colors=Makie.wong_colors();
    markers = [:circle, :diamond];
    sa_labels = ["BESS", "EV"];
    # loop over seasons
    axs=[]
    for (j,season) ∈ enumerate([:summer, :winter])
        ax1=Axis(f[j,1];
            xlabel=L"$t$ [hr]", ylabel=L"$FEC$",
            )
        ax2=Axis(f[j,2];
            xlabel=L"$FEC$", ylabel=L"$Q_{\textrm{loss}}$ [%]",
            )
        axes = [ax1, ax2]
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            concRes = concatResultsFlex(results_DA[season][Symbol(DAtype)]["sim"], CT_MPC())
            t = convert.(Float64, concRes["t"]); Δt = t[2] - t[1];
            for sa ∈ 1:length(i_key)
                isa = convert.(Float64, concRes[i_key[sa]])
                Qsa = convert.(Float64, concRes[Q_key[sa]])
                Qsa0 = Qsa[1]
                FECsa = cumsum(abs.(isa)) .* Δt / 3600 / (2 * Qsa0) 
                lines!(axes[1], t./3600, FECsa,
                    color=colors[i], linewidth=2)
                scatter!(axes[1], t[1:step:end]./3600, FECsa[1:step:end],
                    color=(colors[i], 0.5), marker=markers[sa],
                    strokewidth = 1, strokecolor = colors[i], markersize = 12)
                lines!(axes[2], FECsa, 100*(Qsa .- Qsa0)/Qsa0,
                    color=colors[i], linewidth=2, label="$DAtype")
                scatter!(axes[2], FECsa[1:step:end], 100*((Qsa .- Qsa0)/Qsa0)[1:step:end],
                    color=(colors[i], 0.5), marker=markers[sa],
                    strokewidth = 1, strokecolor = colors[i], markersize = 12, label="$(sa_labels[sa])")
            end
            limits!(ax1, t[1]/3600, t[end]/3600, nothing, nothing);
        end
        # j == 1 ? axislegend(ax2, position=:rt, merge = true) : nothing;
        push!(axs, axes)
    end
    linkxaxes!(axs[1][1], axs[2][1])
    linkxaxes!(axs[1][2], axs[2][2])
    hidexdecorations!.(axs[1][:], grid = false)
    
    #add a text annotation
    textlabel!(f[1,1], 200, 25, text= "Summer", fontsize = 20)
    textlabel!(f[2,1], 200, 25, text= "Winter", fontsize = 20)
    group_type = [PolyElement(color = colors[j], strokecolor = :transparent) for j ∈ 1:3]
    group_sa = [MarkerElement(marker = mkr, color = (:gray, 0.5), strokecolor = :gray, markersize = 15) for mkr ∈ markers]
    # println([group_type[:] group_sa[:]])
    legends = Legend(f[1,2],
        [group_type, group_sa],
        [["BNoDeg", "CEmpDeg", "CPBDeg"], ["BESS", "EV"]],
        ["Type","Device"],
        tellheight = false, tellwidth = false,
        halign = :right, valign = :top,
        # nbanks = 2,
        # orientation = :horizontal,
        margin = (5,5,5,5),
        )

    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function compareFEC_wLoss(results_DA,
        i_key,
        Q_key;
        backend::String="GLMakie",
        filename::String="fig.png",
        step::Int=850,
        )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    set_theme!(theme_latexfonts(), fontsize = 20)
    f=Figure(size=(900, 700))
    # colors=ColorSchemes.tab10.colors;
    colors=Makie.wong_colors();
    markers = [:circle, :diamond];
    sa_labels = ["BESS", "EV"];
    # loop over seasons
    axs=[]
    for (j,season) ∈ enumerate([:summer, :winter])
        ax1=Axis(f[j,1];
            xlabel=L"$t$ [hr]", ylabel=L"$FEC$",
            )
        ax2=Axis(f[j,2:3];
            xlabel=L"$FEC$", ylabel=L"$Q_{\textrm{loss}}$ [%]",
            yaxisposition = :right,
            )
        axes = [ax1, ax2]
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            concRes = concatResultsFlex(results_DA[season][Symbol(DAtype)]["sim"], CT_MPC())
            t = convert.(Float64, concRes["t"]); Δt = t[2] - t[1];
            for sa ∈ 1:length(i_key)
                isa = convert.(Float64, concRes[i_key[sa]])
                Qsa = convert.(Float64, concRes[Q_key[sa]])
                Qsa0 = Qsa[1]
                FECsa = cumsum(abs.(isa)) .* Δt / 3600 / (2 * Qsa0) 
                lines!(axes[1], t./3600, FECsa,
                    color=(colors[3], 1/i), linewidth=2)
                scatter!(axes[1], t[1:step:end]./3600, FECsa[1:step:end],
                    color=(colors[3], 0.5), marker=markers[sa],
                    strokewidth = 1, strokecolor = (colors[3], 1/i), markersize = 12)
                lines!(axes[2], FECsa, 100*(Qsa .- Qsa0)/Qsa0,
                    color=(colors[3], 1/i), linewidth=2, label="$DAtype")
                scatter!(axes[2], FECsa[1:step:end], 100*((Qsa .- Qsa0)/Qsa0)[1:step:end],
                    color=(colors[3], 0.5), marker=markers[sa],
                    strokewidth = 1, strokecolor = (colors[3], 1/i), markersize = 12, label="$(sa_labels[sa])")
            end
            limits!(ax1, t[1]/3600, t[end]/3600, nothing, nothing);
        end
        # j == 1 ? axislegend(ax2, position=:rt, merge = true) : nothing;
        push!(axs, axes)
    end
    linkxaxes!(axs[1][1], axs[2][1])
    linkxaxes!(axs[1][2], axs[2][2])
    hidexdecorations!.(axs[1][:], grid = false)
    
    #add a text annotation
    textlabel!(f[1,1], 200, 25, text= "Summer", fontsize = 20)
    textlabel!(f[2,1], 200, 25, text= "Winter", fontsize = 20)
    group_type = [PolyElement(color = (colors[3], 1/j), strokecolor = :transparent) for j ∈ 1:3]
    group_sa = [MarkerElement(marker = mkr, color = (:gray, 0.5), strokecolor = :gray, markersize = 15) for mkr ∈ markers]
    legends = Legend(f[1,2:3],
        [group_type, group_sa],
        [["1", "0.1", "0.01"], ["BESS", "EV"]],
        [L"w_{\text{loss}}","Device"],
        tellheight = false, tellwidth = false,
        halign = :right, valign = :top,
        # nbanks = 2,
        # orientation = :horizontal,
        margin = (5,5,5,5),
        )

    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function create_SPfig(results_DA, key;
    backend::String="GLMakie",
    filename::String="fig.png",
    cmaps::Vector{Symbol}=[:dense, :amp],
    cutoff::Float64=3.0, # limit between resting and cycling
    lims=([-7.5, 7.5], [0.2, 0.8]) # xmin xmax ymin ymax
    )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    set_theme!(theme_latexfonts())
    f=Figure(size=(1000, 400))
    axes=[];
    co=[]; co_ex=[];
    
    # setup keys
    SoC_key = "SoC$key";
    key == "ev[1]" ? key = "evTot[1]" : nothing;
    P_key = "P$key";
    
    for (j, season) ∈ enumerate(keys(results_DA))
        results = results_DA[season]
        # Initialize the histograms
        # Resting 2d-histogram
        histogram = Vector{Matrix{Float64}}(undef, length(results));
        xedges = Vector{Vector{Float64}}(undef, length(results));
        yedges = Vector{Vector{Float64}}(undef, length(results));
        # Cycling 2d-histogram
        histogram_ex = Vector{Matrix{Float64}}(undef, length(results));
        xedges_ex = Vector{Vector{Float64}}(undef, length(results));
        yedges_ex = Vector{Vector{Float64}}(undef, length(results));

        # create the histograms
        for (i, DAtype) ∈ enumerate(labels)
            concRes = concatResultsFlex(results[Symbol(DAtype)]["sim"], CT_MPC())
            Psa = convert.(Float64, concRes[P_key])
            SoCsa = convert.(Float64, concRes[SoC_key])
            if key == "evTot[1]"
                # exclude the points where the EV is not connected
                γ_cont = concRes[:"γf"][:]
                Psa = Psa[γ_cont .== 1];
                SoCsa = SoCsa[γ_cont .== 1];
            end
            # Create a 2D histogram
            histogram[i], xedges[i], yedges[i] = hist2d(Psa, SoCsa; nbins=30, normalize=true)
            histogram_ex[i], xedges_ex[i], yedges_ex[i] = hist2d(Psa[abs.(Psa) .≥ cutoff], SoCsa[abs.(Psa) .≥ cutoff]; nbins=30, normalize=true)
            # append!(histogram, hist);
        end
        
        # get the global limits
        # Resting 2d-histogram
        extremas = map(extrema, histogram) # get the extrema of each histograms
        global_min = minimum(t->first(t), extremas)
        global_max = maximum(t->last(t), extremas)
        clims = (global_min, global_max) # these limits have to be shared by the maps and the colorbar
        # Cycling 2d-histogram
        extremas_ex = map(extrema, histogram_ex) # get the extrema of each histograms
        global_min_ex = minimum(t->first(t), extremas_ex)
        global_max_ex = maximum(t->last(t), extremas_ex)
        clims_ex = (global_min_ex, global_max_ex) # these limits have to be shared by the maps and the colorbar
        
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            ax = Axis(f[j,i+1];
                xlabel=L"$P$ [kW]", ylabel=L"$SoC$ [p.u.]",
                )
            j == 1 ? ax.title = "$DAtype" : nothing;
            # Create a 2D histogram
            # histogram, xedges, yedges = hist2d(Psa, SoCsa; nbins=60, normalize=true)
            cf = contourf!(ax, xedges[i][2:end], yedges[i][2:end], histogram[i],
                        levels=0.01:0.1:1, mode = :relative,
                        colormap=cmaps[1], 
                        # colorrange=clims,
                        label = "$DAtype")
            cf_ex = contourf!(ax, xedges_ex[i][2:end], yedges_ex[i][2:end], histogram_ex[i],
                        levels=0.01:0.1:1, mode = :relative,
                        colormap=cmaps[2],
                        # colorrange=clims_ex,
                        label = "$DAtype")
            limits!(ax, lims[1][1], lims[1][2], lims[2][1], lims[2][2]);
            # limits!(ax, -7.5, 7.5, 0.2, 0.8);
            j == 1 ? hidexdecorations!(ax, grid = false) : nothing;
            # j == 2 ? Colorbar(f[j+1,i], co, vertical = false) : nothing;
            # Colorbar(f[2,i], co, vertical = false, limits=clims)
            i == 1 ? nothing : hideydecorations!(ax, grid = false);
            push!(axes, ax)
            push!(co, cf)
            push!(co_ex, cf_ex)
        end
    end
    
    # add annotation with seasons for each row
    textlabel!(f[1,4], 1, 0.25, text= "Summer", fontsize = 17)
    textlabel!(f[2,4], 1, 0.25, text= "Winter", fontsize = 17)
    Colorbar(f[1:2,1], co_ex[end], label = "Cycling", flipaxis = false)
    Colorbar(f[1:2,length(results_DA[:winter])+2], co[end], label = "Resting")
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function costSummaryPlot(results_DA, costSummary_DA;
    backend::String="GLMakie",
    filename::String="fig.png",
    )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end
    f=Figure(size=(600, 400))
    # colors=ColorSchemes.Dark2_7.colors;
    colors=Makie.wong_colors();
    set_theme!(theme_latexfonts())
    ax1=Axis(f[1,1];
            title = L"\textrm{Summer}",
            xlabel=L"$t$ [hr]", ylabel=L"$Cost$ [€]",
            )
    # ax2=Axis(f[2,1];
    #         title = L"\textrm{Winter}",
    #         xlabel=L"$t$ [hr]", ylabel=L"$Cost$ [€]",
    #         )
    # axes = [ax1, ax2]
    # # loop over the season
    # for (j,season) ∈ enumerate([:summer, :winter])
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            concRes = concatResultsRH(results_DA[Symbol(DAtype)]; typeOpt="day-ahead")
            t = convert.(Float64, concRes["t"])
            stairs!(ax1, t./3600, costSummary_DA[Symbol(DAtype)].wCgrid, color=colors[i],
                step=:post, linewidth=2, label="$DAtype")
            stairs!(ax1, t./3600, costSummary_DA[Symbol(DAtype)].wCloss, color=colors[i],
                step=:post, linewidth=2, linestyle=:dash, label="$DAtype")
            limits!(ax1, t[1]/3600, t[end]/3600, nothing, nothing);
        end
    # end
    axislegend(ax1, position=:lt, merge=true, orientation=:vertical)
    # hidexdecorations!(ax1, grid=false)
    # Add arrows/annotations to ax2
    textlabel!(ax1, 1300, 0.45, text= L"$w_{\textrm{grid}} C_{\textrm{grid}}$", fontsize = 20)
    textlabel!(ax1, 1500, 0.25, text= L"$w_{\textrm{loss}} C_{\textrm{loss}}$", fontsize = 20)
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function dailyCgridPlot(dailyCgrid;
    backend::String="GLMakie",
    filename::String="fig.png",
    )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end
    
    f=Figure(size=(600, 400))
    colors=Makie.wong_colors();
    set_theme!(theme_latexfonts())

    # loop over the DAtype
    for (j,season) ∈ enumerate([:summer, :winter])
        ax=Axis(f[j,1];
            xlabel=L"$t$ [hr]", ylabel=L"$Cost$ [€]",
            )
        lessExpensive = [findmin([dailyCgrid[season][Symbol(DAtype)][day] for DAtype ∈ labels]) for day ∈ 1:29]
        lessExpensive_labels = [labels[lessExpensive[i][2]] for i ∈ 1:29]
        println("In, ", season, " the cheapest DAtype is:")
        println("BNoDeg wins in ", sum(lessExpensive_labels .== "BNoDeg_W10"), " days")
        println("CEmpDeg wins in ", sum(lessExpensive_labels .== "CEmpDeg_W11"), " days")
        println("CPBDeg_w11 wins in ", sum(lessExpensive_labels .== "CPBDeg_W11"), " days")
        println("CPBDeg_w11e4 wins in ", sum(lessExpensive_labels .== "CPBDeg_W11e4"), " days")
        for st ∈ 1:29
            vspan!(ax, st-0.5, st+0.5; color=(colors[lessExpensive[st][2]], 0.4), label = lessExpensive_labels[st])
        end
        for DAtype ∈ labels
            lines!(ax, 1:29, dailyCgrid[season][Symbol(DAtype)],
                linewidth=2, label=DAtype)
        end
        Legend(f[j,2], ax, position=:lt, merge=true)
    end
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function dailyCgridBoxplot(dailyCgrid, s::modelSettings;
                        days::Array{DateTime,1}=DateTime.(Date(2023,1,1):Day(1):Date(2023,12,31)),
                        backend::String="GLMakie",
                        filename::String="fig.png",
                        labels=["BNoDeg_W10", "CEmpDeg_W11e4", "CPBDeg_W11", "CPBDeg_W11e4", "CPBDeg_W11_wst", "CPBDeg_W10"],
                        )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    # extract the months of the days array
    months = month.(days);
    # cut to the first 180 days
    days = days[1:s.steps]; months = months[1:s.steps];
    # create the data for the boxplot
    dCg = vcat([dailyCgrid[Symbol(DAtype)] for DAtype ∈ labels]...);
    months = repeat(months, outer = length(labels));
    dodges = repeat(1:length(labels), inner = s.steps);
    colors = Makie.wong_colors()[1:length(labels)]; colors = repeat(colors, inner = s.steps);

    # create a boxplot
    f=Figure(size=(600, 400))
    ax = Axis(f[1,1]; xlabel=L"\textrm{month}", ylabel=L"$c_{g}$ [€/day]")
    boxplot!(ax, months, float.(dCg), dodge = dodges, show_notch = true, color = colors, label = labels)
    # build the marker elements
    elems = [[MarkerElement(color = col, marker=:circle, markersize = 15,
          strokecolor = :black)] for col in unique(colors)]
    axislegend(ax, elems, labels; position=:rb)
    
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function compareSimCtrl(results_DA::Dict, SoC_key::AbstractString, label::AbstractString;
    backend::String="GLMakie", filename::String="fig.pdf"
    )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    set_theme!(theme_latexfonts(), fontsize = 20)
    f=Figure(size=(800, 800/1.618))
    # colors=ColorSchemes.tab10.colors;
    colors=Makie.wong_colors();
    
    #loop over seasons
    for (j,season) ∈ enumerate([:summer, :winter])
        ax1=Axis(f[1,j];
            title = L"BNoDeg",
            ylabel=label,
            xticks = 0:24*3:24*29,
            xautolimitmargin = (0, 0), yautolimitmargin = (0.1, 0.15),
            )
        ax2=Axis(f[2,j];
                title = L"CEmpDeg",
                ylabel=label,
                xticks = 0:24*3:24*29,
                xautolimitmargin = (0, 0), yautolimitmargin = (0.1, 0.15),
                )
        ax3=Axis(f[3,j];
                title = L"CPBDeg",
                xlabel=L"$t$ [hr]", ylabel=label,
                xticks = 0:24*3:24*29,
                xautolimitmargin = (0, 0), yautolimitmargin = (0.1, 0.15),
                )

        axes = [ax1 ax2 ax3];
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(labels)
            concRes = concatResultsRH(results_DA[season][Symbol(DAtype)]["sim"]; typeOpt="day-ahead")
            concCtrl = concatResultsRH(results_DA[season][Symbol(DAtype)]["ctr"]; typeOpt="day-ahead")
            t = convert.(Float64, concRes["t"])/3600/24
            SoCsa_pos = convert.(Float64, concRes[SoC_key])
            SoCsa_pre = convert.(Float64, concCtrl[SoC_key])
            γ_cont = concRes[:"γ_cont"][:]
            stairs!(axes[i], t, SoCsa_pre, color=(colors[i],0.5), label=L"$\tilde{S}_{b,t}$", linewidth=2, step=:post, linestyle=:dash)
            stairs!(axes[i], t, SoCsa_pos, color=colors[i], label=L"$S_{b,t}$", linewidth=2, step=:post)
            
            iArr = findall(x -> x == 1, diff(γ_cont))
            # departure when γ_cont changes from 1 to 0
            iDep = findall(x -> x == -1, diff(γ_cont))
            # check if they are the same length.
            iDep = Int.(iDep); iArr = Int.(iArr);
            vspan!(axes[1], t[iDep],t[iArr], ymax=100, color = (:grey, 0.15))
            vspan!(axes[2], t[iDep],t[iArr], ymax=100, color = (:grey, 0.15))
            vspan!(axes[3], t[iDep],t[iArr], ymax=100, color = (:grey, 0.15))
        end
        linkxaxes!(ax1, ax2, ax3);
        hidexdecorations!.([ax1, ax2], grid = false, ticks = false);
        axislegend(ax1, position=:rb, size = 12)
        #add a text annotation
        textlabel!(axes[1], 450/24, .1, text= String(season), fontsize = 20)
        # textlabel!(axes[1,2], 150, .35, text= "Winter", fontsize = 20)
        # textlabel!(axes[2,1], 290, 0.83, text= "Days 12-13", fontsize = 18)
        # textlabel!(axes[2,2], 24*15, 0.14, text= "Days 15-16", fontsize = 18)

    end
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end
end

function calcDailyCost(result, data, s)
    t=result["t"]; # time
    Δt=t[2]-t[1]; # time step [s]
    t0=t[1];
    # make new time window
    it0 = Int((t0/Δt) + 1);
    itend = it0+length(t)-1;
    λbuy = data["grid"].λ[it0:itend,1]; # buy price [€/kWs]
    λsell = data["grid"].λ[it0:itend,2]; # sell price [€/kWs]
    Pgpos = result[:"PgPos"]; # positive grid power
    Pgneg = result[:"PgNeg"]; # negative grid power

    dailyCost = []
    for day ∈ 1:s.steps
        # samples per day
        samples = Int.(length(t) / s.steps);
        # time index for the day
        idt = (1:samples) .+ (day-1)*samples;
        cg = sum(Pgpos[idt].*λbuy[idt] - Pgneg[idt].*λsell[idt]).*Δt;
        push!(dailyCost, cg)
    end
    return dailyCost
end

function plotDispatchAnalysis(results_DA::Dict,
        inputs::Dict;
        backend::String="CairoMakie",
        filename::String="dispatchAnalysis.pdf"
        )
        @assert backend ∈ ["CairoMakie", "GLMakie"] "Invalid backend. Choose CairoMakie or GLMakie."
        if backend == "CairoMakie"
                CairoMakie.activate!(type="svg")
        else
                GLMakie.activate!()
        end

        set_theme!(theme_latexfonts())
        fig = Figure(size=(1000, 650), fontsize=20)
        colors=ColorSchemes.Dark2_7.colors;
        axs = Matrix{Axis}(undef, 3, 2)
        for (j, season) in enumerate(keys(results_DA))
                # Unwrap the model
                results = concatResultsFlex(results_DA[season][Symbol("CPBDeg_W11e-2")]["sim"], CT_MPC())
                data = inputs[season]
                t = results[:"t"]
                t0 = t[1]
                Δt = t[2] - t[1]
                it0 = round(Int, (t0 / Δt))
                itend = it0 + length(t) - 1
                # Inputs
                spvModel = data[:"SPV"];
                gridModel = data[:"grid"];
                # Results
                Pg = convert.(Float64, results[:"Pg"])
                haskey(results, "Phpe") ? Phpe = convert.(Float64, results[:"Phpe"]) : nothing
                Pbess = convert.(Float64, results[:"Pbess"])
                SoCbess = convert.(Float64, results[:"SoCbess"])
                nEV = length(data["EV"])
                if length(data["EV"]) !=1
                        γ_cont = results[:"γf"]
                        Pev = [results[:"Pev[1]"], results[:"Pev[2]"]]
                        SoCev = [results[:"SoCev[1]"], results[:"SoCev[2]"]]
                else
                        γ_cont = results[:"γf"][:]
                        Pev = convert.(Float64, results[:"Pev[1]"])
                        SoCev = convert.(Float64, results[:"SoCev[1]"])
                end
                Ptess = convert.(Float64, results[:"Ptess"])
                SoCtess = convert.(Float64, results[:"SoCtess"])
                
                # Create the axes
                eb_ax = Axis(fig[1, j]);
                tb_ax = Axis(fig[2,j]);
                eess_ax = Axis(fig[3, j];
                        xlabel=L"$t$ [hr]",        
                        );
                axs[1,j] = eb_ax; axs[2,j] = tb_ax; axs[3,j] = eess_ax;
                # inset axes, zoom in the summer for day 12 and 23 in the winter
                j == 1 ? day = 12 : day = 23;
                it0_zoom = day*4*24
                itend_zoom = it0_zoom + 24*4 - 1;
                eb_ax_inset = Axis(fig[1, j],
                        width=Relative(0.55/1.618),
                        height=Relative(0.55),
                        xticks = t[it0_zoom]:6:t[itend_zoom],
                        halign=.3, valign = -0.25,
                        yaxisposition = :right,
                        xautolimitmargin = (0, 0), yautolimitmargin = (0.05, 0.05),
                        )
                tb_ax_inset = Axis(fig[2, j],
                        width=Relative(0.55/1.618),
                        height=Relative(0.55),
                        xticks = t[it0_zoom]:6:t[itend_zoom],
                        halign=0.8, valign = 1.05,
                        yaxisposition = :right,
                        xautolimitmargin = (0, 0), yautolimitmargin = (0.05, 0.05),
                        )
                eess_ax_inset = Axis(fig[3, j],
                        width=Relative(0.55/1.618),
                        height=Relative(0.55),
                        xticks = t[it0_zoom]:6:t[itend_zoom],
                        halign=.8, valign = 0.1,
                        yaxisposition = :right,
                        xautolimitmargin = (0, 0), yautolimitmargin = (0.05, 0.05),
                        )
                # Energy Balance plots
                stairs!(eb_ax, t/3600, gridModel.loadE[it0:itend], color=colors[1], linewidth=1.5, label=L"P_{\textrm{load}}^{e}", step=:post)
                stairs!(eb_ax, t/3600, spvModel.MPPTData[it0:itend], color=colors[2], linewidth=1.5, label=L"P_{\textrm{PV}}", step=:post)
                stairs!(eb_ax, t/3600, Pg, color=colors[3], linewidth=1.5, label=L"P_{\textrm{g}}", step=:post)
                if haskey(results, "Phpe")
                        stairs!(eb_ax, t/3600, Phpe, color=colors[4], linewidth=1.5, label=L"P_{\textrm{HP}}^{e}", step=:post)
                end
                stairs!(eb_ax, t/3600, Pbess, color=colors[5], linewidth=1.5, label=L"P_{\textrm{BESS}}", step=:post)
                if length(data["EV"]) != 1
                        [stairs!(eb_ax, t/3600, Pev[n] .* γ_cont[n], color=colors[5+n], linewidth=1.5, label=L"P_{\textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
                else
                        stairs!(eb_ax, t/3600, Pev .* γ_cont, color=colors[6], linewidth=1.5, label=L"P_{\textrm{EV}}", step=:post)
                end
                limits!(eb_ax, t[1]/3600, t[end]/3600, nothing, nothing);
                # zoom In
                stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, gridModel.loadE[it0_zoom:itend_zoom], color=colors[1], linewidth=1.5, step=:post)
                stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, spvModel.MPPTData[it0_zoom:itend_zoom], color=colors[2], linewidth=1.5, step=:post)
                stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, Pg[it0_zoom:itend_zoom], color=colors[3], linewidth=1.5, step=:post)
                if haskey(results, "Phpe")
                        stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, Phpe[it0_zoom:itend_zoom], color=colors[4], linewidth=1.5, step=:post)
                end
                stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, Pbess[it0_zoom:itend_zoom], color=colors[5], linewidth=1.5, step=:post)
                if length(data["EV"]) != 1
                        [stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, Pev[n][it0_zoom:itend_zoom] .* γ_cont[n][it0_zoom:itend_zoom], color=colors[5+n], linewidth=1.5, step=:post) for n ∈ 1:nEV]
                else
                        stairs!(eb_ax_inset, t[it0_zoom:itend_zoom]/3600, Pev[it0_zoom:itend_zoom] .* γ_cont[it0_zoom:itend_zoom], color=colors[6], linewidth=1.5, step=:post)
                end

                # Thermal Balance plot
                stairs!(tb_ax, t/3600, gridModel.loadTh[it0:itend], color=colors[1], linewidth=1.5, label=L"P_{\textrm{load}}^{th}", step=:post)
                stairs!(tb_ax, t/3600, spvModel.MPPTData[it0:itend]*data["ST"].η, color=colors[2], linewidth=1.5, label=L"P_{\textrm{ST}}", step=:post)
                stairs!(tb_ax, t/3600, Ptess, color=colors[3], linewidth=1.5, label=L"P_{\textrm{TESS}}", step=:post)
                if haskey(results, "Phpe")
                        stairs!(tb_ax, t/3600, Phpe*data["HP"].η, color=colors[4], linewidth=1.5, label=L"P_{\textrm{HP}}^{t}", step=:post)
                end
                limits!(tb_ax, t[1]/3600, t[end]/3600, nothing, nothing);
                # zoom In
                stairs!(tb_ax_inset, t[it0_zoom:itend_zoom]/3600, gridModel.loadTh[it0_zoom:itend_zoom], color=colors[1], linewidth=1.5, step=:post)
                stairs!(tb_ax_inset, t[it0_zoom:itend_zoom]/3600, spvModel.MPPTData[it0_zoom:itend_zoom]*data["ST"].η, color=colors[2], linewidth=1.5, step=:post)
                stairs!(tb_ax_inset, t[it0_zoom:itend_zoom]/3600, Ptess[it0_zoom:itend_zoom], color=colors[3], linewidth=1.5, step=:post)
                if haskey(results, "Phpe")
                        stairs!(tb_ax_inset, t[it0_zoom:itend_zoom]/3600, Phpe[it0_zoom:itend_zoom]*data["HP"].η, color=colors[4], linewidth=1.5, step=:post)
                end
                
                # Electric Storage plot
                stairs!(eess_ax, t/3600, SoCbess*100, color=colors[1], linewidth=1.5, label=L"SoC_{\textrm{BESS}}", step=:post)
                if length(data["EV"]) != 1
                        [stairs!(eess_ax, t/3600, SoCev[n] .* γ_cont[n] * 100, 
                        color=colors[1+n], linewidth=1.5, label=L"SoC_{\textrm{EV}, %$n}", step=:post) for n ∈ 1:nEV]
                else
                        stairs!(eess_ax, t/3600, SoCev .* γ_cont .* 100, 
                        color=colors[1+1], linewidth=1.5, label=L"SoC_{\textrm{EV}}", step=:post)
                        # add a vspan for the EV availability
                        # find the indeces were γ_cont changes from 0 to 1 and viceversa
                        # arrival when γ_cont changes from 0 to 1
                        iArr = findall(x -> x == 1, diff(γ_cont))
                        # departure when γ_cont changes from 1 to 0
                        iDep = findall(x -> x == -1, diff(γ_cont))
                        # check if they are the same length.
                        iDep = Int.(iDep); iArr = Int.(iArr);
                        vspan!(eess_ax, t[iDep]/3600,t[iArr]/3600, ymax=100, color = (:grey, 0.2))
                end
                stairs!(eess_ax, t/3600, SoCtess * 100, color=colors[3], linewidth=1.5, label=L"SoC_{\textrm{TESS}}", step=:post)
                # zoom In
                stairs!(eess_ax_inset, t[it0_zoom:itend_zoom]/3600, SoCbess[it0_zoom:itend_zoom]*100, color=colors[1], linewidth=1.5, step=:post)
                stairs!(eess_ax_inset, t/3600, SoCev .* γ_cont .* 100, color=colors[2], linewidth=1.5, step=:post)
                stairs!(eess_ax_inset, t[it0_zoom:itend_zoom]/3600, SoCtess[it0_zoom:itend_zoom] * 100, color=colors[3], linewidth=1.5, step=:post)
                # add a vspan for the EV availability
                # find the indeces were γ_cont changes from 0 to 1 and viceversa
                # arrival when γ_cont changes from 0 to 1
                iArr = findall(x -> x == 1, diff(γ_cont))
                # departure when γ_cont changes from 1 to 0
                iDep = findall(x -> x == -1, diff(γ_cont))
                # check if they are the same length.
                iDep = Int.(iDep); iArr = Int.(iArr);
                vspan!(eess_ax_inset, t[iDep]/3600,t[iArr]/3600, ymax=100, color = (:grey, 0.2))

                if j == 2
                        Legend(fig[1, j+1], eb_ax)
                        Legend(fig[2, j+1], tb_ax)
                        Legend(fig[3, j+1], eess_ax)
                        hideydecorations!(eb_ax, grid=false)
                        hideydecorations!(tb_ax, grid=false)
                        hideydecorations!(eess_ax, grid=false)
                else
                        eb_ax.ylabel = L"$P$ [kW]"
                        tb_ax.ylabel = L"$P$ [kWt]"
                        eess_ax.ylabel = L"[%]"
                end
                limits!(eess_ax, t[1]/3600, t[end]/3600, 0, 100);
                linkxaxes!(eb_ax, eess_ax); linkxaxes!(eess_ax, tb_ax);
                hidexdecorations!.([eb_ax, tb_ax], grid=false, ticks=false)
                # limits for the inset axes
                limits!(eb_ax_inset, t[it0_zoom]/3600, t[itend_zoom]/3600, nothing, nothing);
                limits!(tb_ax_inset, t[it0_zoom]/3600, t[itend_zoom]/3600, nothing, nothing);
                limits!(eess_ax_inset, t[it0_zoom]/3600, t[itend_zoom]/3600, nothing, nothing);
                # move inset to front
                translate!.([eb_ax_inset.blockscene, tb_ax_inset.blockscene, eess_ax_inset.blockscene], 0, 0, 150)
                hideydecorations!.([eb_ax_inset, tb_ax_inset, eess_ax_inset], grid=false, ticks = false)
                if j == 1
                        textlabel!(eess_ax_inset, t[it0_zoom]/3600+13, 15, text= "Day $(day)", fontsize = 17)
                else
                        textlabel!(eess_ax_inset, t[it0_zoom]/3600+12, 15, text= "Day $(day)", fontsize = 17)
                end
        end
        linkyaxes!(axs[1,1], axs[1,2])
        linkyaxes!(axs[2,1], axs[2,2])
        linkyaxes!(axs[3,1], axs[3,2])
        textlabel!(axs[1,1], 530, -8, text= "Summer", fontsize = 18)
        textlabel!(axs[1,2], 530, -8, text= "Winter", fontsize = 18)
        textlabel!(axs[2,1], 260, 5.5, text= "Day 12", fontsize = 17)
        textlabel!(axs[2,2], 260, 5.5, text= "Day 23", fontsize = 17)
        if backend == "CairoMakie"
                save(filename, fig)
                return fig
        else
                return display(GLMakie.Screen(), fig)
        end
end

### Compilation Time plot for each DAtype and season
# lets plot all the compTime for each season and each DAtype in 2 subplots

# save the compTime for each DAtype in a dictionary
compTime = Dict(:summer => Dict(), :winter => Dict())
for season ∈ [:summer, :winter]
    for DAtype ∈ labels
        compTime[season][Symbol(DAtype)] = concatResultsFlex(results_DA[season][Symbol(DAtype)]["sim"], CT_MPC())["compTime"]
        println("The mean comptime of planner $DAtype in $(String(season)) is: ", mean(compTime[season][Symbol(DAtype)]))
    end
end

# plot the compTime for each DAtype, for each season use a different marker
GLMakie.activate!()
CairoMakie.activate!()
f=Figure(size=(800, 500))
# colors=ColorSchemes.tab10.colors;
colors=Makie.wong_colors();
set_theme!(theme_latexfonts(), fontsize = 20)
ax1=Axis(f[1,1];
        ylabel = L"$f(t_{\textrm{comp}})$",
        xscale = log10,
        yautolimitmargin = (0.1, 0.1)
        # yticks = (2:2:6, labels[1:3])
        )
ax2=Axis(f[2,1];
        xlabel=L"$log(t_{\textrm{comp}})$ [s]", ylabel = L"$f(t_{\textrm{comp}})$",
        xscale = log10,
        yautolimitmargin = (0.1, 0.1)
        # yticks = (2:2:6, labels[1:3])
        )
for (i, DAtype) ∈ enumerate(labels)
    logbins = 10 .^ range(-1., stop=1, length=30)
    # Makie.scatter!(ax1, 1:length(compTime[:summer][Symbol(DAtype)]), compTime[:summer][Symbol(DAtype)], 
    #         color=(colors[i], 0.5), strokewidth = 1, strokecolor = colors[i], label=DAtype * " " * "summer", marker=:circle, markersize=20, )
    # Makie.scatter!(ax1, 1:length(compTime[:winter][Symbol(DAtype)]), compTime[:winter][Symbol(DAtype)], 
    #         color=(colors[i], 0.5), strokewidth = 1, strokecolor = colors[i], label=DAtype * " " * "winter", marker=:diamond, markersize=20, )
    Makie.hist!(ax1, compTime[:summer][Symbol(DAtype)],
            offset = i*2, direction = :y, scale_to = 1, bins = logbins,
            color=(colors[i], 0.8), strokewidth = 1, strokecolor = colors[i])
    Makie.hist!(ax2, compTime[:winter][Symbol(DAtype)],
            offset = i*2, direction = :y, scale_to = 1, bins = logbins,
            color=(colors[i], 0.8), strokewidth = 1, strokecolor = colors[i], label = split(DAtype, "_")[1])
    # Makie.density!(ax1, compTime[:summer][Symbol(DAtype)],
    #         # offset = i*2,
    #         color=(colors[i], 0.8), strokewidth = 1, strokecolor = colors[i], label=DAtype,)
    # Makie.density!(ax2, compTime[:winter][Symbol(DAtype)],
    #         # offset = i*2,
    #         color=(colors[i], 0.8), strokewidth = 1, strokecolor = colors[i], label=DAtype,)
    mu_winter=round(mean(compTime[:winter][Symbol(DAtype)]),digits=2)
    mu_summer=round(mean(compTime[:summer][Symbol(DAtype)]),digits=2)
    textlabel!(ax2, 0.11, i*2+0.5, text = "μ=" * "$mu_winter", fontsize = 20)
    textlabel!(ax1, 0.11, i*2+0.5, text = "μ=" * "$mu_summer", fontsize = 20)
end
axislegend(ax2, position = :rt)
textlabel!(ax2, 5, 2.3, text = L"\textbf{Winter}", fontsize = 20)
textlabel!(ax1, 5, 2.3, text = L"\textbf{Summer}", fontsize = 20)
linkxaxes!(ax1, ax2); hidexdecorations!(ax1, grid = false); hideydecorations!.([ax1 ax2], grid = false, label = false)
# limits!(ax1, 10, 1400, 1.5, 7.5);
# limits!(ax2, 10, 1400, 1.5, 7.5)
f
folder = "../images/AppEnergy/cs1/"
save(folder * "CS1_compTime.pdf", f)
## Raincloud plot
folder = "../images/AppEnergy/cs1/"
plotCompT(compTime[:winter]; backend="GLMakie")
plotCompT(compTime; backend="CairoMakie", filename = folder * "CS1_compTime.pdf")

### BESS analysis
# lets plot all the power for each season and each DAtype in 2 subplots
label=L"$P_{\textrm{BESS}}$ [kW]"
P_key = "Pbess";
filename="CS1_Pbess_d29.pdf"
create_Pfig(results_DA, P_key, label; backend="GLMakie")
create_Pfig(results_DA, P_key, label; backend="CairoMakie", filename = folder*filename)
label = L"$SoC_{\textrm{BESS}}$ [p.u.]"
SoC_key = "SoCbess";
filename="CS1_SoCbess_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="GLMakie")
create_Pfig(results_DA, SoC_key, label; backend="CairoMakie", filename=folder*filename)

### Compare the power of the EV for each DAtype and season
# lets plot all the currents for each season and each DAtype in 2 subplots
# label = L"$i_{\textrm{EV}}$ [A]"
# i_key = "iev[1]"
# create_Pfig(results_DA, i_key, label)

# lets plot all the power for each season and each DAtype in 2 subplots
label = L"$P_{\textrm{EV, tot}}$ [kW]"
P_key = "PevTot[1]"; filename="CS1_PevTot_d29.pdf";
create_Pfig(results_DA, P_key, label; backend="GLMakie")
create_Pfig(results_DA, P_key, label; backend="CairoMakie", filename=folder*filename)

### Compare the SoC of the EV for each DAtype and season
# lets plot all the power for each season and each DAtype in 2 subplots
label = L"$SoC_{\textrm{EV}}$ [p.u.]"
SoC_key = "SoCev[1]"; filename="CS1_SoCev_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="GLMakie")
create_Pfig(results_DA, SoC_key, label; backend="CairoMakie", filename=folder*filename)

### Compare the power of the TESS for each DAtype and season
# lets plot all the power for each season and each DAtype in 2 subplots
label = L"$SoC_{\textrm{TESS}}$ [p.u.]"
SoC_key = "SoCtess"; filename="CS1_SoCtess_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="CairoMakie", filename=folder*filename)

### Compare ageing for each DAtype and season
label = L"$Q_{\textrm{BESS}}$ [p.u.]"
SoC_key = "Qbess"; filename="CS1_Qbess_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="GLMakie")
create_Pfig(results_DA, SoC_key, label; backend="CairoMakie", filename=folder*filename)
label = L"$Q_{\textrm{EV}}$ [p.u.]"
SoC_key = "Qev[1]"; filename="CS1_Qev_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="GLMakie")
create_Pfig(results_DA, SoC_key, label; backend="CairoMakie", filename=folder*filename)

### Compare controller state prediction and realized state
label = L"$SoC_{\textrm{BESS}}$ [p.u.]"
SoC_key = "SoCbess"; filename="CS1_SoCbess_prediction_error_d29.pdf";
compareSimCtrl(results_DA, SoC_key, label; backend="GLMakie")
# save the RMSE error between the prediction and the realization for each DAtype and season in a dataframe
using DataFrames
rmse_dict = Dict(); rej_dict = Dict();
for SoC_key ∈ ["SoCbess", "SoCev[1]", "Qbess", "Qev[1]"]
    rmse_dict[SoC_key] = DataFrame();
    rej_dict[SoC_key] = DataFrame();
    # RMSE error between the prediction and the realization
    rmse_df = DataFrame(DAtype = ["BNoDeg_W10", "CEmpDeg_W11", "CPBDeg_W11"], summer = zeros(3), winter = zeros(3))
    # Rejection rate of the 
    rej_df = DataFrame(DAtype = ["BNoDeg_W10", "CEmpDeg_W11", "CPBDeg_W11"], summer = zeros(3), winter = zeros(3))
    for season ∈ [:summer, :winter]
        for DAtype ∈ labels
            try
                concRes = concatResultsRH(results_DA[season][Symbol(DAtype)]["sim"]; typeOpt="day-ahead")
                concCtrl = concatResultsRH(results_DA[season][Symbol(DAtype)]["ctr"]; typeOpt="day-ahead")
                SoCsa_pos = convert.(Float64, concRes[SoC_key])
                SoCsa_pre = convert.(Float64, concCtrl[SoC_key])
                # RMSE error between the prediction and the realization
                rmse = sqrt(mean((SoCsa_pos .- SoCsa_pre).^2))
                rmse_df[rmse_df.DAtype .== DAtype, season] .= rmse
                # Rejection rate of the controller when SoC⁻ is at the lower bound (SoCmin) and SoC⁺ is not at the lower bound
                SoCmin = 0.15
                idx_sim_min = (SoCsa_pos .- SoCmin .≤ 1e-3)# idx in which SoC⁺ is in lower bound
                idx_rej = SoCsa_pre[idx_sim_min] .!= SoCsa_pos[idx_sim_min] # idx in which SoC⁻ is not in lower bound
                rej_df[rej_df.DAtype .== DAtype, season] .= typeof(idx_rej) == BitVector ? sum(idx_rej)/length(SoCsa_pos)*100 : nothing
            catch e
                println("Error occurred while processing $SoC_key for $DAtype in $season: $e")
            end
        end
    end
    rmse_dict[SoC_key] = rmse_df
    rej_dict[SoC_key] = rej_df
end

### Compare ageing for each DAtype and season
label = L"$Q_{\textrm{EV}}$ [p.u.]"
SoC_key = "Qev[1]"; filename="CS1_Qev_d29.pdf";
create_Pfig(results_DA, SoC_key, label; backend="GLMakie")

compareFEC(results_DA, ["iev[1]"], ["Qev[1]"];
        backend="GLMakie",step=200)
compareFEC(results_DA, ["ibess", "iev[1]"], ["Qbess", "Qev[1]"];
        backend="GLMakie",step=200)
compareFEC(results_DA, ["ibess", "iev[1]"], ["Qbess", "Qev[1]"];
        step = 250, backend="CairoMakie", filename=folder*"CS1_FECQloss_d29.pdf")
compareFEC_wLoss(results_DA, ["ibess", "iev[1]"], ["Qbess", "Qev[1]"];
        step = 250, backend="CairoMakie", filename=folder*"CS1_FECQloss_d29_CPBDeg.pdf")

# Compare eESS utilization with SoC vs P plots
folder = "../images/AppEnergy/cs1/"
kk = ["bess", "ev[1]", "tess"]
titles = [L"BESS", L"EV_1", L"TESS"]
lims = [(-7.0,7.0,0.1,1.0),
        (-12.5,12.5,0.1,1.0),
        (-6.0,6.0,0.1,1.0)]
countourSP(concatResultsFlex(results_DA[:CPBDeg_W11e4], CT_MPC()), kk, titles; 
    lims = lims)
create_SPfig(results_DA, "tess"; backend="CairoMakie",
             filename=folder*"CS1_SPtess_summerd29.pdf",
             cmaps=[:viridis, :turbid],
            #  cmaps=[:thermal, :acton],
             lims=([-6., 2.0], [0.3, 1.0]))
create_SPfig(results_DA, "bess"; backend="CairoMakie",
             filename=folder*"CS1_SPbessd29.pdf",
             cmaps=[:seaborn_icefire_gradient,:heat],
            #  cmaps=[:seaborn_rocket_gradient, :reds],
            #  cmaps=[:berlin, :managua],
             lims=([-10., 10.], [0.15, 1.0]))
create_SPfig(results_DA, "ev[1]"; backend="CairoMakie",
             filename=folder*"CS1_SPevd29.pdf",
            #  cmaps=[:seaborn_icefire_gradient,:heat],
             cmaps=[:delta,:heat],
             cutoff = 3.,
             lims=([-12.5, 12.5], [0.15, 1.]))

### Grid histogram
gridHist(results_DA, inputs)
gridHist(results_DA, inputs; backend="CairoMakie",
             filename=folder*"CS1_gridHist_d365.pdf")
# violin plots
function violinPlt(results_DA, key;
    backend::String="GLMakie",
    filename::String="fig.png",
    # cmaps::Symbol=:balance,
    # lims=([-7.5, 7.5], [0.2, 0.8]) # xmin xmax ymin ymax
    )
    @assert backend ∈ ["CairoMakie", "GLMakie"] "backend must be either CairoMakie or GLMakie"
    if backend == "CairoMakie"
        CairoMakie.activate!(type = "svg")
    else
        GLMakie.activate!()
    end

    set_theme!(theme_latexfonts())
    f=Figure(size=(700, 300))
    colors=ColorSchemes.Dark2_7.colors;
    axP = Axis(f[1,1]; ylabel=L"$P$ [kW]")
    axSoC = Axis(f[2,1]; ylabel=L"$SoC$ [p.u.]")
    # setup keys
    SoC_key = "SoC$key";
    key == "ev[1]" ? key = "evTot[1]" : nothing;
    P_key = "P$key";
    
    # loop over the season
    for (j,season) ∈ enumerate([:summer, :winter])
        # loop over the DAtype
        for (i, DAtype) ∈ enumerate(["BNoDeg_W10", "CEmpDeg_W11e4", "CPBDeg_W11e4", "CPBDeg_W10"])
            concRes = concatResultsRH(results_DA[season][Symbol(DAtype)]; typeOpt="day-ahead")
            Psa = convert.(Float64, concRes[P_key])
            SoCsa = convert.(Float64, concRes[SoC_key])
            if key == "evTot[1]"
                # exclude the points where the EV is not connected
                γ_cont = concRes[:"γ_cont"][:]
                Psa = Psa[γ_cont .== 1];
                SoCsa = SoCsa[γ_cont .== 1];
            end

            # Create a 2D histogram
            hist!(axP, Psa; bins=10, normalization= :pdf, scale_to=0.4*(-1)^j, label="$season",
                color=(colors[j], 0.6), offset=i, direction=:x, strokewidth = 1, strokecolor = colors[j])
            hist!(axSoC, SoCsa; bins=10, normalization= :pdf, scale_to=0.4*(-1)^j,
                color=(colors[j],0.6), offset=i, direction=:x, strokewidth = 1, strokecolor = colors[j])
        end
    end
    # change the xticks to the DAtype
    axSoC.xticks = ([1, 2, 3], ["BNoDeg_W10", "CEmpDeg_W11e4", "CPBDeg_W11e4"])
    # add legend
    axislegend(axP, position=:rt, merge=true, orientation = :horizontal)
    hidexdecorations!(axP, grid = false)
    if backend == "CairoMakie"
        save(filename, f)
        return f
    else
        return display(GLMakie.Screen(), f)
    end 
end

violinPlt(results_DA, "bess"; backend="CairoMakie",
             filename=relpath*"CS1_violinBESS.pdf")
violinPlt(results_DA, "ev[1]"; backend="CairoMakie",
             filename=relpath*"CS1_violinEV.pdf")

### Compare the cost for each DAtype and season
# costSummary = calcFullObj(concac_RHCPBDeg, EMSData, s)
costSummary_DA = Dict{Symbol, Dict}(:summer=>Dict(), :winter=>Dict())
dailyCgrid = Dict{Symbol, Dict}(:summer=>Dict(), :winter=>Dict())
# costSummary_DA = Dict(); 
Wgrid = 1;
Wloss = 1;
# Wloss = 1;
W=[Wgrid 1000 Wloss 0];

for season ∈ [:summer, :winter]
    for DAtype ∈ labels
        concRes = concatResultsFlex(results_DA[season][Symbol(DAtype)]["sim"], CT_MPC())
        t0 = concRes["t"][1]; 
        Δt = concRes["t"][2] - concRes["t"][1];
        s=modelSettings(nEV=1, t0=t0,Tw=48-1/4, Δt=Δt, steps=29, costWeights=W,
                season=String(season), profType="monthly", loadType="GV", year=2023);
        costSummary_DA[season][Symbol(DAtype)] = calcFullObj(concRes, inputs[season], s);
        # dailyCgrid[season][Symbol(DAtype)] = calcDailyCost(concRes, inputs[season], s);
    end
end

s = modelSettings(nEV=1, t0=0,Tw=48-1/4, Δt=900, steps=29, costWeights=W,
    season="winter", profType="yearly", loadType="GV", year=2023);

dailyCgridBoxplot(dailyCgrid[:winter], s; backend="GLMakie", labels=labels)
dailyCgridBoxplot(dailyCgrid, s; backend="CairoMakie",
    filename=folder*"CS1_dailyCgridBoxplot_d365.pdf", labels=labels)

dailyCgridPlot(dailyCgrid; backend="GLMakie", filename=folder*"CS1_dailyCgridPlot_d29.pdf")

# for each day create sort who is the most expensive DAtype
mostExpensive = [findmax([dailyCgrid[:winter][Symbol(DAtype)][day] for DAtype ∈ labels]) for day ∈ 1:29]
mostExpensive_labels = [labels[mostExpensive[i][2]] for i ∈ 1:29]

# save in a DataFrame
dailyCgrid_DF = DataFrame(day=1:s.steps,
                        month = month.(DateTime.(Date(2023,1,1):Day(1):Date(2023,12,31))),
                        BNoDeg_W10 = [dailyCgrid[:BNoDeg_W10][day] for day ∈ 1:s.steps],
                        CEmpDeg_W11 = [dailyCgrid[:CEmpDeg_W11][day] for day ∈ 1:s.steps],
                        CPBDeg_W11e4 = [dailyCgrid[:CPBDeg_W11e4][day] for day ∈ 1:s.steps],
                        CPBDeg_W11e4_wst = [dailyCgrid[:CPBDeg_W11e4_wst][day] for day ∈ 1:s.steps],
                        MostExpCtr = mostExpensive_labels,
                        value1 = [mostExpensive[i][1] for i ∈ 1:s.steps],
                        LessExpCtr = lessExpensive_labels,
                        value2 = [lessExpensive[i][1] for i ∈ 1:s.steps])
# make a heatmap of the dailyCgrid for each DAtype
GLMakie.activate!()
f = Figure(size=(800, 400)); ax = Axis(f[1,1], xlabel="day", yticks=(1:length(labels),labels))
# take the mean per month
dailyCgrid_DF_month = combine(groupby(dailyCgrid_DF, :month), :BNoDeg_W10 => mean, :CEmpDeg_W11 => mean,
    :CPBDeg_W11e4 => mean, :CPBDeg_W11e4_wst => mean)
dCg_m = Matrix(dailyCgrid_DF_month[!,2:end])
dCg_m = dCg_m .- mean(dCg_m, dims=2)
dCg_m = dCg_m ./ maximum(dCg_m)
heatmap!(ax, 1:12, 1:length(labels), dCg_m, colormap = :balance,
    )
# normalize and center
# dCg = dCg .- mean(dCg, dims=2)
# dCg = dCg ./ maximum(dCg)
# heatmap!(ax, 1:s.steps, 1:length(labels), dCg, colormap = :balance, colorrange = (0,1),
#     )
Colorbar(f[1,2],  colormap = :balance)
f
plotname ="DACPBDeg_W11_y2023d29_EMSplots_detail.pdf"

makeEMSplots(concatResultsRH(results_DA[:summer][:CPBDeg_W11][:"sim"], typeOpt ="day-ahead"),
            inputs[:summer]; backend="GLMakie", plot_ageing=false)
makeEMSplots(rhDict, inputs; backend="CairoMakie", filename= relpath * plotname, plot_ageing=false)
folder = "../images/AppEnergy/cs1/"
plotDispatchAnalysis(results_DA, inputs; backend="CairoMakie", filename=folder*"CS1_dispatchAnalysis.pdf")
makeInputsplot([inputs[:summer], inputs[:winter]]; backend="CairoMakie", filename=folder*"inputs.pdf")

costSummaryPlot(results_DA, costSummary_DA; backend="GLMakie")
costSummaryPlot(results_DA, costSummary_DA; backend="CairoMakie",
    filename=folder*"CS1_costSummary.pdf")
# make a DataFrame with the costSummary for each DAtype and season summarizing the grid cost
costSummary_DF = DataFrame();

for season ∈ [:summer, :winter]
    for DAtype ∈ keys(results_DA[season])
        costSummary_DF = vcat(costSummary_DF, DataFrame(DAtype=DAtype, season=season,
            wCgrid=costSummary_DA[season][Symbol(DAtype)].wCgrid[end],
            Qloss=costSummary_DA[season][Symbol(DAtype)].wCloss[end]./Wloss*1e3/1.2))
    end
end

for DAtype ∈ keys(results_DA)
    costSummary_DF = vcat(costSummary_DF, DataFrame(DAtype=DAtype,
        wCgrid = costSummary_DA[Symbol(DAtype)].wCgrid[end],
        wCloss = costSummary_DA[Symbol(DAtype)].wCloss[end],
        pSoCDep = costSummary_DA[Symbol(DAtype)].wpDep,
        ))
end

# Rejections for each DAtype
SoCmin = 0.2
rej_DF = DataFrame()
for (j, season) ∈ enumerate(keys(results_DA))
    for (da, DA) ∈ enumerate(keys(results_DA[season]))
        for (d, key) ∈ enumerate(["SoCbess", "SoCev[1]"])
            sim = concatResultsRH(results_DA[season][DA]["sim"]; typeOpt = "day-ahead")
            ctrl = concatResultsRH(results_DA[season][DA]["ctr"]; typeOpt = "day-ahead")
            idx_sim_min = (sim[key] .- SoCmin .≤ 1e-4)# idx in which SoC⁺ is in lower bound
            idx_rej = ctrl[key][idx_sim_min] .!= sim[key][idx_sim_min] # idx in which SoC⁻ is not in lower bound
            println("The % of rejection (SoC⁻ ≠ SoC⁺ == SoCmin) for $(String(season)) $(String(DA)) $(String(key)) is:")
            println(typeof(idx_rej) == BitVector ? sum(idx_rej)/length(sim[key])*100 : nothing)
            # save results in a dataframe
            if d == 1 && da == 1 && j == 1
                rej_DF = DataFrame(DAtype=DA, season=season, key=key, rej=sum(idx_rej)/length(sim[key])*100)
            else
                typeof(idx_rej) == BitVector ? push!(rej_DF, (DA, season, key, sum(idx_rej)/length(sim[key])*100)) : nothing;
            end
        end
    end
end

# Sensitivity analysis
CairoMakie.activate!(type = "svg")
set_theme!(theme_latexfonts(), fontsize = 20)
f = Figure(size = (700, 700/1.618))
ax = Axis(f[1,1], xlabel = L"C_{\text{grid}}", ylabel = L"Q_{\text{loss}}",
    yautolimitmargin = (0.1, 0.1), xautolimitmargin = (0.1, 0.2))
# Define mapping for consistency
seasons = [:summer, :winter]
colors = Makie.wong_colors()
color_map = Dict(:summer => colors[2], :winter => colors[1])

for s in seasons
    # Filter once per loop iteration
    df_sub = subset(costSummary_DF, :season => x -> x .== s) 
    # Plot points
    scatter!(ax, df_sub.wCgrid, df_sub.Qloss,
        markersize = 20, 
        strokewidth = 2.0, 
        strokecolor = color_map[s], 
        color = (color_map[s], 0.5), 
        label = uppercasefirst(string(s)))

    # Add labels
    # textlabel!(ax, df_sub.wCgrid, df_sub.Qloss, 
    #     text = string.(df_sub.DAtype), 
    #     text_align = (:left, :bottom), 
    #     offset = (20, -15), 
    #     fontsize = 13)
end
axislegend(ax, position = :rt)
f
folder = "../images/AppEnergy/cs1/"
save(folder * "CS1_Cgrid_v_Qloss.svg", f)