# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    GaussianApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of K-Means algorithms through explicit procedural and 
declarative approaches, as well as with optimized libraries and JuMP solvers.

Inputs: 
* routes: List of Transports Routes, i.e. Vector of Vector of Stops (Lat, Lon).
* traveltimes: Matrix associated with times for travelling from one stop to
               another (all combinations), such as expected tt or static tt.

[Optional]
* params: Relevant parameters to configure Segmentation process.

Outputs:
* routes: List of updated routes (i.e. stops with zones generic identifiers).
* dzones: List of zones (i.e. Segmented Stops)
* zones_quant_arr: List of Zones Ammount per route used on stops.

Note: In this step, all zone generic identifiers are direct (e.g. 1, 2, 3, ...)
"""

using CSV, DataFrames, Dates
using Clustering, StatsBase, NamedArrays, GaussianMixtures
using JuMP, GLPK, CPLEX, Gurobi

# Parallelization Opportunity:
# using ClusterManagers
# ClusterManagers.addprocs_sge(20)                                        
# @everywhere using GaussianMixtures

struct GaussianApproach <: AbstractSegmentationAlgorithm end

struct DZones
    ids::Vector{String}
    stops::Vector{Stop}
end

function get_seg_score(segments)
    seg_score = 0
    seg_zones = Vector{String}()
    seg_cant_zones = Vector{Int64}()
    for S in segments
        current_zone = zone(S)
        if(!(current_zone in seg_zones))
            push!(seg_zones, zone(S))
            push!(seg_cant_zones, 1)
        else
            zone_pos = findfirst(x->x==zone(S), seg_zones)
            seg_cant_zones[zone_pos] += 1
        end
    end
    return maximum(seg_cant_zones)/sum(seg_cant_zones)
end

function compare_zones_seg(predicted_zones)
    current_score = 0
    for PZ in predicted_zones
        zone_score = get_seg_score(PZ)
        current_score = current_score + zone_score
    end
    return current_score/length(predicted_zones)
end

function get_stops_coords(route)
    stops_mat = Matrix{Float64}(undef, 2, length(unique(route.stops)))
    stops_count = 1
    for S in unique(route.stops)
        s_lat = latitude(S)
        s_lng = longitude(S)
        stops_mat[1, stops_count] = s_lng
        stops_mat[2, stops_count] = s_lat
        stops_count = stops_count + 1
    end
    return stops_mat
end

function gaussian_library(stops_cost, dparams)
    #clusters = kmeans(stops_mat, dparams[2]; maxiter=dparams[3], tol=dparams[4])

    # Definition of Gaussian Model
    if(dparams[7] == 1)
        if(dparams[8] == 1)
            gauss_model = GMM(dparams[5], stops_cost; method=:kmeans, kind=:diag, nInit=dparams[6], nIter=dparams[3])
        else
            gauss_model = GMM(dparams[5], stops_cost; method=:kmeans, kind=:full, nInit=dparams[6], nIter=dparams[3])
        end
    else
        if(dparams[8] == 1)
            gauss_model = GMM(dparams[5], stops_cost; method=:split, kind=:diag, nInit=dparams[6], nIter=dparams[3])
        else
            gauss_model = GMM(dparams[5], stops_cost; method=:split, kind=:full, nInit=dparams[6], nIter=dparams[3])
        end 
    end

    # Application of Gaussian Model
    if(dparams[9] == 1)
        prob_clust = em!(gauss_model, stops_cost; nIter=dparams[3])[1]
    elseif(dparams[9] == 2)
        prob_clust = llpg(gauss_model, stops_cost)[1]
    elseif(dparams[9] == 3)
        prob_clust = avll(gauss_model, stops_cost)[1]
    else
        prob_clust = gmmposterior(gauss_model, stops_cost)[1]
    end
    clusters = [indmax(prob_clust[i,:]) for i=1:size(stops_cost,1)] # Tolerance dparams[4] can be added
    cluster_sizes = [length(findall(x->x==i,clusters)) for i in 1:length(clusters)]
    return clusters, cluster_sizes, clust_centers
end

function solve(alg::GaussianApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #mixture_quantity = dParams[5] # {3, 7, 10, dparams[2]}
    #max_iter_em = dParams[6] # {10, 20, 30}
    #clust_method = dParams[7] # {1, 2} = {KMeans, Split}
    #covar_matrix = dParams[8] # {1, 2} = {Diag, Full}
    #train_func = dParams[9] # {1, 2, 3, 4} = {Expectation Maximization, Log Likelihood, Average Log Likelihood, Posterior}
    dParams = [2, 21, 20, 1, 3, 10, 1]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::GaussianApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    solver = JuMP.Gurobi
    
    for R in routes
        stops_cost = get_stops_cost_matrices(R, traveltimes, params)
        #stops_mat = get_stops_coords(R)
        ass_stops, cl_sizes = gaussian_library(stops_cost, params)

        stops = unique(R.stops, 1)
        for S in stops
            curr_ass = string(ass_stops[findfirst(x->x.id==S.id, stops)])
            S.zone = Zone("", curr_ass)
            #prev_zone = findfirst(x->x.ids==curr_ass, zones)
            #if(isnothing(prev_zone))
            push!(zones.ids, curr_ass)
            push!(zones.stops, S)
        end
    end
    return routes, dZones, zones_quant_arr
end