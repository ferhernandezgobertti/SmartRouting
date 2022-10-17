# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    KMeansApproach <: AbstractSegmentationAlgorithm

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
using Clustering, StatsBase, NamedArrays
using JuMP, GLPK, CPLEX, Gurobi

struct KMeansApproach <: AbstractSegmentationAlgorithm end

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

function kmeans_procedural(stops_mat, dparams)
    X = transpose(stops_mat)
    cluster_no = dparams[2]
    max_iter = dparams[3]
    threshold = dparams[4]
    centroids = X[:, sample(1:size(X,2), cluster_no, replace = false)]
    new_centroids = copy(centroids)
    cluster_ass = zeros(Float64, size(X, 2))

    for _ in 1:max_iter             # For each iteration until iteration_quantity
        for col_idx in 1:size(X, 2) # Iteration over each point
            p = X[:, col_idx]
            point_difference = mapslices(x -> x .- p, centroids, dims=[1])
            distances = mapslices(sum, point_difference .^ 2, dims=[1])
            cluster_ids[col_idx] = findmin(distances)[2][2]
            # println("p: $p diff: $point_difference dist: $distances $cluster_ids")
        end
        # println("old: $centroids new: $new_centroids")

        for cluster_id in 1:size(centroids, 2) # Iteration over each centroid
            mask = [i for (i, m) in enumerate(cluster_id .== cluster_ids) if m]
            new_centroids[:, cluster_id] = mapslices(mean, X[:, mask], dims=[2])
        end
        # println("old_centroids: $centroids new_centroids: $new_centroids point assignemnts: $cluster_ids")

        center_change = sum(mapslices(x -> sum(x.^2), new_centroids .- centroids, dims=[2]))
        centroids = copy(new_centroids)
        if center_change < threshold
            break
        end
    end

    cluster_sizes = [length(findall(x->x==i,cluster_ids)) for i in 1:length(centroids)]

    return cluster_ids, cluster_sizes, centroids 
end

function kmeans_library(stops_mat, dparams)
    clusters = kmeans(stops_mat, dparams[2]; maxiter=dparams[3], tol=dparams[4])
    #@assert nclusters(R1) == dparams[2] # verify the number of clusters
    assignment_stops = assignments(clusters) # get the assignments of points to clusters
    clust_sizes = counts(clusters) # get the cluster sizes
    clust_centers = clusters.centers # get the cluster centers
    return assignment_stops, clust_sizes, clust_centers
end

function apply_kmeans_vmeasure(data_route, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
        data_zones_mod = Vector{Vector{Vector{Stop}}}()
        stops_mat = Matrix{Float64}(undef, 2, length(data_route[route_num].stops))
        stops_count = 1
        for S in data_route[route_num].stops
            s_lat = latitude(S)
            s_lng = longitude(S)
            stops_mat[1, stops_count] = s_lng
            stops_mat[2, stops_count] = s_lat
            stops_count = stops_count + 1
        end

        R = kmeans(stops_mat, zones_cant; maxiter=200, tol=tolerance) # NO DIF: length(zones(data_route[route_num])); maxiter=200, tol=tolerance)  

        a = assignments(R) # get the assignments of points to clusters
        c = counts(R) # get the cluster sizes
        M = R.centers # get the cluster centers

        stops_zones = Vector{Vector{Stop}}(undef, zones_cant) # NO DIF: length(zones(data_route[route_num]))) 
        
        vec_compare = Vector{Int64}()
        zones_route = zones(data_route[route_num])
        for S in data_route[route_num].stops
            zone_pos = findfirst(x->x==zone(S), zones_route)
            push!(vec_compare, zone_pos)
        end
        
        score_reached = vmeasure(R, vec_compare)
        
        push!(scores_arr, score_reached)
        #println("SCORE of ROUTE "*string(route_num)*": "*string(score_reached))
    end
    return scores_arr
end

function apply_kmeans(data_route, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
        data_zones_mod = Vector{Vector{Vector{Stop}}}()
        stops_mat = Matrix{Float64}(undef, 2, length(data_route[route_num].stops))
        stops_count = 1
        for S in data_route[route_num].stops
            s_lat = latitude(S)
            s_lng = longitude(S)
            stops_mat[1, stops_count] = s_lng
            stops_mat[2, stops_count] = s_lat
            stops_count = stops_count + 1
        end

        R = kmeans(stops_mat, zones_cant; maxiter=200, tol=tolerance) # NO DIF: length(zones(data_route[route_num])); maxiter=200, tol=tolerance)

        a = assignments(R) # get the assignments of points to clusters
        c = counts(R) # get the cluster sizes
        M = R.centers # get the cluster centers

        stops_zones = Vector{Vector{Stop}}(undef, zones_cant) # NO DIF: length(zones(data_route[route_num]))) 
        cluster_unique = unique(a)

        for i in 1:length(a)
            stop_zone = Vector{Stop}()
            cluster_pos = a[i]
            current_stop = data_route[route_num].stops[i]
            current_pos = findfirst(x->x==cluster_pos,cluster_unique)
            if(!isassigned(stops_zones, cluster_pos)) #cluster_unique[current_pos]))
                push!(stop_zone, current_stop)
                insert!(stops_zones, cluster_pos, stop_zone)
                deleteat!(stops_zones, length(stops_zones))
            else
                push!(stops_zones[cluster_pos], current_stop)    
            end
        end
        stops_zones = [stops_zones[i] for i in 1:length(stops_zones) if isassigned(stops_zones, i)]
        score_reached = compare_zones_seg(stops_zones)
        #new_score = compare_zones_seg(stops_zones)
        #println("NEW SCORE: "*string(new_score))
        #println("CANT ZONES ROUTE "*string(route_num)*": "*string(length(zones(data_route_apply[route_num]))))
        #score_reached = new_score/length(zones(data_route_apply[route_num]))
        push!(scores_arr, score_reached)
        #println("SCORE of ROUTE "*string(route_num)*": "*string(score_reached))
    end
    return scores_arr
end

function preparation_kmeans(data_route)
    iter_kmeans = 10
    zones_cant = 21
    tolerance = 1
    mean_score = 0
    for iter in 1:iter_kmeans
        scores_arr = apply_kmeans(data_route, zones_cant, tolerance)
        #println("SCORES: "*string(scores_arr))
        mean_score = mean_score + mean(scores_arr)
        println("MEAN SCORE: "*string(mean(scores_arr)))
    end
    println("MEAN SCORE TOTAL: "*string(mean_score/iter_kmeans))
end

function kmeans_jump(stops_mat, dparams, opt_solver=Gurobi)
    
    # JuMP PreProcessing
    n = 2
    threshold = dparams[4]
    m = length(stops_mat[1])
    k = dparams[2]
    W = zeros(m, m)
    for i in 1:m
        for j in i+1:m
            W[i, j] = W[j, i] = exp(-LinearAlgebra.norm(stops_mat[i] - stops_mat[j]) / 1.0)
        end
    end
    
    # JuMP Modelling
    model = Model(opt_solver.Optimizer)
    #set_silent(model)
    @variable(model, Z[1:m, 1:m], PSD)
    fix(Z[1, 1], 1; force = true) # reducing symetry
    @constraint(model, Z .>= 0) # Z >= 0, PSD
    I = Matrix(1.0 * LinearAlgebra.I, m, m)
    @objective(model, Min, LinearAlgebra.tr(W * (I - Z))) # min Tr(W(I-Z))
    @constraint(model, Z * ones(m) .== ones(m)) # Z e = e
    @constraint(model, LinearAlgebra.tr(Z) == k) # Tr(Z) = k
    optimize!(model)
    Z_val = value.(Z)
    
    # JuMP PostProcessing
    ass_cluster = zeros(Int, m)
    num_clusters = 0
    for i in 1:m
        if Z_val[i, i] <= threshold
            continue
        elseif ass_cluster[i] == 0
            num_clusters += 1
            ass_cluster[i] = num_clusters
            for j in i+1:m
                if LinearAlgebra.norm(Z_val[i, j] - Z_val[i, i]) <= threshold
                    ass_cluster[j] = num_clusters
                end
            end
        end
    end
    cluster_sizes = [length(findall(x->x==i,ass_cluster)) for i in 1:length(ass_cluster)]
    return ass_cluster, cluster_sizes, _
end

function solve(alg::KMeansApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #procedure_no = dParams[5] # {1, 2, 3} = {Procedural, Library, JuMP}
    dParams = [2, 21, 20, 1, 2]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::KMeansApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    solver = JuMP.Gurobi
    
    for R in routes
        #stops_cost = get_stops_cost_matrices(R, traveltimes, params)
        stops_mat = get_stops_coords(R)

        if(dParams[5] == 1)
            ass_stops, cl_sizes, cl_centers = kmeans_procedural(stops_mat, params)
        end
        if(dParams[5] = 2)
            ass_stops, cl_sizes, cl_centers = kmeans_library(stops_mat, params)
        end
        if(dParams[5] = 3)
            ass_stops, cl_sizes, cl_centers = kmeans_jump(stops_mat, params, opt_solver=solver)
        end

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