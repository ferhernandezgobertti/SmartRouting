# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    SpectralApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of Spectral clustering algorithms through explicit 
procedural approaches.

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

using CSV, DataFrames, Dates, Stop
using Clustering, StatsBase, NamedArrays, LinearAlgebra

struct SpectralApproach <: AbstractSegmentationAlgorithm end

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

function get_combinations_l1(route)
    L = length(route.stops)
    distances = Matrix{Float64}(undef, L, L)
    for i in 1:L
        stop_init = route.stops[i]
        for j in 1:L
            stop_final = route.stops[j]
            distances[i][j] = l1_distance(stop_init, stop_final)
        end
    end
    return distances
end

function get_combinations_l2(route)
    L = length(route.stops)
    distances = Matrix{Float64}(undef, L, L)
    for i in 1:L
        stop_init = route.stops[i]
        for j in 1:L
            stop_final = route.stops[j]
            distances[i][j] = l2_distance(stop_init, stop_final)
        end
    end
    return distances
end

function get_stops_cost_matrices(R, routedistance, traveltimes, velc, vela, veld, dParams)
    stops_cost = Matrix{Float64}(undef, size(traveltimes)[1], size(traveltimes)[1])
    lambda = [0.005, 0.007, 0.004, 1.23, 1.55, 1.69, 0.02, 0.01, 0.04]
    if(dParams[1] == 1)
        stops_cost = lambda[2]*traveltimes+lambda[4]*routedistance
    end
    if(dParams[1] == 2)
        mandistance = get_combinations_l1(R)
        eucdistance = get_combinations_l2(R)
        stops_cost = lambda[2]*traveltimes+lambda[5]*eucdistance+lambda[6]*mandistance
    end
    if(dParams[1] == 3)
        mandistance = get_combinations_l1(R)
        eucdistance = get_combinations_l2(R)
        stops_cost = lambda[4]*routedistance+lambda[5]*eucdistance+lambda[6]*mandistance
    end
    if(dParams[1] == 4)
        eucdistance = get_combinations_l2(R)
        stops_cost = lambda[1]*traveltimes+lambda[5]*eucdistance+lambda[7]*velc+lambda[8]*vela+lambda[9]*veld
    end
    if(dParams[1] == 5)
        mandistance = get_combinations_l1(R)
        eucdistance = get_combinations_l2(R)
        stops_cost = lambda[1]*traveltimes+lambda[5]*eucdistance+lambda[6]*mandistance+lambda[8]*vela+lambda[9]*veld
    end
    return stops_cost
end

function get_stops_coords(route)
    stops_cost = Matrix{Float64}(undef, 2, length(unique(route.stops)))
    stops_count = 1
    for S in unique(route.stops)
        s_lat = latitude(S)
        s_lng = longitude(S)
        stops_cost[1, stops_count] = s_lng
        stops_cost[2, stops_count] = s_lat
        stops_count = stops_count + 1
    end
    return stops_cost
end

function get_adjacency(stops_cost, dparams)
    adjacency_cost = Matrix{Float64}(undef, size(stops_cost)[1], size(stops_cost)[2])
    for i in 1:size(stops_cost,1)
        for j in 1:size(stops_cost,2)
            if(stops_cost[i,j] <= dparams[6])
                adjacency_cost[i,j] = 1
            else
                adjacency_cost[i,j] = 0
            end
        end
    end
    return adjacency_cost
end

function get_degrees(stops_cost)
    degrees_cost = Matrix{Float64}(undef, size(stops_cost)[1], size(stops_cost)[2])
    for i in 1:size(stops_cost,1)
        for j in 1:size(stops_cost,2)
            if(i=j)
                degrees_cost[i,j] = count(x->x!=0, stops_cost[i,:])
            else
                degrees_cost[i,j] = 0
            end
        end
    end
    return degrees_cost
end

function spectral_procedural(stops_cost, dparams)
    num_stops = size(stops_cost,1)
    degree_cost = get_degrees(stops_cost)

    # Laplacian calculation
    if(dparams[9] == 1)
        lap_cost = degree_cost - stops_cost
    else
        lap_cost = transpose(degree_cost)*stops_cost*degree_cost
    end

    # Eigenvalues and Eigenvectors of Laplacian
    eig_vals = eigvals(lap_cost)
    eig_vecs = eigvecs(lap_cost)

    # Sort based on Eigenvalues
    eig_vecs_sorted = eig_vecs[sortperm(eig_vals)]

    # Centroid-based Clusterings on Eigenvectors
    if(dparams[10] == 1)
        clust_result = KMeans(eig_vecs_sorted[:,1:dparams[2]], dparams[2]; maxiter=dparams[3], tol=dparams[4]) # KMeans on Eigenvectors
    else
        clust_result = KMedoids(eig_vecs_sorted[:,1:dparams[2]], dparams[2]; maxiter=dparams[3], tol=dparams[4]) #KMedoids on Eigenvectors
    end
    
    clusters = assignments(clust_result)
    cluster_sizes = [length(findall(x->x==i,clusters)) for i in 1:length(clusters)]
    return clusters, cluster_sizes
end

function solve(alg::SpectralApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3, 4, 5}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #process_no = dParams[5] # {1, 2} = {Procedural Adjacency, Procedural Affinity}
    #max_distance = dParams[6] # {1.21, 2.64, 3.82}
    #main_matrix_scale = dParams[7] # {1, 5, 9}
    #degree_matrix_scale = dParams[8] # {1, 2, 4}
    #laplacian_type = dParams[9] # {1, 2} = {Differential (D-M), Expected Value (DT)(M)(D)}
    #centroid_clust = dParams[10] # {1, 2} = {KMeans, KMedoids}
    dParams = [2, 21, 20, 1, 2, 2, 5, 3, 1]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::SpectralApproach, routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    
    for R in routes
        #stops_cost = get_stops_coords(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)

        if(dParams[5] == 1)
            adjacency_stops = get_adjacency(stops_cost, params)
            ass_stops, cl_sizes = spectral_procedural(adjacency_stops, params)
        else
            ass_stops, cl_sizes = spectral_procedural(stops_cost, params)
        end
        push!(zones_quant_arr, length(cl_centers))
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