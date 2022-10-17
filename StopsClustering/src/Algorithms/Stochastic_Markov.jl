# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    MarkovApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of Markov-based algorithms through explicit procedural and 
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

using CSV, DataFrames, Dates, Stop
using Clustering, StatsBase, NamedArrays
using JuMP, GLPK, CPLEX, Gurobi

struct MarkovApproach <: AbstractSegmentationAlgorithm end

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

function process_final_matrix(stops_cost, zero_tolerance)
    # removal of rows containing only zero elements and convert into a mask of nonzero elements
    clusters_mask = stops_cost[dropdims(sum(stops_cost, dims=2), dims=2) .> zero_tolerance, :] .> zero_tolerance

    # associate cluster indexes to each node
    # cluster index is the index of the first TRUE in a given column
    _ms = mapslices(some_mask->isempty(some_mask) ? 0 : argmax(some_mask), clusters_mask, dims=1)
    clusters = dropdims(_ms, dims=1)
    cluster_sizes = zeros(Int, size(clusters_mask, 1))
    @inbounds for clu_ix in clusters
        (clu_ix > 0) && (cluster_sizes[clu_ix] += 1)
    end
    @inbounds for i in eachindex(clusters)
        clu_ix = clusters[i]
        if clu_ix > 0 && cluster_sizes[clu_ix] == 1
            clusters[i] = 0
            cluster_sizes[clu_ix] = 0
        end
    end

    # change clusters numbers to be in 1:N range (or 0:N if there's collapsed cluster)
    clu_id_map = zeros(Int, length(cluster_sizes))
    next_clu_ix = 0
    @inbounds for i in eachindex(clusters)
        old_clu_ix = clusters[i]
        if old_clu_ix > 0
            new_clu_ix = clu_id_map[old_clu_ix]
            clusters[i] = new_clu_ix == 0 ? clu_id_map[old_clu_ix] = (next_clu_ix += 1) : new_clu_ix
        end
    end
    old_cluster_sizes = cluster_sizes
    cluster_sizes = zeros(Int, next_clu_ix)
    for (old_clu_ix, new_clu_ix) in enumerate(clu_id_map)
        if new_clu_ix > 0
            cluster_sizes[new_clu_ix] = old_cluster_sizes[old_clu_ix]
        end
    end

    return clusters, cluster_sizes
end

function markov_procedural(stops_cost, dparams)

    niter = 0
    converged = false
    rel_delta = NaN
    while !converged && niter < dparams[3]
        expanded = stops_cost^dparams[7]

        @inbounds for i in eachindex(expanded)
            next_stops_cost[i] = expanded[i]^dparams[8]
        end

        for i in 1:size(next_stops_cost,2)
            c = view(next_stops_cost, :, i)
            mean_prune = mean(c)*dparams[9]
            @inbounds for j in eachindex(c)
                c[j] = ifelse(c[j] >= mean_prune, c[j], 0.0)
            end
        end
        if(issparse(next_stops_cost))
            dropzeros!(next_stops_cost)
        end

        # normalization
        rmul!(next_stops_cost, Diagonal(map(x -> x != 0.0 ? 1.0/x : x, dropdims(sum(next_stops_cost, dims=1), dims=1))))
        next_mcl_norm = norm(next_stops_cost)

        if !isfinite(next_mcl_norm)
            break
        end
        rel_delta = euclidean(next_stops_cost, stops_cost)/mcl_norm
        if(converged = rel_delta <= dparams[4])
            break
        end
        # update (swap) MCL adjacency
        niter += 1
        stops_cost, next_stops_cost = next_stops_cost, stops_cost
        mcl_norm = next_mcl_norm
        if(mcl_norm < dparams[4])
            break # matrix is zero
        end
    end
    clusters, cluster_sizes = process_final_matrix(stops_cost, dparams[4]/length(stops_cost[1]))

    return clusters, cluster_sizes 
end

function markov_library(stops_cost, dparams)
    clusters = mcl(stops_cost; add_loops=true, expansion=dparams[7], inflation=dparams[8], save_final_matrix=false, prune_tol=dparams[9], maxiter=dparams[3], tol=dparams[4])
    #@assert nclusters(R1) == dparams[2] # verify the number of clusters
    assignment_stops = assignments(clusters) # get the assignments of points to clusters
    clust_sizes = counts(clusters) # get the cluster sizes
    return assignment_stops, clust_sizes
end

function apply_markov_vmeasure(data_route, data_route_act, tt_data, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
    
        TT = tt_data[route_num].travel_times
        #TT = TT[Not(name(station(data_route_act[route_num]))),Not(name(station(data_route_act[route_num])))]
        stops_0 = data_route_act[route_num].stops[findall(x->zone(x)=="0", data_route_act[route_num].stops)]
        for S in stops_0
            TT = TT[Not(name(S)),Not(name(S))]
        end
        
        R = mcl(TT; add_loops=true, expansion=2, inflation=1, save_final_matrix=false, prune_tol=1, maxiter=200, tol=tolerance)
        a = assignments(R) # get the assignments of points to clusters
        #c = counts(R) # get the cluster sizes
        #M = R.centers # get the cluster centers

        stops_zones = Vector{Vector{Stop}}(undef, zones_cant)
        cluster_unique = unique(a)
        stops_zones = Vector{Vector{Stop}}(undef, zones_cant) # NO DIF: length(zones(data_route[route_num]))) 
        
        vec_compare = Vector{Int64}()
        zones_route = zones(data_route[route_num])
        for S in data_route[route_num].stops
            zone_pos = findfirst(x->x==zone(S), zones_route)
            push!(vec_compare, zone_pos)
        end
        
        #println("DIMENSION A: "*string(length(a)))
        #println("DIMENSION VEC: "*string(length(vec_compare)))
        score_reached = vmeasure(a, vec_compare)
        push!(scores_arr, score_reached)
        println("SCORE of ROUTE "*string(route_num)*": "*string(score_reached))
    end
    return scores_arr
end

function apply_markov(data_route, data_route_act, tt_data, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
    
        TT = tt_data[route_num].travel_times
        TT = TT[Not(name(station(data_route_act[route_num]))),Not(name(station(data_route_act[route_num])))]
        
        R = mcl(TT; add_loops=true, expansion=2, inflation=1, save_final_matrix=false, prune_tol=1, maxiter=200, tol=tolerance)
        a = assignments(R) # get the assignments of points to clusters
        #c = counts(R) # get the cluster sizes
        #M = R.centers # get the cluster centers

        stops_zones = Vector{Vector{Stop}}(undef, zones_cant)
        cluster_unique = unique(a)
        #println("LENGTH STOPS: "*string(length(data_route[route_num].stops)))
        #println("LENGTH A: "*string(length(a)))
        for i in 1:length(data_route[route_num].stops)
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
    
function preparation_markov(data_route, data_route_act)
    #data_route = data_mod_apply #data_mod_apply
    #data_route_act = data_route_apply
    iter_kmeans = 20
    zones_cant = 21
    tolerance = 1
    mean_score = 0
    for iter in 1:iter_kmeans
        scores_arr = apply_markov(data_route, data_route_act, data_travel_times, zones_cant, tolerance)
        #println("SCORES: "*string(scores_arr))
        mean_score = mean_score + mean(scores_arr)
        println("MEAN SCORE: "*string(mean(scores_arr)))
    end
    println("MEAN SCORE TOTAL: "*string(mean_score/iter_kmeans))
end

function solve(alg::MarkovApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3, 4, 5}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #procedure_no = dParams[5] # {1, 2} = {Procedural, Library}
    #solver_id = dParams[6] # {1, 2, 3} = {GLPK, CPLEX, Gurobi}
    #expansion_constant = dParams[7] # {1, 2, 5}
    #inflation_constant = dParams[8] # {1, 2, 5}
    #prune_tolerance = dParams[9] # {1, 3}
    dParams = [2, 21, 20, 1, 2, 2, 1]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::MarkovApproach, routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    
    for R in routes
        #stops_cost = get_stops_coords(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)

        if(dParams[5] == 1)
            ass_stops, cl_sizes = markov_procedural(stops_cost, params)
        else
            ass_stops, cl_sizes = markov_library(stops_cost, params)
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