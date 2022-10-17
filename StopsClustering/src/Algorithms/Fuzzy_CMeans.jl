# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    FuzzyApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of K-Medoids algorithms through explicit procedural and 
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

struct FuzzyApproach <: AbstractSegmentationAlgorithm end

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

function get_fuzzy_clusters(weights)
    num_stops = size(weights,2)
    clusters = Array{Int64}(undef, num_stops)
    for i in 1:num_stops
        _, clusters[i] = findmax(weights[i,:])
    end
    return clusters
end

function fuzzy_procedural(stops_cost, dparams)

    nrows, ncols = size(stops_cost)
    weights = [[mean(stops_cost[:,i]) for i in 1:ncols] for j in 1:dparams[2]] #rand(Float64, ncols, dparams[2])
    weights ./= sum(weights, dims=2)
    fuzziness = dparams[7]
    centers = zeros(nrows, (nrows, dparams[2]))
    prev_centers = identity.(centers)

    delta = Inf
    iter = 0
    while iter < dparams[3] && delta > dparams[4]
        
        wnrows, wncols = size(weights)
        T = eltype(centers)
        for j in 1:wncols
            num = zeros(T, size(stops_cost,1))
            den = zero(T)
            for i in 1:wnrows
                deltam = weights[i,j]^fuzziness
                num += deltam * stops_cost[:,i]
                den += deltam
            end
            centers[:,j] = num/den
        end

        pow = 2.0/(fuzziness-1)
        dists = pairwise(Euclidean(), stops_cost, centers, dims=2)
        for i in 1:wnrows
            for j in 1:wncols
                den = 0.0
                for k in 1:wncols
                    den += (dists[i,j]/dists[i,k])^pow
                end
                weights[i,j] = 1.0/den
            end
        end

        delta = maximum(colwise(Euclidean(), prev_centers, centers))
        copyto!(prev_centers, centers)
        iter += 1
    end

    clusters = get_fuzzy_clusters(weights)
    cluster_sizes = [length(findall(x->x==i,clusters)) for i in 1:length(clusters)]
    return clusters, cluster_sizes, centers 
end

function fuzzy_library(stops_cost, dparams)
    clusters = fuzzy_cmeans(stops_cost, dparams[2], dparams[7]; maxiter=dparams[3], tol=dparams[4])
    #@assert nclusters(R1) == dparams[2] # verify the number of clusters
    assignment_stops = assignments(clusters) # get the assignments of points to clusters
    clust_sizes = counts(clusters) # get the cluster sizes
    clust_centers = clusters.centers
    return assignment_stops, clust_sizes, clust_centers
end

function apply_fuzzy(data_route, data_route_act, tt_data, zones_cant, tolerance)
    scores_arr = Vector{Float64}()
    fuzziness = 2
    for route_num in 1:length(data_route)
    
        TT = tt_data[route_num].travel_times
        TT = TT[Not(name(station(data_route_act[route_num]))),Not(name(station(data_route_act[route_num])))]
        
        R = fuzzy_cmeans(TT, zones_cant, fuzziness; maxiter=200, tol=tolerance)
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
    
function preparation_fuzzy(data_route, data_route_act)
    #data_route = data_mod_apply #data_mod_apply
    #data_route_act = data_route_apply
    iter_kmeans = 20
    zones_cant = 21
    tolerance = 1
    mean_score = 0
    for iter in 1:iter_kmeans
        scores_arr = apply_fuzzy(data_route, data_route_act, data_travel_times, zones_cant, tolerance)
        #println("SCORES: "*string(scores_arr))
        mean_score = mean_score + mean(scores_arr)
        println("MEAN SCORE: "*string(mean(scores_arr)))
    end
    println("MEAN SCORE TOTAL: "*string(mean_score/iter_kmeans))
end

function fuzzy_jump(stops_cost, dparams, opt_solver=Gurobi)
    # Not needed
    return [], [], []
end

function selectSolver(solver_id)
    solver = JuMP.Gurobi
    if(solver_id == 1)
        solver = JuMP.GLPK
    end
    if(solver_id == 2)
        solver = JuMP.CPLEX
    end
    return solver
end

function solve(alg::FuzzyApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3, 4, 5}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #procedure_no = dParams[5] # {1, 2, 3} = {Procedural, Library, JuMP}
    #solver_id = dParams[6] # {1, 2, 3} = {GLPK, CPLEX, Gurobi}
    #fuzziness = dParams[7] # {2, 4, 6}
    dParams = [2, 21, 20, 1, 2, 3, 2]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::FuzzyApproach, routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    
    for R in routes
        #stops_cost = get_stops_coords(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)

        if(dParams[5] == 1)
            ass_stops, cl_sizes, cl_centers = fuzzy_procedural(stops_cost, params)
        end
        if(dParams[5] = 2)
            ass_stops, cl_sizes, cl_centers = fuzzy_library(stops_cost, params)
        end
        if(dParams[5] = 3)
            solver = selectSolver(dParams[6])
            ass_stops, cl_sizes, cl_centers = fuzzy_jump(stops_cost, params, opt_solver=solver)
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