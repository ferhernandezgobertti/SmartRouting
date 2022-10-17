# 1. Procedural - OK
# 2. Library - OK
# 3. JuMP - OK
# 4. Declarative

"""
    HierarchicalApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of Hierarchical Clustering algorithms with Aggregation 
or Divisive procedures through explicit procedural and declarative approaches, 
as well as with optimized libraries and JuMP solvers.

Inputs: 
* routes: List of Transports Routes, i.e. Vector of Vector of Stops (Lat, Lon).
* traveltimes: Matrix associated with times for travelling from one stop to
               another (all combinations), such as expected tt or static tt.

[Optional]
* params: Relevant parameters to configure Identification process.

Outputs:
* routes: List of updated routes (i.e. stops with zones generic identifiers).
* dzones: List of zones (i.e. Segmented Stops)
* zones_quant_arr: List of Zones Ammount per route used on stops.

Note: In this step, all zone generic identifiers are direct (e.g. 1, 2, 3, ...)
"""

using CSV, DataFrames, Dates, Stop
using Clustering, StatsBase, NamedArrays, LinearAlgebra
using JuMP, GLPK, CPLEX, Gurobi

struct HierarchicalApproach <: AbstractSegmentationAlgorithm end

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

function process_hierarchical_model(stops_cost, dparams)
    clusters = Vector{Int}[]
    nnodes = length(hclust.labels)
    nodes = [[i::Int] for i=1:nnodes]
    N = nnodes - dparams[2]
    i = 1
    while i<=N && hclust.height[i] <= dparams[9]
        both = vec(hclust.merge[i,:])
        new_cluster = Int[]
            for x in both
                if x<0
                    push!(new_cluster, -x)
                    nodes[-x] = []
                else
                    append!(new_cluster, clusters[x])
                    clusters[x] = []
                end
            end
        push!(clusters, new)
        i += 1
    end
    all_clusters_ass = vcat(clusters, nodes)
    all_clusters_ass = all_clusters_ass[map(length, all_clusters_ass) .> 0]
    ass_clusters = Array(Int, nnodes)
    for (i,cl) in enumerate(all_clusters_ass)
        ass_clusters[cl] = i
    end
    cluster_sizes = [length(findall(x->x==i,ass_clusters)) for i in 1:length(ass_clusters)]

    return ass_clusters, cluster_sizes 
end

function op_mindist(mindist, cost, i, j)
    return mindist>cost[i,j]
end

function op_avedist(avedist, cost, i, j)
    return avedist>(sum(cost[i,:])+sum(cost[j,:]))/2
end

function op_maxdist(maxdist, cost, i, j)
    return maxdist>cost[i,j]
end

function div_hierarchical(stops_cost, dparams, func_dist)
    clusters = ones(size(stops_cost,1))
    iter_count = 0
    refdist = Inf
    bigger_cluster = 0
    tol = dparams[4]
    while iter_count<dparams[3] && length(unique(clusters))<dparams[2]
        num_clust = length(unique(clusters))
        for i,j in 1:size(stops_cost,1)
            if(i!=j && func_dist(refdist, stops_cost, i, j) && clusters[i]==clusters[j])
                refdist = stops_cost[i,j]
                bigger_cluster = clusters[i]
            end
        end
        for i,j in 1:size(stops_cost,1)
            if(i!=j && func_dist(refdist+tol, stops_cost, i, j) && clusters[i] == bigger_cluster)
                clusters[i] = num_clust + 1
            end
        end
        iter_count = iter_count + 1
    end
    cluster_sizes = [length(findall(x->x==i,clusters)) for i in 1:length(clusters)]
    return clusters, cluster_sizes
end

function mindist(stops_cost, cluster1, cluster2)
    mindist_val = Inf
    for i in cluster1 
        for j in cluster2
            if stops_cost[i,j] < mindist_val
                mindist_val = stops_cost[i,j]
            end
        end 
    end
    return mindist_val
end

function avedist(stops_cost, cluster1, cluster2)
    avedist_val = zero(size(stops_cost,1))
    for i in cluster1 
        for j in cluster2
            s += d[i,j]
        end 
    end
    return avedist_val / (length(cluster1)*length(cluster2)) 
end

function maxdist(stops_cost, cluster1, cluster2)
    maxdist_val = -Inf
    for i in cluster1 
        for j in cluster2
            if stops_cost[i,j] > maxdist_val
                maxdist_val = stops_cost[i,j]
            end
        end 
    end
    return maxdist_val
end

function merge_hierarchical(stops_cost, dparams, func_dist)
    nc = size(stops_cost,1)             # cluster number
    mr = Array(Int, nc-1)               # min row
    mc = Array(Int, nc-1)               # min col
    h = Array(T, nc-1)                  # height
    clusters = [[x] for x in 1:nc]      # clusters
    merges = -[1:nc]
    next = 1
    i = 1
    N = Array(Int, nc+1)
    N[1] = 1                            # arbitrary choice
    while nc > 1
        found = false
        mindist = Inf
        while !found
            i += 1
            mi = 0
            mindist = Inf
            Nim1 = N[i-1]
            ## c[i] = nearest neigbour c[i-1]
            for j = 1:nc 
                if Nim1 != j
                    distance = func_dist(stops_cost, clusters[Nim1], clusters[j])
                    if distance < mindist
                        mindist = distance
                        mi = j
                    end
                end 
            end
            N[i] = mi           # N[i+1] is nearest neigbor to N[i]
            found = i > 2 && N[i] == N[i-2]
        end
        ## merge c[i] and nearest neigbor c[i], i.e., c[i-1]
        if N[i-1] < N[i]
            lo, high = N[i-1], N[i]
        else
            lo, high = N[i], N[i-1]
        end
        ## Store of result
        mr[next] = merges[lo]
        mc[next] = merges[high]
        h[next] = mindist
        merges[lo] = next
        merges[high] = merges[nc]
        next += 1
        ## Cluster Merge
        clusters[lo] = vcat(clusters[lo], clusters[high])
        clusters[high] = clusters[nc]
        if i>3
            i -= 3
        else
            i = 1
        end
        ## replace any nearest neighbor referring to cluster number
        for k=1:i
            if N[k] == nc
                N[k] = high
            end
        end
        nc -= 1
    end
    cluster_sizes = [length(findall(x->x==i,clusters)) for i in 1:length(clusters)]
    return clusters, cluster_sizes
end

function hierarchical_procedural(stops_cost, dparams)
    if(dparams[5] == 1)
        if dparams[7] == 1
            clusters, cl_sizes = merge_hierarchical(stops_cost, dparams, mindist)
        elseif dparams[7] == 2
            clusters, cl_sizes = merge_hierarchical(stops_cost, dparams, avedist)
        else
            clusters, cl_sizes = merge_hierarchical(stops_cost, dparams, maxdist)
        end
    else
        if dparams[7] == 1
            clusters, cl_sizes = div_hierarchical(stops_cost, dparams, op_mindist)
        elseif dparams[7] == 2
            clusters, cl_sizes = div_hierarchical(stops_cost, dparams, op_avedist)
        else
            clusters, cl_sizes = div_hierarchical(stops_cost, dparams, op_maxdist)
        end
    end
    return clusters, cl_sizes
end

function select_linkage(link_type)
    link = :complete # MaxDist
    if(link_type == 1)
        link = :single  # MinDist
    end
    if(link_type == 2)
        link = :average # AveDist
    end
    return link
end

function hierarchical_library(stops_cost, dparams)
    cluster_linkage = select_linkage(dparams[7])
    if(dparams[8] = 1)
        tree_model = hclust(stops_cost, linkage=cluster_linkage, uplo=:U)
    end
    if(dparams[8] = 2)
        tree_model = hclust(stops_cost, linkage=cluster_linkage, uplo=:L)
    end
    if(dparams[8] = 3)
        tree_model = hclust(stops_cost, linkage=cluster_linkage)
    end
    if(dparams[9] = 0)
        clusters = cutree(R, k=dparams[2]; maxiter=dparams[3], tol=dparams[4])
    else
        clusters = cutree(R, k=dparams[2], h=dparams[9]; maxiter=dparams[3], tol=dparams[4])
    end
    
    #@assert nclusters(R1) == dparams[2] # verify the number of clusters
    assignment_stops = assignments(clusters) # get the assignments of points to clusters
    clust_sizes = counts(clusters) # get the cluster sizes
    return assignment_stops, clust_sizes
end

function apply_hierarchical_vmeasure(data_route, data_route_act, tt_data, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
        
        zones_cant = length(zones(data_route[route_num]))
    
        TT = tt_data[route_num].travel_times
        #TT = TT[Not(name(station(data_route_act[route_num]))),Not(name(station(data_route_act[route_num])))]
        
        #DEL stops_0 = data_route_act[route_num].stops[findall(x->zone(x)=="0", data_route_act[route_num].stops)]
        #DEL for S in stops_0
        #DEL    TT = TT[Not(name(S)),Not(name(S))]
        #DEL end
        
        R = hclust(TT, linkage=:single, uplo=:U)
        a = cutree(R, k=zones_cant)

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
        #println("SCORE of ROUTE "*string(route_num)*": "*string(score_reached))
    end
    return scores_arr
end

function apply_hierarchical(data_route, data_route_act, tt_data, zones_cant, tolerance)
    scores_arr = Vector{Float64}()

    for route_num in 1:length(data_route)
    
        TT = tt_data[route_num].travel_times
        #TT = TT[Not(name(station(data_route_act[route_num]))),Not(name(station(data_route_act[route_num])))]
        stops_0 = data_route_act[route_num].stops[findall(x->zone(x)=="0", data_route_act[route_num].stops)]
        for S in stops_0
            TT = TT[Not(name(S)),Not(name(S))]
        end
        
        R = hclust(TT, linkage=:single, uplo=:L)
        a = cutree(R, k=zones_cant)

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
    
function preparation_hierarchical(data_route, data_route_act)
    #data_route = data_mod_apply #data_mod_apply
    #data_route_act = data_route_apply
    iter_kmeans = 20
    zones_cant = 21
    tolerance = 1
    mean_score = 0
    for iter in 1:iter_kmeans
        scores_arr = apply_hierarchical(data_route, data_route_act, data_travel_times, zones_cant, tolerance)
        #println("SCORES: "*string(scores_arr))
        mean_score = mean_score + mean(scores_arr)
        println("MEAN SCORE: "*string(mean(scores_arr)))
    end
    println("MEAN SCORE TOTAL: "*string(mean_score/iter_kmeans))
end

function haversine(lat1, long1, lat2, long2, r = 6372.8)
    lat1, long1 = deg2rad(lat1), deg2rad(long1)
    lat2, long2 = deg2rad(lat2), deg2rad(long2)
    hav(a, b) = sin((b - a) / 2)^2
    inner_term = hav(lat1, lat2) + cos(lat1) * cos(lat2) * hav(long1, long2)
    d = 2 * r * asin(sqrt(inner_term))
    # Round distance to nearest kilometer.
    return round(Int, d)
end

function hierarchical_jump(stops_cost, dparams, opt_solver=Gurobi)
    
    # JuMP PreProcessing
    if(dparams[8] == 1)
        dm = LinearAlgebra.UpperTriangular(stops_cost)
            #[haversine(cities.lat[i], cities.lon[i], cities.lat[j], cities.lon[j]) for i in 1:n, j in 1:n])
    else
        dm = LinearAlgebra.LowerTriangular(stops_cost)
            #[haversine(cities.lat[i], cities.lon[i], cities.lat[j], cities.lon[j]) for i in 1:n, j in 1:n])
    end

    # JuMP Modelling
    model = Model(opt_solver.Optimizer)
    @variable(model, cluster_ass[1:size(stops_cost,1), 1:dparams[2]], Bin)
    @constraint(model, [i = 1:n], sum(cluster_ass[i, :]) == 1)
    fix(cluster_ass[1, 1], 1; force = true) # reducing symetry
    @variable(model, -3 <= cost_diff[1:dparams[2]] <= 3)
    @constraint(model, cost_diff .== cluster_ass' * stops_cost .- dparams[4])
    @variable(model, total_distance[i = 1:n, j = 1:i], Bin)
    for k in 1:dparams[2], i in 1:size(stops_cost,1), j in 1:i
        @constraint(model, total_distance[i, j] >= cluster_ass[i, k] + cluster_ass[j, k] - 1)
    end
    @objective(model, Min, sum(dm[i, j] * total_distance[i, j] for i in 1:n, j in 1:i))
    optimize!(model)

    # JuMP PostProcessing
    cluster_sizes = [length(findall(x->x==i,cluster_ass)) for i in 1:length(cluster_ass)]
    return cluster_ass, cluster_sizes
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

function solve(alg::HierarchicalApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3, 4, 5}
    #expected_cluster_no = dParams[2] # {12, 21, 48}
    #iterations_quantity = dParams[3] # {10, 15, 20}
    #tolerance = dParams[4] # {0, 1, 5}
    #procedure_no = dParams[5] # {1, 2, 3, 4} = {Procedural Agglomerative, Procedural Divisive, Library, JuMP}
    #solver_id = dParams[6] # {1, 2, 3} = {GLPK, CPLEX, Gurobi}
    #cluster_linkage = dParams[7] # {1, 2, 3} = {MinDist, AveDist, MaxDist}
    #symetry_matrix = dParams[8] # {1, 2, 3} = {Upper Triangular, Lower Triangular, Symetrical}
    #max_tree_height = dParams[9] # {0, N} = {None, [175, 221, ...]} Dendogram
    dParams = [2, 21, 20, 1, 2, 3, 2, 1, 175]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::HierarchicalApproach, routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    
    for R in routes
        #stops_cost = get_stops_coords(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)

        if(dParams[5] == 4)
            ass_stops, cl_sizes = hierarchical_jump(stops_cost, params)
        elseif(dParams[5] == 3)
            ass_stops, cl_sizes = hierarchical_library(stops_cost, params)
        else
            ass_stops, cl_sizes = hierarchical_procedural(stops_cost, params)
        end
        push!(zones_quant_arr, length(unique(ass_stops)))
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