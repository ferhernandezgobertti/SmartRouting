"""

    ConstrainedApproach <: AbstractSegmentationAlgorithm

Algorithm that segments the already filtered stops per route while using 
customized versions of Constrained Clustering algorithms based on Clusters Number
and Optimization Criteria considerations through declarative approaches and diverse
JuMP solvers.

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
using JuMP, GLPK, CPLEX, Gurobi

struct ConstrainedApproach <: AbstractSegmentationAlgorithm end

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

function get_cluster_size(clust_ass, stops_coords, clust_num)
    lat_min = minimum(stops_coords[findall(x->x==clust_num, clust_ass),:])
    lat_max = maximum(stops_coords[findall(x->x==clust_num, clust_ass),:])
    lon_min = minimum(stops_coords[:, findall(x->x==clust_num, clust_ass)])
    lon_max = maximum(stops_coords[:, findall(x->x==clust_num, clust_ass)])
    return ((lat_max-lat_min)+(lon_max-lon_min))/2
end

function min_distance_clusters(clust_ass, stops_cost, clust_num)
    stops_on_clust = findall(x->x==clust_num, clust_ass)
    min_distance = Inf
    for stop in stops_on_clust
        for pos in 1:size(stops_cost,2)
            if(!(pos in stops_on_clust) && min_distance>stops_cost[stop][pos])
                min_distance = stops_cost[stop][pos]
            end
        end
    end
    return min_distance
end

function get_distance(coords_a, coords_b, is_euclidean)
    distance_ab = 0
    if(is_euclidean)
        distance_ab = sqrt((coords_b[1]-coords_a[1])^2+(coords_b[2]-coords_a[2])^2)
    else 
        distance_ab = abs(coords_b[1]-coords_a[1])+abs(coords_b[2]-coords_a[2])
    end
    return distance_ab 
end

function stops_quantity_per_radius(clust_ass, stops_coords, stop_pos, radius)
    stops_on_clust = findall(x->x==clust_ass[stop_pos], clust_ass)
    stops_quantity = 0
    for stop_clust in stops_on_clust
        if(get_distance(stops_coords[stop_pos,:],stops_coords[stop_clust,:], true)<=radius)
            stops_quantity = stops_quantity + 1
        end
    end
    return stops_quantity-1
end

function constrained_jump_size(stops_cost, stops_coords, dparams, opt_solver=Gurobi)
    
    # JuMP Preprocessing
    num_stops = size(stops_cost,1)

    # JuMP Modelling
    model = Model(opt_solver.Optimizer)
    @variable(model,clust_ass[1:num_stops],Int64)
    fix(clust_ass[1, 1], 1; force = true) # reducing symetry
    @objective(model,Min,sum([get_cluster_size(clust_ass, stops_coords, clust_num) for clust_num in 1:length(unique(clust_ass))]))

    # Cluster Number Constraints
    @constraint(model,length(unique(clust_ass))>=dParams[2])
    @constraint(model,length(unique(clust_ass))<=dParams[3])

    # Cluster Stops Constraints
    for clust_id in 1:length(unique(clust_ass))
        @constraint(model,count(x->x=clust_id, clust_ass)>=dParams[6])
        @constraint(model,count(x->x=clust_id, clust_ass)<=dParams[7])
    end

    # Delta, Must/Cannot Link Constraints
    for clust_id in 1:length(unique(clust_ass))
        #stops_on_clust = findall(x->x==clust_id, clust_ass)
        @constraint(model,min_distance_clusters(clust_ass, stops_cost, clust_id)>=dParams[8])
        @constraint(model,min_distance_clusters(clust_ass, stops_cost, clust_id)<=dParams[9])
    end

    # Epsilon Constraints
    for i in 1:size(stops_cost,1)
        @constraint(model,stops_quantity_per_radius(clust_ass, stops_coords, i, dParams[10])>=dParams[11])
    end

    optimize!(model)

    # JuMP PostProcessing
    cluster_sizes = [length(findall(x->x==i,clust_ass)) for i in 1:length(clust_ass)]

    return clust_ass, cluster_sizes
end

function constrained_jump_distance(stops_cost, stops_coords, dparams, opt_solver=Gurobi)
    
    # JuMP Preprocessing
    num_stops = size(stops_cost,1)

    # JuMP Modelling
    model = Model(opt_solver.Optimizer)
    @variable(model,clust_ass[1:num_stops],Int64)
    fix(clust_ass[1, 1], 1; force = true) # reducing symetry
    @objective(model,Min,sum([min_distance_clusters(clust_ass, stops_cost, clust_num) for clust_num in 1:length(unique(clust_ass))]))

    # Cluster Number Constraints
    @constraint(model,length(unique(clust_ass))>=dParams[2])
    @constraint(model,length(unique(clust_ass))<=dParams[3])

    # Cluster Sizes Constraints
    for clust_id in 1:length(unique(clust_ass))
        @constraint(model,get_cluster_size(clust_ass, stops_coords, clust_num)>=dParams[4])
        @constraint(model,get_cluster_size(clust_ass, stops_coords, clust_num)<=dParams[5])
    end

    # Cluster Stops Constraints
    for clust_id in 1:length(unique(clust_ass))
        @constraint(model,count(x->x=clust_id, clust_ass)>=dParams[6])
        @constraint(model,count(x->x=clust_id, clust_ass)<=dParams[7])
    end

    # Epsilon Constraints
    for i in 1:size(stops_cost,1)
        @constraint(model,stops_quantity_per_radius(clust_ass, stops_coords, i, dParams[10])>=dParams[11])
    end

    optimize!(model)

    # JuMP PostProcessing
    cluster_sizes = [length(findall(x->x==i,clust_ass)) for i in 1:length(clust_ass)]

    return clust_ass, cluster_sizes
end

function constrained_jump_inertia(stops_cost, stops_coords, dparams, opt_solver=Gurobi)
    
    # JuMP Preprocessing
    num_stops = size(stops_cost,1)

    # JuMP Modelling
    model = Model(opt_solver.Optimizer)
    @variable(model,clust_ass[1:num_stops],Int64)
    fix(clust_ass[1, 1], 1; force = true) # reducing symetry
    @objective(model,Min,sum([get_cluster_size(clust_ass, stops_coords, clust_num)/min_distance_clusters(clust_ass, stops_cost, clust_num) for clust_num in 1:length(unique(clust_ass))]))

    # Cluster Number Constraints
    @constraint(model,length(unique(clust_ass))>=dParams[2])
    @constraint(model,length(unique(clust_ass))<=dParams[3])

    # Cluster Stops Constraints
    for clust_id in 1:length(unique(clust_ass))
        @constraint(model,count(x->x=clust_id, clust_ass)>=dParams[6])
        @constraint(model,count(x->x=clust_id, clust_ass)<=dParams[7])
    end

    # Epsilon Constraints
    for i in 1:size(stops_cost,1)
        @constraint(model,stops_quantity_per_radius(clust_ass, stops_coords, i, dParams[10])>=dParams[11])
    end

    optimize!(model)

    # JuMP PostProcessing
    cluster_sizes = [length(findall(x->x==i,clust_ass)) for i in 1:length(clust_ass)]

    return clust_ass, cluster_sizes
end

function solve(alg::ConstrainedApproach, routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
    #cost_no= dParams[1] # {1, 2, 3, 4, 5}
    #min_cluster_num = dParams[2] # {12, 21, 48}
    #max_cluster_num = dParams[3] # {12, 21, 48}
    #min_cluster_size = dParams[4] # {0.5, 0.75, 1}
    #max_cluster_size = dParams[5] # {1.25, 1.5, 2}
    #min_cluster_stops = dParams[6] # {1, 2, 4}
    #max_cluster_stops = dParams[7] # {6, 8, 10}
    #min_cluster_split = dParams[8] # {0.03, 0.05, 0.08} aka Delta Constraint
    #max_cluster_split = dParams[9] # {0.17, 0.22, 0.25} aka Must/Cannot Link Constraint
    #radius_near_stops_same_cluster = dParams[10] # {0.1, 0.3, 0.5} aka Epsilon Constraint
    #radius_near_stops_quantity = dParams[11] # {0, 1, 2}
    #tolerance = dParams[12] # {0, 0.02, 0.05}
    #procedure_no = dParams[13] # {1, 2, 3} = {JuMP Cluster Size, JuMP Cluster Distance, JuMP Size/Distance Inertia}
    #solver_id = dParams[14] # {1, 2, 3} = {GLPK, CPLEX, Gurobi}
    #triangular_matrix = dParams[15] # {1, 2} = {Upper Triangular, Lower Triangular}
    dParams = [2, 21, 20, 1, 2, 3, 5, 20, 1, 1, 1]
    return solve(alg, routes, traveltimes; params=dParams)
end

function solve(alg::ConstrainedApproach, routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    
    for R in routes
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)
        solver = selectSolver(dParams[14])

        if(dParams[13] == 3)
            ass_stops, cl_sizes = constrained_jump_inertia(stops_cost, stops_coords, params, opt_solver=solver)
        elseif(dParams[13] == 2)
            ass_stops, cl_sizes = constrained_jump_distance(stops_cost, stops_coords, params, opt_solver=solver)
        else
            stops_coords = get_stops_coords(R)
            ass_stops, cl_sizes = constrained_jump_size(stops_cost, stops_coords, params, opt_solver=solver)
        end
        push!(zones_quant_arr, length(unique(ass_stops)))
        stops = unique(R.stops, 1)
        for S in stops
            curr_ass = string(ass_stops[findfirst(x->x.id==S.id, stops)])
            S.zone = Zone("", curr_ass)
            push!(zones.ids, curr_ass)
            push!(zones.stops, S)
        end
    end
    return routes, dZones, zones_quant_arr
end