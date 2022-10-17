"""

    minizinc_preprocessing

Algorithm that prepares the stops segmentation data model using the cost matrices 
and relevant modelling paramenters in order to be used in MiniZinc application 
through declarative processing.

Inputs: 
* routes: List of Routes (which contain determined lists of stops).
* route_distance: List of Matrices (one per Route) with approximate route distance [km] 
between stops (considering all combinations Stop_A -> Stop_B ).
* travel_times: List of Matrices (one per Route) with static or expected travel times 
between stops (considering all combinations Stop_A -> Stop_B ).

[Optional]
* params: Relevant parameters to configure Segmentation process.

Output:
* File with .dzn format compatible with MiniZinc containing data model specification.

--------------------------------------------------------------------------------------

    minizinc_postprocessing

Algorithm that further processes the stops segmentation output collected from MiniZinc 
optimization solvers, in order to obtain and define the corresponding stops groups (zones).

Inputs: 
* opt_solver_output: Numerical vector that differentiates stops into several groups.
* routes: List of Routes (which contain determined lists of stops).
* route_no: Model number related to the current route group being processed.

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

# Generic Vector Formulation formatted as custom String
function get_vector_string(vector)
    vector_str = "["
    for i in 1:(length(vector)-1)
        vector_str = vector_str * string(vector[i]) * ", "
    end
    return vector_str*string(vector[length(vector)])*"]"
end

function get_vector_for_matrix(vector)
    vector_str = ""
    for i in 1:(length(vector)-1)
        vector_str = vector_str * string(vector[i]) * ", "
    end
    return vector_str*string(vector[length(vector)])*"|"
end

function get_matrix_string(matrix)
    matrix_str = "["
    for i in 1:(size(matrix,1))
        matrix_str = matrix_str * get_vector_for_matrix(matrix[i,:]) * " "
    end
    return matrix_str*"]"
end

function get_params_string(params)
    params_str = "minnum = "*string(params[2])*"\n maxnum = "*string(params[3])*"\n minsize = "*string(params[4])
    params_str = params_str * "\n maxsize = "*string(params[5])*"\n minstops = "*string(params[6])
    params_str = params_str * "\n maxstops = "*string(params[7])*"\n minsplit = "*string(params[8])
    params_str = params_str * "\n maxsplit = "*string(params[9])*"\n radiuseps = "*string(params[10])
    params_str = params_str * "\n radiusquant = "*string(params[11])*"\n tolerance = "*string(params[12])
    return params_str
end

# Model Content obtention and Data recollection
function get_model_content(stops_cost, stops_coords, dParams)
    num_stops = size(stops_cost,1)
    cost_matrix_str = get_matrix_string(stops_cost)
    lat_vector_str = get_vector_string([stops_coords[i,:] for i in 1:num_stops])
    lon_vector_str = get_vector_string([stops_coords[:,j] for j in 1:num_stops])
    dparams_str = get_params_string(dParams)
    model_str = " nstops = "*string(size(stops_cost,1))*"\n nparams = "*string(length(dParams))
    model_str = model_str*"\n costmatrix = "*string(cost_matrix_str)*"\n lat = "*string(lat_vector_str)
    model_str = model_str*"\n lon = "*string(lon_vector_str)*"\n "*string(dparams_str)
    return model_str
end

# Auxiliary function that saves a DZN file
function save_dzn_file(filepath, str_to_save)
    open(filepath, "w") do dzn_file
        write(dzn_file, str_to_save)
    end
end

function minizinc_preprocessing(routes::Vector{Route}, traveltimes::Vector{Matrix{Float64}})
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
    #radius_near_stops_quantity = dParams[11] # {1, 3, 5} aka Epsilon Constraint
    #tolerance = dParams[11] # {0, 0.02, 0.05}
    dParams = [2, 21, 48, 1, 2, 2, 8, 0.05, 0.25, 0.5, 1, 0.02]
    return solve(alg, routes, traveltimes; params=dParams)
end

function minizinc_preprocessing(routes::Vector{Route}, routedistance::Vector{Matrix{Float64}}, traveltimes::Vector{Matrix{Float64}}; params=dParams)
    model_count = 1
    for R in routes
        stops_coords = get_stops_coords(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        stops_cost = get_stops_cost_matrices(R, routedistance[route_num], traveltimes[route_num], velc[route_num], vela[route_num], veld[route_num], params)
        model_str = get_model_content(stops_cost, stops_coords, params)
        save_dzn_file("./StopsClust_MiniZinc/data/stop_clust_"*string(model_count)*".dzn", model_str)
        model_count = model_count + 1
    end
end

function minizinc_postprocessing(opt_solver_output::Vector{Int64}, routes::Vector{Route}, model_no::Int)

    zones = DZones(Vector{String}(),Vector{Stop}())
    zones_quant_arr = Vector{Int64}()
    R = routes[model_no]
    push!(zones_quant_arr, length(unique(opt_solver_output)))
    stops = unique(R.stops, 1)
    for S in stops
        curr_ass = string(opt_solver_output[findfirst(x->x.id==S.id, stops)])
        S.zone = Zone("", curr_ass)
        push!(zones.ids, curr_ass)
        push!(zones.stops, S)
    end
    return routes, dZones, zones_quant_arr
end