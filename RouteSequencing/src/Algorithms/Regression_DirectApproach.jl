

"""
MinTimeDirect <: AbstractSequencingAlgorithm

Algorithm that decides the stops sequence based on auxiliary minimization optimization
algorithm guided by complete zones sequencing modelling through regression using specified 
training route data with filtering capabilities.

Inputs: 
* trainroutes: List of Training Transports Routes (unsequenced).
# trainsequence: List of Training Transport Sequences.
* valroutes: List of Validation/Testing Transports Routes (unsequenced).

or 
* ord_trainroutes: List of Ordered Training Transports Routes (sequenced).
* valroutes: List of Validation/Testing Transports Routes (unsequenced).

[Optional]
* params: Relevant parameters to configure Sequencing process.

Outputs: routes, scores
* routes: List of updated routes (i.e. ordered stops by applying this method).
* scores: List of sequences scores (i.e. evaluation of indicates routes through 
competition scoring algorithm)

"""

using CSV, DataFrames, Dates, RouteSequencing
using LazySets, IOCapture, Combinatronics
using StatsBase, StaticArrays, NamedArrays, LinearAlgebra

struct MinTimeDirect{T} <: AbstractSequencingAlgorithm
    algz::T
end

MinTimeDirect() = MinTimeDirect(MinTime())

struct ClusterOrdering
    alpha_ord::Vector{String}
    beta_ord::Vector{String}
    gamma_ord::Vector{String}
    delta_ord::Vector{String}
    city_id::Int64
    repetitions::Int64
end

ClusterOrdering() = ClusterOrdering(Vector{String}(),Vector{String}(),Vector{String}(),Vector{String}(),0,0)

function value_in_between(value, min_val, max_val)
    return value >= min_val && value <= max_val
end

function stops_in_region(stops, lat_range, lng_range)
    return !isnothing(findfirst(x-> value_in_between(x.lat, lat_range[1], lat_range[2]) && value_in_between(x.lng, lng_range[1], lng_range[2]), stops))
end

function identify_city(route)
    route_city_id = 0
    if(stops_in_region(x.stops, [30.134, 30.517], [-97.951, -97.531]))          # Austin
        route_city_id = 1
    elseif(stops_in_region(x.stops, [42.235, 42.397], [-71.199, -70.967]))      # Boston
        route_city_id = 2
    elseif(stops_in_region(x.stops, [41.653, 42.015], [-87.866, -87.487]))      # Chicago
        route_city_id = 3
    elseif(stops_in_region(x.stops, [33.748, 34.312], [-118.639, -118.152]))    # Los Angeles
        route_city_id = 4
    elseif(stops_in_region(x.stops, [47.502, 47.732], [-122.435, -122.244]))    # Seattle
        route_city_id = 5
    elseif(stops_in_region(x.stops, [-34.922, -34.705], [-56.433, -56.014]))    # Montevideo
        route_city_id = 6
    end
    return route_city_id
end

function is_segment_repeated(ord_zones_segment, clust_ord, repetition_tolerance)
    clust_ord_segment = [clust_ord.alpha_ord[i]*"-"*clust_ord.beta_ord[i]*"."*clust_ord.gamma_ord[i]*clust_ord.delta_ord[i] for i in 1:length(clust_ord.delta_ord)]
    return length(findall(x->x in clust_ord_segment, ord_zones_segment)) >= repetition_tolerance
end

function get_repeated_zones_segment(ord_zones_segment, clusters_ordering, repetition_tolerance)
    return findfirst(x->is_segment_repeated(ord_zones_segment, x, repetition_tolerance), clusters_ordering)
end

function add_repetition_value(clusters_ordering, pos_repeated, route_score, repetition_calculation)
    if(repetition_calculation == 1)
        clusters_ordering[pos_repeated].repetitions = clusters_ordering[pos_repeated].repetitions * route_score + route_score
    elseif(repetition_calculation == 2)
        clusters_ordering[pos_repeated].repetitions = clusters_ordering[pos_repeated].repetitions * route_score
    else
        clusters_ordering[pos_repeated].repetitions = clusters_ordering[pos_repeated].repetitions + route_score
    end
    return clusters_ordering
end

function are_zoneids_differences_appropiate(ord_zones_segment, difference_calculation, diff_macro_ref, diff_micro_ref)
    if(difference_calculation == 1)
        difference_macro = prod([cmp(split(ord_zones_segment[i+1], ".")[1], split(ord_zones_segment[i], ".")[1]) for i in 1:(length(ord_zones_segment)-1)])
        difference_micro = prod([cmp(split(ord_zones_segment[i+1], ".")[2], split(ord_zones_segment[i], ".")[2]) for i in 1:(length(ord_zones_segment)-1)])
    elseif(difference_calculation == 2)
        difference_macro = mean([cmp(split(ord_zones_segment[i+1], ".")[1], split(ord_zones_segment[i], ".")[1]) for i in 1:(length(ord_zones_segment)-1)])
        difference_micro = mean([cmp(split(ord_zones_segment[i+1], ".")[2], split(ord_zones_segment[i], ".")[2]) for i in 1:(length(ord_zones_segment)-1)])
    else
        difference_macro = sum([cmp(split(ord_zones_segment[i+1], ".")[1], split(ord_zones_segment[i], ".")[1]) for i in 1:(length(ord_zones_segment)-1)])
        difference_micro = sum([cmp(split(ord_zones_segment[i+1], ".")[2], split(ord_zones_segment[i], ".")[2]) for i in 1:(length(ord_zones_segment)-1)])
    end
    return difference_macro < diff_macro_ref & difference_micro < diff_micro_ref
end

function process_training_routes(ord_routes, params)
    clusters_ordering = Vector{ClusterOrdering}()
    for R in ord_routes
        route_city = identify_city(R)
        ord_zones = zones(R)
        ord_zones_segment = ord_zones[params[3]:params[4]]
        for i in 1:fld(length(ord_zones), params[4])
            pos_repeated = get_repeated_zones_segment(ord_zones_segment, clusters_ordering, params[5])
            if(!isnothing(pos_repeated))
                clusters_ordering = add_repetition_value(clusters_ordering, pos_repeated, score(R), params[6])
            elseif(are_zoneids_differences_appropiate(ord_zones_segment, params[7], params[8], params[9]))
                alpha_ord = [ord_zones_segment[i][1:2] for i in 1:params[4]]
                beta_ord = [split(ord_zones_segment[i][3:end],".")[1] for i in 1:params[4]]
                gamma_ord = [split(ord_zones_segment[i],".")[2][1:2] for i in 1:params[4]]
                delta_ord = [split(ord_zones_segment[i],".")[2][3:4] for i in 1:params[4]]
                push!(clusters_ordering, ClusterOrdering(alpha_ord, beta_ord, gamma_ord, delta_ord, route_city, 0))
            end
            ord_zones_segment = ord_zones[params[4]*i:params[4]*(i+1)]
        end
    end
    return clusters_ordering
end

function get_zone_value(zone_id, params)
    return params[15]*zone_id[1:2] + params[16]*split(zone_id[3:end],".")[1] + params[17]*split(zone_id,".")[2][1:2] + params[18]*split(zone_id,".")[2][3:4]
end

function get_initial_sequencing(route_zones, params)
    ord_route_zones = Vector{String}()
    zones_values = [get_zone_value(route_zones[i], params) for i in 1:length(route_zones)]
    for i in 1:length(route_zones)
        _, min_pos = findmin(zones_values)
        push!(ord_route_zones, route_zones[min_pos])
        deleteat!(route_zones, min_pos)
    end
    return ord_route_zones
end

function get_all_repeated_zones_segment(ord_zones_segment, clusters_ordering, repetition_tolerance)
    return findall(x->is_segment_repeated(ord_zones_segment, x, repetition_tolerance), clusters_ordering)
end

function get_zones_order(route_zones, clusters_ordering, params)
    zones_sequence = Vector{String}()
    route_zones = get_initial_sequencing(route_zones, params)
    route_seg_zones = route_zones[params[10]:params[11]]
    for i in 1:fld(length(route_zones),params[11])
        repeated_zones = get_all_repeated_zones_segment(route_seg_zones, clusters_ordering, params[14])
        quant_coincidents, _ = findmax([clusters_ordering[i].repetitions for i in repeated_zones])
        if(quant_coincidents >= params[12])
            rem_zones = get_remaining_zones(repeated_zones, route_seg_zones)
            repeated_zones2 = get_all_repeated_zones_segment(rem_zones, clusters_ordering, params[14])
            quant_coincidents2, _ = findmax([clusters_ordering[i].repetitions for i in repeated_zones2])
            if(quant_coincidents2 >= params[13]) 
                [push!(repeated_zones, repeated_zones2[i]) for i in 1:length(repeated_zones2)]
            end
            [push!(zones_sequence, repeated_zones[i]) for i in 1:length(repeated_zones)]
        end
        route_seg_zones = route_zones[params[11]*i:params[11]*(i+1)]
    end
    return zones_sequence
end

function get_remaining_zones(zones_route, zones_sequence)
    return zones_route[findall(x->!(x in zones_sequence), zones_route)]
end

function regression_cluster_ordering(training_routes::Vector{Route}, training_sequences::Vector{Sequence}, validation_routes::Vector{Route})
    # training_routes_quantity = dParams[1] {0, 1000, 3000, 5000} = {All, N Routes}
    # validation_routes_quantity = dParams[2] {0, 20, 35, 50} = {All, N Routes}
    # train_zones_segment_offset = dParams[3] {1, 2, 3}
    # train_zones_segment_limit = dParams[4] {7, 9, 10}
    # train_repetition_tolerance = dParams[5] {2, 4, 5}
    # repetition_calculation = dParams[6] {1, 2, 3} = {Multiplication + Sum, Multiplication, Sum}
    # difference_calculation = dParams[7] {1, 2, 3} = {Prod Diff, Mean Diff, Sum Diff}
    # diff_macro_ref = dParams[8] {1, 2, 3}
    # diff_micro_ref = dParams[9] {3, 5, 7}
    # val_zones_segment_offset = dParams[10] {1, 2, 3}
    # val_zones_segment_limit = dParams[11] {7, 9, 10}
    # first_round_compare_val = dParams[12] {12, 14, 18}
    # second_round_compare_val = dParams[13] {7, 9, 10}
    # val_repetition_tolerance = dParams[14] {2, 4, 5}
    # alpha_param = dParams[15] {8, 10, 12}
    # beta_param = dParams[16] {5, 6, 8}
    # gamma_param = dParams[17] {2, 4, 5}
    # delta_param = dParams[18] {0.8, 1, 2}
    dParams = [0, 0, 1, 7, 4, 3, 3, 1, 3, 1, 7, 14, 9, 2, 10, 8, 4, 2]
    [apply_sequence!(training_routes[i], training_sequences[i]) for i in 1:length(training_routes)]
    regression_cluster_ordering(training_routes, validation_routes; params=dParams)
end

function regression_cluster_ordering(ord_training_routes::Vector{Route}, validation_routes::Vector{Route}; params)
    routes_zones_sequence = Vector{Vector{String}}()
    routes_rem_zones_sequence = Vector{Vector{String}}()
    params[1] != 0 && ord_training_routes = ord_training_routes[1:params[1]]
    params[2] != 0 && validation_routes = validation_routes[1:params[2]]
    clusters_ord = process_training_routes(ord_training_routes, params)
    for R in validation_routes
        route_city = identify_city(R)
        zones_sequence = get_zones_order(zones(R), clusters_ord[findall(x->x.city_id == route_city, clusters_ord)], params)
        push!(routes_zones_sequence, zones_sequence)
        rem_zones_sequence = get_remaining_zones(zones(R), zones_sequence)
        push!(routes_rem_zones_sequence, rem_zones_sequence)
    end
    return routes_zones_sequence, routes_rem_zones_sequence
end

function solve(alg::MinTimeDirect, valroute::AbstractRoute, times::TravelTimes, clustord::Vector{ClusterOrdering}; initial_stop=get_station(route), params)
    T = times.travel_times
    return solve(alg, valroute, T, zones_list=get_zones_order(zones(R), clustord[findall(x->x.city_id == route_city, clustord)], params), initial_stop=initial_stop)
end

function solve(alg::MinTimeDirect, route::AbstractRoute, T::AbstractMatrix; zones_list=get_zones_order(route), initial_stop=get_station(route))
    route_id = name(route)
    #T = times.travel_times
    algz = alg.algz
    # zonas que quedan por recorrer
    waiting_list = zones_list
    idx = findfirst(x -> x == "0", waiting_list)
    sequence = Sequence(route_id, Vector{BareStop}())
    while length(waiting_list) > 0
        # filtrar solo las paradas de la zona actual
        actual_zone = waiting_list[idx]
        idx_stops = findall(s -> split(zone(s), ".")[1] == actual_zone, route.stops)
        # obtener las paradas y la matriz de tiempos asociada a esa zona
        subroute = BareRoute(route_id, route.stops[idx_stops])
        # indexamos la matriz de tiempos por nombres de parada
        name_stops = [name(s) for s in route.stops[idx_stops]]
        subtimes = view(T, name_stops, name_stops)
        # calcular la secuencia de acuerdo a los tiempos minimos
        initial_stop = _next_stop(alg, subroute)
        subsequence = solve(algz, subroute, subtimes, initial_stop=initial_stop)
        # append stops to the result
        for s in subsequence.stops
            push!(sequence, s)
        end
        deleteat!(waiting_list, idx)
    end
    return sequence
end

# simply return the first element of the list
function _next_stop(::MinTimeDirect, subroute)
    return subroute.stops[1]
end

function write_json_sequence(data; path=joinpath("..", "data", "model_apply_outputs", "proposed_sequences.json"))
    # dictionary to write
    dict = sequence_to_dict(data)
    # pass data as a json string (how it shall be displayed in a file)
    stringdata = JSON.json(dict)
    # write the file with the stringdata variable information
    open(path, "w") do f
        write(f, stringdata)
    end;
    return nothing
end

function sequence_to_dict(S::Sequence)
    dict = Dict{String, Any}()
    route_id = "RouteID_" * name(S)
    seq = Dict{String, Int}()
    for (i, stop) in enumerate(stops(S)) # loop over stops, assumed that the stops are sorted
        stop_id = name(stop)
        seq[stop_id] = i-1
    end
    dict[route_id] = Dict("proposed"=>seq)
    return dict
end

function sequence_to_dict(data::Vector{<:Sequence})
    dict = Dict{String, Any}()
    for S in data # loop over sequences
        route_id = "RouteID_" * name(S)
        seq = Dict{String, Int}()
        for (i, stop) in enumerate(stops(S)) # loop over stops, assumed that the stops are sorted
            stop_id = name(stop)
            seq[stop_id] = i-1
        end
        dict[route_id] = Dict("proposed"=>seq)
    end
    return dict
end

function print_score(S::Sequence; src="evaluar.py")
    # convert the sequence to JSON format
    write_json_sequence(S)
    route_id = "RouteID_" * name(S)
    score = run(`python3 $src $route_id`)
    return score
end

function get_score2(S::Sequence)
    out = IOCapture.capture() do
                 print_score(S)
          end
    out_str = split(out.output[1:end-1], "\n")
    result = [parse(Float64, x) for x in out_str]
    return result[1]
end

get_score2(S::Vector{<:Sequence}) = get_score2.(S)

function make_sequencing_process(alg, routes, data_times)
    scores_arr = Vector{Vector{Float64}}()
    sequences_rand = Vector{Sequence}(undef, length(routes))
    for i in [1:1:length(routes);]
        #alg = MinTimeDirect()
        #alg = MinTime()
        # calcular la secuencia
        #sequences_rand[i] = solve(data_route_apply[i], data_times[i].travel_times)
        sequences_rand[i] = solve(alg, routes[i], data_times[i].travel_times)
        # aplicar el algoritmo
        #sequences_rand[i] = solve_MinTimeDirect(data_route_apply[i], data_times[i].travel_times) #initial_zone = "H24.2C") #get_stop(route,"LZ"))
        # obtener el score
        push!(scores_arr, get_score2(sequences_rand))
        #println(get_score2(sequences_rand[i]))
    end
    return sequences_rand, scores_arr
end

function make_compact_sequencing_process(alg, routes, data_times)
    #alg = MinTimeDirect()
    # calcular la secuencia
    return [solve(alg, routes[i], data_times[i].travel_times) for i in 1:length(routes)];
end