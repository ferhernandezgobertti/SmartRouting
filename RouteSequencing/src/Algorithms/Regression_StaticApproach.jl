"""
MinTimeMacro <: AbstractSequencingAlgorithm

Algorithm that decides the stops sequence based on auxiliary minimization optimization
algorithm guided by static macrozone sequencing modelling using specified training route
data with filtering capabilities.

Inputs: 
* trainroutes: List of Training Transports Routes (unsequenced).
* valroutes: List of Validation/Testing Transports Routes (unsequenced).

or 
* valroutes: List of Validation/Testing Transports Routes (unsequenced).
* zones_order: List of MacroZones Sequences from Training Routes Processing.


Outputs: routes, scores
* routes: List of updated routes (i.e. ordered stops by applying this method).
* scores: List of sequences scores (i.e. evaluation of indicates routes through 
competition scoring algorithm)

"""

using CSV, DataFrames, Dates, RouteSequencing
using LazySets, IOCapture, Combinatronics
using StatsBase, StaticArrays, NamedArrays, LinearAlgebra

struct MinTimeMacro{T} <: AbstractSequencingAlgorithm
    algz::T
end

MinTimeMacro() = MinTimeMacro(MinTime())

struct StationsInfo
    station_lng_arr::Vector{Float64}
    station_lat_arr::Vector{Float64}
end

StationsInfo() = StationsInfo(Vector{Float64}(),Vector{Float64}())

struct PackagesInfo
    rejected_arr::Vector{Int64}
    attempted_arr::Vector{Int64}
    delivered_arr::Vector{Int64}
end

PackagesInfo() = PackagesInfo(Vector{Int64}(),Vector{Int64}(),Vector{Int64}())

function get_zones_seq_str_arr(R)
    zones_list = zones(R)
    zones_list = zones_list[2:length(zones_list)]
    zones_str_arr = Vector{String}()
    for i in 1:zones_max
        if(i<=length(zones_list))
            push!(zones_str_arr, zones_list[i])
        else
            push!(zones_str_arr, "N/A")
        end
    end
    return zones_str_arr 
end

function get_macrozones(zones_list)
    macrozones_list = Vector{String}()
    for Z in zones_list
        push!(macrozones_list, split(Z, ".")[1])
    end
    #println(length(unique(macrozones_list)))
    #if(length(unique(macrozones_list))>40)
    #    println("WARNING")
    #end
    return unique(macrozones_list)
end

function get_macrozones_seq_str_arr(R)
    zones_list = zones(R)
    zones_list = zones_list[2:length(zones_list)]
    zones_list = get_macrozones(zones_list)
    zones_str_arr = Vector{String}()
    for i in 1:zones_max
        if(i<=length(zones_list))
            push!(zones_str_arr, zones_list[i])
        else
            push!(zones_str_arr, "N/A")
        end
    end
    return zones_str_arr
end

function get_station_coord(R)
    station = get_station(R)
    return "("*string(station.lng)*","*string(station.lat)*")"
end

function get_zones_list(zones_arr)
    list_str = ""
    for Z in zones_arr
        list_str = list_str*","*Z
    end
    return list_str
end

function get_packages_info(routes)
    packages_info = PackagesInfo() 
    for R in routes
        num_rejected = 0
        num_attempted = 0
        num_delivered = 0
        for S in R.stops
            for P in S.packages
                if(is_rejected(P))
                    num_rejected += 1
                end
                if(is_attempted(P))
                    num_attempted += 1
                end
                if(is_delivered(P))
                    num_delivered += 1
                end
            end
        end
        push!(packages_info.rejected_arr, num_rejected)
        push!(packages_info.attempted_arr, num_attempted)
        push!(packages_info.delivered_arr, num_delivered)
    end
    return packages_info
end

function get_zones_seq(R, zones_arr)
    zones_seq_str = ""
    zones_list = zones(R)
    zones_list = zones_list[2:length(zones_list)]
    zones_seq = Int.(zeros(length(zones_arr)))
    for i in 1:length(zones_list)
        Z = zones_list[i]
        zones_pos = findfirst(x->x==Z, zones_arr)
        zones_seq[zones_pos] = i
    end
    for S in zones_seq
        zones_seq_str = zones_seq_str*","*string(S)
    end
    return zones_seq_str 
end

function get_csv_data(routes, zones_arr, rejected_arr, attempted_arr, delivered_arr)
    csv_data = "Station,RouteID,Score,NumPackRejected,NumPackAttempted,NumPackDelivered"*get_zones_list(zones_arr)
    for i in 1:length(routes)
        #println(i)
        R = routes[i]
        station_coord = get_station_coord(R)
        csv_data = csv_data*"\n"*string(station_coord)*","*name(R)*","*string(score(R))*","*string(rejected_arr[i])
        csv_data = csv_data*","*string(attempted_arr[i])*","*string(delivered_arr[i])*string(get_zones_seq(R,zones_arr))
    end
    return csv_data
end

function get_zones_max(routes)
    zones_max = 0 # = findmax([length(zones(R)) for R in routes])
    for R in routes
        zones_quantity = length(zones(R))
        if(zones_quantity > zones_max)
            zones_max = zones_quantity
        end
    end
    return zones_max
end

function get_stations_info(routes)
    stations_info = StationsInfo()
    for R in routes
        station = get_station(R)
        station_lng = station.lng
        station_lat = station.lat
        if(!(station_lng in stations_info.station_lng_arr) & !(station_lat in stations_info.station_lat_arr))
            push!(stations_info.station_lng_arr, station_lng)
            push!(stations_info.station_lat_arr, station_lat)
        end
    end
    return stations_info
end

function get_col_arr(df_route)
    col_pos = 1
    col_arr = Vector{String}()
    col_name_arr = names(df_route)
    col_data_arr = eachcol(df_route)
    for i in 1:length(col_name_arr)
        col_name = col_name_arr[i]
        col_data = col_data_arr[i]
        if(sum(col_data)!=0)
            push!(col_arr, col_name)
        end
    end
    return col_arr
end

function write_csv_dataframe(routes, df_routes)
    routes_processed = Vector{String}()
    for i in 1:length(routes)
        #println(i)
        station = get_station(routes[i])
        station_lng = station.lng
        station_lat = station.lat
        df_route = df_routes[(df_routes.StationLng .== station_lng) .& (df_routes.StationLat .== station_lat), :]
        col_arr = Vector{String}()
        col_name_arr = names(df_route)
        col_data_arr = eachcol(df_route)
        for i in 1:length(col_name_arr)
            col_name = col_name_arr[i]
            col_data = col_data_arr[i]
            is_col_empty = true
            for C in col_data
                if(string(C)!="N/A")
                    is_col_empty = false
                    break
                end
            end
            if(!is_col_empty)
                push!(col_arr, col_name)
            end
        end
        df_route_red = df_route[:, filter(x->x in col_arr, names(df_route))]
        route_name = name(data_route_apply[i])
        CSV.write("6-Training_Dataframe_"*string(route_name)*".csv", df_route_red) 
        push!(routes_processed, route_name)
    end
    return routes_processed
end

function generate_routes_dataframe(routes)
    df_routes = DataFrame(StationNum=[0], StationLng=[0.0], StationLat=[0.0], Score=[0.0], NumPackRejected=[0.0], NumPackAttempted=[0.0], NumPackDelivered=[0.0])
    seq_arr = Vector{String}()
    for i in 1:get_zones_max(routes)
        push!(seq_arr,"Z"*string(i))
    end
    for S in seq_arr
        df_routes[!,S] = [""]
    end
    for i in 1:length(routes)
        arr_addition = Vector{String}()
        println(i)
        R = routes[i]
        #station_coord = get_station_tuple(R)
        station = get_station(R)
        station_lng = station.lng
        station_lat = station.lat
        route_score = score(R)
        macrozones_seq_str_arr = get_macrozones_seq_str_arr(R)
        station_lng_num = findfirst(x->x==station_lng, station_lng_arr)
        station_lat_num = findfirst(x->x==station_lat, station_lat_arr)
        #arr_to_add = [string(station_lng_num), string(station_lng), string(station_lat)]
        push!(arr_addition, station_lng_num)
        push!(arr_addition, station_lng)
        push!(arr_addition, station_lat)
        push!(arr_addition, route_score)
        push!(arr_addition, rejected_arr[i])
        push!(arr_addition, attempted_arr[i])
        push!(arr_addition, delivered_arr[i])
        for Z in macrozones_seq_str_arr
            push!(arr_addition, Z)
        end
        #println(length(names(df9)))
        #println(length(arr_addition))
        push!(df_route, arr_addition)
    end
    return df_routes
end

function get_reduced_dataframe(route_obj, df_route)
    station = get_station(route_obj)
    station_lng = station.lng
    station_lat = station.lat
    df_route = df_route[(df_route.StationLng .== station_lng) .& (df_route.StationLat .== station_lat), :]
    col_arr = Vector{String}()
    col_name_arr = names(df_route)
    col_data_arr = eachcol(df_route)
    for i in 1:length(col_name_arr)
        col_name = col_name_arr[i]
        col_data = col_data_arr[i]
        is_col_empty = true
        for C in col_data
            if(string(C)!="N/A")
                is_col_empty = false
                break
            end
        end
        if(!is_col_empty)
            push!(col_arr, col_name)
        end
    end
    df_route_red = df_route[:, filter(x->x in col_arr, names(df_route))]
    return df_route_red
end

function filter_dataframe_with_zones(df_route, filter_zones)
    col_arr = Vector{String}()
    col_name_arr = names(df_route)
    col_data_arr = eachcol(df_route)
    for i in 1:length(col_name_arr)
        col_name = col_name_arr[i]
        col_data = col_data_arr[i]
        is_col_zone = false
        if(occursin("Z",col_name))
            for C in col_data
                if(string(C) in filter_zones)
                    is_col_zone = true
                    break
                end
            end
            if(is_col_zone)
                push!(col_arr, col_name)
            end
        else
            push!(col_arr, col_name)
        end
    end
    df_route_red = df_route[:, filter(x->x in col_arr, names(df_route))]
    return df_route_red
end

function get_zones_order(route_obj)
    zones_obj = get_macrozones(zones(route_obj))
    zones_obj_red = zones_obj[1:length(zones_obj)-1]
    df_route_obj = get_reduced_dataframe(route_obj)
    #show(df_route_obj, allcols=true)
    #df_red_route_obj = df_route_obj[(df_route_obj.Z1 .== zones_obj_red[1]) .| (df_route_obj.Z1 .== zones_obj_red[2]), :]
    df_red_route_obj = filter_dataframe_with_zones(df_route_obj, zones_obj_red)
    #show(df_red_route_obj, allcols=true)
    zones_obj_red_perms = permutations(zones_obj_red) |> collect
    zones_obj_red_cants = zeros(length(zones_obj_red_perms))
    row_data_arr = eachrow(df_red_route_obj)
    for i in 1:length(zones_obj_red_perms)
        current_perm = zones_obj_red_perms[i]
        zones_perm_cant = 0.0
        for row in row_data_arr
            row_zones = row |> collect
            are_zones_in_perm = check_zones_perm(row_zones, current_perm)
            if(are_zones_in_perm)
                zones_perm_cant += 1
            end
        end
        zones_obj_red_cants[i] = zones_perm_cant
    end
    best_zones_seq = zones_obj_red_perms[argmax(zones_obj_red_cants)]
    pushfirst!(best_zones_seq, "0")
    return best_zones_seq
end

function check_zones_perm(row_zones, current_perm)
    are_zones_in_perm = true
    current_index = 0
    for Z in current_perm
        zone_pos = findfirst(x->x==Z, row_zones)
        if(!isnothing(zone_pos) && zone_pos>current_index)
            current_index = zone_pos
        else
            are_zones_in_perm = false
            break
        end
    end
    return are_zones_in_perm
end

function evaluate_best_zones(routes)
    best_zones_seq_arr = Vector{Vector{String}}()
    for route_num in 1:length(routes)
        route_obj = routes[route_num]
        zones_obj = get_macrozones(zones(route_obj))
        zones_obj_red = zones_obj[1:length(zones_obj)-1]
        df_route_obj = get_reduced_dataframe(route_obj)
        #show(df_route_obj, allcols=true)
        #df_red_route_obj = df_route_obj[(df_route_obj.Z1 .== zones_obj_red[1]) .| (df_route_obj.Z1 .== zones_obj_red[2]), :]
        df_red_route_obj = filter_dataframe_with_zones(df_route_obj, zones_obj_red)
        #show(df_red_route_obj, allcols=true)
        zones_obj_red_perms = permutations(zones_obj_red) |> collect
        zones_obj_red_cants = zeros(length(zones_obj_red_perms))
        row_data_arr = eachrow(df_red_route_obj)
        for i in 1:length(zones_obj_red_perms)
            current_perm = zones_obj_red_perms[i]
            zones_perm_cant = 0.0
            for row in row_data_arr
                row_zones = row |> collect
                are_zones_in_perm = check_zones_perm(row_zones, current_perm)
                if(are_zones_in_perm)
                    zones_perm_cant += 1
                end
            end
            zones_obj_red_cants[i] = zones_perm_cant
        end
        best_zones_seq = zones_obj_red_perms[argmax(zones_obj_red_cants)]
        push!(best_zones_seq_arr, best_zones_seq)
        #println("BEST ZONES SEQ (ROUTE "*string(name(route_obj))*"): "*string(best_zones_seq))
    end
    return best_zones_seq_arr
end

function solve(alg::MinTimeMacro, route::AbstractRoute, times::TravelTimes; zones_order=get_zones_order(route), initial_stop=get_station(route))
    T = times.travel_times
    return solve(alg, route, T, zones_list=zones_order, initial_stop=initial_stop)
end

function solve(alg::MinTimeMacro, route::AbstractRoute, T::AbstractMatrix; zones_list=get_zones_order(route), initial_stop=get_station(route))
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
function _next_stop(::MinTimeMacro, subroute)
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
        #alg = MinTimeMacro()
        #alg = MinTime()
        # calcular la secuencia
        #sequences_rand[i] = solve(data_route_apply[i], data_times[i].travel_times)
        sequences_rand[i] = solve(alg, routes[i], data_times[i].travel_times)
        # aplicar el algoritmo
        #sequences_rand[i] = solve_mintimemacro(data_route_apply[i], data_times[i].travel_times) #initial_zone = "H24.2C") #get_stop(route,"LZ"))
        # obtener el score
        push!(scores_arr, get_score2(sequences_rand))
        #println(get_score2(sequences_rand[i]))
    end
    return sequences_rand, scores_arr
end

function make_compact_sequencing_process(alg, routes, data_times)
    #alg = MinTimeMacro()
    # calcular la secuencia
    return [solve(alg, routes[i], data_times[i].travel_times) for i in 1:length(routes)];
end