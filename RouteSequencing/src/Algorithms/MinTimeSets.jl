"""
    MinTimeSets <: AbstractSequencingAlgorithm

Algorithm that decides the sequence based on minimizing the time from the current
stop to the next stop.
"""

using GLPK

struct MinTimeSets{F} <: AbstractSequencingAlgorithm
    minfunc::F
end

MinTimeSets() = MinTimeSets(_min_length_path_exhaustive)

function solve(alg::MinTimeSets, route::AbstractRoute, times::TravelTimes; initial_stop=station(route))
    T = times.travel_times
    return solve(alg, route, T, initial_stop=initial_stop)
end

function solve(alg::MinTimeSets, route::AbstractRoute, T::AbstractMatrix;
               initial_stop=station(route))

    S = stops(route)
    num_stops = length(S)

    # salida
    route_id = name(route)
    zones_list = _min_zone(route)
    
    seq = Vector{BareStop}()
    out = Sequence(route_id, seq)
    
    x = initial_stop |> name
    push!(out, BareStop(x))
    println("ROUTE: "*string(name(route)))

    for Z in zones_list
        println("ZONES: "*Z)
        stops_per_zone = stops(route)[findall(x->zone(x)==Z, stops(route))]
        if length(stops_per_zone) == 1
            push!(out, BareStop(name(stops_per_zone[1])))
        else
            insert!(stops_per_zone, 1, stops(route)[findfirst(x->x.id==initial_stop.id, stops(route))])
            stops_str = [name(stops) for stops in stops_per_zone]
            time_table = T[stops_str, stops_str]

            _, result = _min_length_path_mip(time_table, stops_per_zone, optimizer=GLPK.Optimizer)  # Esto se puede optimizar eligiendo entre MIP y Exhaustive

            for r in result[2:length(result)]
                push!(out, BareStop(name(r)))
            end
            initial_stop = stops(route)[findfirst(x->x.id==result[length(result)].id, stops(route))]
        end
    end
    println(out)
    return out
end

function _min_length_path_exhaustive(T::TravelTimes, stops::Vector{<:AbstractStop})
    _min_length_path_exhaustive(T.travel_Times, stops)
end

# sdlve the minimum length hamiltonian path by brute-force search
# computes the minimum path of the given stops, starting from the first stop given
# in `stops` and finishing in any other other stops
function _min_length_path_exhaustive(T::AbstractMatrix, stops::Vector{<:AbstractStop})
    p0 = stops[1]

    # FIXME refactor
    stop_names = [name(s) for s in stops]
    M = view(T, stop_names, stop_names)

    perm = permutations(stops[2:end])

    num_stops = length(stops)

    minvalue = Inf
    local result
    for p in perm
        value = M[name(p0), name(p[1])]
        for i in 1:num_stops-2
            value += M[name(p[i]), name(p[i+1])]
        end

        if value < minvalue
            minvalue = value
            result = vcat(p0, p)
        end
    end

    return (minvalue, result)
end

# given an nxn matrix, return an (n+m)×(n+m) matrix with m new zero rows
# and m new zero columns
function _add_dimension(A::AbstractMatrix, m=1)
    n = size(A, 1)
    return vcat(hcat(A, zeros(n, m)), zeros(m, n+m))
end

# solve the minimum length Hamiltonian problem using a MIP formulation where
# the first and last stops are fixed
function _min_length_path_mip(T::AbstractMatrix, stops::Vector{<:AbstractStop};
                              optimizer)
    # FIXME refactor
    stop_names = [name(s) for s in stops]
    M = view(T, stop_names, stop_names)

    model = Model(optimizer)

    # number of stops
    N = length(stops)

    # binary decision variables: 1 if x_{st} stop s goes to stop t and 0 otherwise
    @variable(model, x[1:N, 1:N], Bin)

    # cost function
    @objective(model, Min, sum(x[s, t] * M[s, t] for s=1:N, t=1:N))

    for s = 1:N
        # no edge reaches the initial stop
        @constraint(model, x[s, 1] == 0)

        # no self-loops
        @constraint(model, x[s, s] == 0)

        # each stop arrives to one and only one stop (unless the last one)
        if s < N
            @constraint(model, sum(x[s, t] for t in 2:N) == 1)
        end
    end

    for t = 1:N
        # no edge after the final stop
        @constraint(model, x[N, t] == 0)

        # each stop receives one and only one edge (unless the first stop)
        if t > 1
           @constraint(model, sum(x[s, t] for s in 1:N) == 1)
        end
    end

    # no repetitions (not required in ppple)
    for  s= 1:N, t = 1:N
        @constraint(model, x[s, t] + x[t, s] <= 1)
    end

    optimize!(model)

    minvalue = JuMP.objective_value(model)

    Mpath = JuMP.getvalue.(x)

    result = similar(stops)
    result[1] = stops[1]
    for i in 1:N-1
        idx = findfirst(==(1,), view(Mpath, i, :))
        result[i+1] = stops[idx]
    end
    result[N] = stops[N]

    return (minvalue, result)
end

function _min_zone(R::AbstractRoute; initial_zone=zone(R.stops[2]), tol=1e-4)
    
    # split stops by zone
    points = point_cloud(R)
    zone_names = zones(R, remove_station=true)
    num_zones = length(zone_names)

    # overapproximate with convex polygons
    P = [overapproximate(ConvexHullArray(pi), HPolygon, 1e-3) for pi in points]
    
    # preallocate output array
    out = Vector{String}() # FIXME use Vector{Zone} ??
    push!(out, initial_zone)

    # compute cost matrix
    M = Matrix{Float64}(undef, num_zones, num_zones)
    for i in 1:num_zones
        for j in 1:num_zones
            M[i, j] = hausdorff_distance(P[i], P[j], ε=tol)
        end
    end
    T = NamedArray(M, (zone_names, zone_names))
    
    # for each zone, find the next one that we didn't visit and that is closer
    x = initial_zone    #  convert(String, R.stops[2].zone_id)

    # iterate over zones
    while length(out) < num_zones
        s = sort(view(T, x, :))

        i = 2 # ignoro la posicion actual (tiempo = 0)
        y = s.dicts[1].keys[i]
        while (y ∈ out) # iterar hasta encontrar uno que no haya visitado
            i += 1
            i > length(s.dicts[1]) && return out
            y = s.dicts[1].keys[i]
        end
        push!(out, y)
        x = y
    end 
    
    return out
end
