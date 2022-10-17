# instantiate project
import Pkg
Pkg.activate(@__DIR__)
Pkg.resolve()
Pkg.instantiate()

# install dependencies
#Pkg.add("JSON")
#Pkg.add("CSV")
#Pkg.add("NamedArrays")
#Pkg.add("JuMP")
# ...

include(joinpath("./RouteSequencing", "src", "RouteSequencing.jl"))

function run_all(BASE_DIR)

    # load travel times data
    TIMES_DIR = joinpath(BASE_DIR, "data", "model_apply_inputs", "new_travel_times.json")
    data_times = RouteSequencing.parse_travel_times(file=TIMES_DIR)

    # load route data
    ROUTE_DIR = joinpath(BASE_DIR, "data", "model_apply_inputs", "new_route_data.json")
    data_route = RouteSequencing.parse_route_data(file=ROUTE_DIR)
    num_routes = length(data_route)

    # choose sequencing algorithm
    alg = RouteSequencing.MinTime()

    # compute sequences
    sequences = [RouteSequencing.solve(alg, data_route[i], data_times[i]) for i in 1:num_routes]

    # write result to json file
    OUT_DIR = joinpath(BASE_DIR, "data", "model_apply_outputs", "proposed_sequences.json")
    RouteSequencing.write_json_sequence(sequences; path=OUT_DIR)
end
