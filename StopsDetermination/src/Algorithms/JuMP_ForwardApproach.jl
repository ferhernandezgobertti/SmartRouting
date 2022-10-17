"""
    ForwardApproachJuMP <: AbstractDeterminationAlgorithm

Algorithm that determines the stops from raw GPS registries while looking forward 
in the chronologically order measurements using JuMP optimization solvers and libraries.

Inputs: 
* filenames: List of Transports GPS registries.
* is_time_american: Bool regarding time format.

[Optional]
* params: Relevant parameters to configure Determination process.

Outputs:
* dStops: Determined Stops (i.e. Filtered Registries).
* stops_quant_arr: List of Stops Ammount per day used on registries.
"""

using CSV, DataFrames, Dates
using JuMP, GLPK, CPLEX, Gurobi

struct ForwardApproachJuMP <: AbstractDeterminationAlgorithm end

struct DStops
    name::Vector{String}
    lat::Vector{Float64}
    lon::Vector{Float64}
    measures::Vector{Int64}
    planned_time::Vector{Float64}
    stop_type::Bool # 0 if Station, 1 if Dropoff
end 

struct Transport
    date_camion::Vector{String}
    time_camion::Vector{String}
    lat_camion::Vector{String} 
    lon_camion::Vector{String}
    km_camion::Vector{String}
    vel_camion::Vector{String}
    event_camion::Vector{String}
end

# Time duration between contiguous registries
function get_duration(timec, datec, is_time_american)
    date_format = "dd/mm/yyyy HH:MM:SS"
    if(is_time_american)
        date_format = "mm/dd/yyyy HH:MM:SS"
    end
    time1 = Dates.DateTime(datec[2]*" "*string(timec[2]), date_format)
    time2 = Dates.DateTime(datec[1]*" "*string(timec[1]), date_format)
    return (abs((time2-time1).value)/1000)
end

# Count of events identified with 'Encendido' or 'Apagado'
function get_count(events_list, str_to_count)
    return count(x->x==str_to_count, events_list)
end

# Model Definition and Application of ForwardApproach (Causal) for Stops Determination
function applicate_model_forward(transport, dParams; solver=OptSolver)

    model = Model(OptSolver.Optimizer)
    n = length(transport.lat_camion)
    @variable(model, x[1:n], Bin)

    # Regulatory Constraints
    @constraint(model, sum(x)<=120)
    @constraint(model, sum(x)>=80)

    for i in 1:n-dParams[1]
        # Events Constraint
        @constraint(model, (get_count(transport.event_camion[1:i],'Encendido')*x[i] - get_count(transport.event_camion[1:i],'Apagado')*x[i]) <= dParams[6])
        # Distance Constraint
        @constraint(model, sum([sqrt((transport.lat_camion[i+j-1]*x[i+j-1]-transport.lat_camion[i+j]*x[i+j])^2 + (transport.lon_camion[i+j-1]*x[i+j-1]-transport.lon_camion[i+j]*x[i+j])^2) for j in 1:dParams[1]]) < dParams[5])
        # Temporal Constraint
        @constraint(model, sum([get_duration(transport.time_camion[i:i+1]*x[i], transport.date_camion[i:i+1]*x[i], is_time_american) for i in 1:length(transport.time_camion)]) < dParams[4])
        # Velocity Constraint
        @constraint(model, sum([transport.vel_camion[i+j-1]*x[i+j-1]-transport.vel_camion[i+j]*x[i+j] for j in 1:dParams[1]]) < dParams[3])
        # Kilometer Constraint
        @constraint(model, sum([transport.km_camion[i+j-1]*x[i+j-1]-transport.km_camion[i+j]*x[i+j] for j in 1:dParams[1]]) < dParams[2]) 
    end

    cost = sum(x)
    @objective(model, Min, cost)
    optimize!(model) # @time optimize!(model); @allocated optimize!(model)

    return value.(x)
end

# Main function for transport examination and stops determination
function examine_transport!(dStops, df_camion_filtered, stops_quant_arr, dParams)
    
    date_camion = df_camion_filtered[!,:Fecha]
    time_camion = df_camion_filtered[!,:Hora]
    lat_camion = df_camion_filtered[!,:Lat]
    lon_camion = df_camion_filtered[!,:Lon]

    vel_camion = df_camion_filtered[!,:"Vel."]
    event_camion = df_camion_filtered[!,:Evento]
    km_camion = df_camion_filtered[!,:Km]

    transport = Transport(date_camion, time_camion, lat_camion, lon_camion, km_camion, vel_camion, event_camion)
    opt_solver_output = applicate_model_forward(transport, dParams, solver=Gurobi)  # Change solver as needed (GLPK, CPLEX, Gurobi)

    for i in 1:length(lat_camion)
        if(opt_solver_output[i] == 1)
            push!(dStops.lat, lat_camion[i])
            push!(dStops.lon, lon_camion[i])
            push!(dStops.planned_time, get_duration(time_camion[i-1:i], date_camion[i-1:i], is_time_american))
        end
    end    
    dStops[1].stop_type = true
end

# Stops identification, for completion. Irrelevant for zones prediction and stops routing procedures.
function generate_stops_names(stops_ammount)
    
    characters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    comb = combinations(characters,2)
    name_pos = Vector{String}()
    for l in comb
        push!(name_pos, l[1]*l[2])
    end
    shuffle!(name_pos)
    return name_pos[1:stops_ammount]
end

function solve(alg::ForwardApproachJuMP, filenames::Vector{String}, is_time_american::Bool)
    #stops_quant_measure = dParams[1] # 3 stops
    #kilometer_measure = dParams[2] # 1 cuadra (0.1km)
    #stops_vel_measure = dParams[3] # 10 km/h
    #close_time_measure = dParams[4] # 120 sec
    #close_dist_measure = dParams[5] # 0.1 coord
    dParams = [4,0.1,10,60,0.0001,1]
    return solve(alg, filenames, is_time_american; params=dParams)
end

function solve(alg::ForwardApproachJuMP, filenames::Vector{String}, is_time_american::Bool;
               params=dParams)

    dStops = DStops(Vector{String}(),Vector{Float64}(), Vector{Float64}(), Vector{Int64}(), Vector{Float64}(), false)
    stops_quant_arr = Vector{Int64}()
    
    for filepath in filenames
        df_camion = DataFrame(CSV.File(filepath))
        df_unique = unique(df_camion, 1)
        sort!(df_unique, :Fecha)
        unique_dates = df_unique[!,:Fecha]
        
        for date in unique_dates
            df_camion_filtered = filter(:Fecha => Fecha -> Fecha == date, df_camion)
            eachcol(df_camion_filtered)
            examine_transport!(dStops, df_camion_filtered, stops_quant_arr, dParams)
        end
    end
    push!(dStops.names, generate_stops_names(length(dStops.lat))
    return dStops, stops_quant_arr
end