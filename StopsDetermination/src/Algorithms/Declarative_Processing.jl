"""
    minizinc_preprocessing

Algorithm that prepares the stops determination data model using the raw GPS registries 
in order to be used in MiniZinc application through declarative processing.

Inputs: 
* filenames: List of Transports GPS registries.
* is_time_american: Bool regarding time format.

[Optional]
* params: Relevant parameters to configure Determination process.

Output:
* File with .dzn format compatible with MiniZinc containing data model specification.

--------------------------------------------------------------------------------------

    minizinc_postprocessing

Algorithm that further processes the stops determination output collected from MiniZinc 
optimization solvers, in order to obtain and define the corresponding stops.

Inputs: 
* opt_solver_output: Binary vector that selects registries associated to stops.
* filenames: List of Transports GPS registries.
* is_time_american: Bool regarding time format.
* model_no: Model number related to the current registires group being processed.

Outputs:
* dStops: Determined Stops (i.e. Filtered Registries).
* stops_quant: Stops Ammount used on registries.
"""

using CSV, DataFrames, Dates

struct DStops
    name::Vector{String}
    lat::Vector{Float64}
    lon::Vector{Float64}
    measures::Vector{Int64}
    planned_time::Vector{Float64}
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

function process_event(event_camion)
    event_id = "O"
    if(event_camion == 'Encendido')
        event_id = "E"
    else
        if(event_camion == 'Apagado')
            event_id = "A"
        end
    end
    return event_id
end

# Event Vector Formulation formatted as custom String
function get_events_string(events_vector)
    vector_str = "["
    for i in 1:(length(events_vector)-1)
        vector_str = vector_str * process_event(events_vector[i]) * ", "
    end
    return vector_str*process_event(events_vector[length(events_vector)])*"]"
end

# Generic Vector Formulation formatted as custom String
function get_vector_string(vector)
    vector_str = "["
    for i in 1:(length(vector)-1)
        vector_str = vector_str * string(vector[i]) * ", "
    end
    return vector_str*string(vector[length(vector)])*"]"
end

# Model Content obtention and Data recollection
function get_model_content(df_camion_filtered, limits, is_time_american, dParams)
    date_general = df_camion_filtered[!,:Fecha]
    time_general = df_camion_filtered[!,:Hora]
    
    time_camion = [get_duration([time_general[i],time_general[i+1]],[date_general[i],date_general[i+1]],is_time_american) for i in 1:length(time_general)-1]
    lat_camion = df_camion_filtered[!,:Lat]
    lon_camion = df_camion_filtered[!,:Lon]
    vel_camion = df_camion_filtered[!,:"Vel."]
    km_camion = df_camion_filtered[!,:Km]
    event_camion = df_camion_filtered[!,:Evento]

    model_str = " nreg = "*string(length(lat_camion))*"\n nparams = "*string(length(dParams))*"\n minstops = "*string(limits[1])*"\n maxstops = "*string(limits[2])*"\n time = "*get_vector_string(time_camion)*"\n lat = "*get_vector_string(lat_camion)*"\n lon = "*get_vector_string(lon_camion)*"\n vel = "*get_vector_string(vel_camion)*"\n km = "*get_vector_string(km_camion)*"\n events = "*string(get_events_string(event_camion))*"\n dParams = "*get_vector_string(dParams)
    return model_str
end

# Auxiliary function that saves a KML file
function save_dzn_file(filepath, str_to_save)
    open(filepath, "w") do dzn_file
        write(dzn_file, str_to_save)
    end
end

function minizinc_preprocessing(filenames::Vector{String}, limits::Vector{Int64}, is_time_american::Bool; params=dParams)

    model_count = 1
    for filepath in filenames
        df_camion = DataFrame(CSV.File(filepath))
        df_unique = unique(df_camion, 1)
        sort!(df_unique, :Fecha)
        unique_dates = df_unique[!,:Fecha]

        for date in unique_dates
            df_camion_filtered = filter(:Fecha => Fecha -> Fecha == date, df_camion)
            eachcol(df_camion_filtered)
            model_str = get_model_content(df_camion_filtered, limits, is_time_american, params)
            save_dzn_file("./StopsDeterm_MiniZinc/data/stop_determ_"*string(model_count)*".dzn", model_str)
            model_count = model_count + 1
        end
    end
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

# Translation of Solver Output to Determined Stops
function load_stops!(dStops, df_camion_filtered, opt_solver_output, is_time_american)
    date_general = df_camion_filtered[!,:Fecha]
    time_general = df_camion_filtered[!,:Hora]
    
    time_camion = [get_duration([time_general[i],time_general[i+1]],[date_general[i],date_general[i+1]],is_time_american) for i in 1:length(time_general)-1]
    lat_camion = df_camion_filtered[!,:Lat]
    lon_camion = df_camion_filtered[!,:Lon]
    vel_camion = df_camion_filtered[!,:"Vel."]
    km_camion = df_camion_filtered[!,:Km]

    for i in 1:length(lat_camion)
        if(opt_solver_output[i] == 1)
            push!(dStops.lat, lat_camion[i])
            push!(dStops.lon, lon_camion[i])
            push!(dStops.planned_time, get_duration(time_camion[i-1:i], date_camion[i-1:i], is_time_american))
        end
    end

end

function minizinc_postprocessing(opt_solver_output::Vector{Int64}, filenames::Vector{String}, is_time_american::Bool, model_no::Int)

    dStops = DStops(Vector{String}(),Vector{Float64}(), Vector{Float64}(), Vector{Int64}(), Vector{Float64}())
    stops_quant = length([opt_solver_output[i]==1 for i in 1:length(opt_solver_output)])
    model_count = 1
    for filepath in filenames
        df_camion = DataFrame(CSV.File(filepath))
        df_unique = unique(df_camion, 1)
        sort!(df_unique, :Fecha)
        unique_dates = df_unique[!,:Fecha]

        for date in unique_dates
            df_camion_filtered = filter(:Fecha => Fecha -> Fecha == date, df_camion)
            eachcol(df_camion_filtered)
            if(model_count == model_no)
                load_stops!(dStops, df_camion_filtered, opt_solver_output, is_time_american)
                break
            else
                model_count = model_count + 1
            end
        end
        model_count>model_no & break
    end
    push!(dStops.names, generate_stops_names(length(dStops.lat))
    return dStops, stops_quant
end