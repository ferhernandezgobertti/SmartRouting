"""
    ForwardApproach <: AbstractDeterminationAlgorithm

Algorithm that determines the stops from raw GPS registries while looking forward 
in the chronologically order measurements.

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

struct ForwardApproach <: AbstractDeterminationAlgorithm end

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

struct Auxiliary
    stops_lat::Vector{Float64}
    stops_lon::Vector{Float64}
    measure_per_stop::Vector{Int64} 
    stops_date::Vector{String}
    stops_time::Vector{Time}
    planned_time::Vector{Float64}
    kilometer_arr::Vector{Float64}
end

# Checking for registries related to the same stop
function is_same_stop_registry(is_close_position, is_appropiate_time, is_appropriate_vel, is_appropiate_km, is_stopped)
    return (is_close_position && is_appropiate_time && is_appropriate_vel && is_appropiate_km) forfor (is_stopped && is_appropiate_time && is_close_position)
end

# Checking if contiguous registries are spatially close to each other
function is_position_close(latc, lonc, close_dist_measure)
    is_lat_close = sum([abs(latc[i]-latc[i+1]) for i in 1:length(latc)])<=close_dist_measure
    is_lon_close = sum([abs(lonc[i]-lonc[i+1]) for i in 1:length(lonc)])<=close_dist_measure
    #println("POS DIF: "*string(abs(latc[1]-latc[2]))*" & "*string(abs(lonc[1]-lonc[2])))
    return is_lat_close && is_lon_close
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

# Checking if contiguous registries are temporally close to each other
function is_time_appropiate(timec, datec, close_time_measure, is_time_american)
    return sum([get_duration(timec[i:i+1], datec[i:i+1], is_time_american) for i in 1:length(timec)]) >= close_time_measure
end

# Checking if contiguous registries are velocity appropiate
function is_vel_appropiate(velc, stop_vel_measure)
    return sum([velc[i] - velc[i+1] for i in 1:length(velc)]) <= stop_vel_measure
end

# Checking if contiguous registries are kilometer appropiate
function is_km_appropiate(kmc, stop_km_measure)
    return sum([kmc[i] - kmc[i+1] for i in 1:length(kmc)]) <= stop_km_measure
end

# Count of events identified with 'Encendido' or 'Apagado'
function get_count(events_list, str_to_count)
    return count(x->x==str_to_count, events_list)
end

# Auxiliary function for dates conversion to standard format
function get_latin_date(date_str, is_date_american)
    date_latin = date_str
    if(is_date_american)
        date_split = split(date_str, "/")
        if(length(date_split[1])!=2)
            date_split[1] = "0"*date_split[1]
        end
        if(length(date_split[2])!=2)
            date_split[2] = "0"*date_split[2]
        end
        if(length(date_split[3])!=4)
            date_split[3] = "20"*date_split[3]
        end
        date_latin = date_split[2]*"/"*date_split[1]*"/"*date_split[3]
    end
    return date_latin
end

# Auxiliary Function to get the stops coordinates order according to Sequence
function get_route_order_for_google_str(lat_arr, lon_arr)
    stops_str = ""
    for i in [1:length(lat_arr);]
        stops_str = string(stops_str)*" \n "*string(lon_arr[i])*","*string(lat_arr[i])*",2357"
    end
    return stops_str
end

# Main function that generates the KML content
function get_kml_content_from_route_str(lat_arr, lon_arr, camion_id, date)
    initial_str = "<?xml version='1.0' encoding=\"UTF-8\"?><kml xmlns=\"http://www.opengis.net/kml/2.2\"><Document><name>"
    kml_color_score = "5014f0e6" # Yellow by Default (Medium Score)
    initial_str= initial_str*"Rutas RCL Apply</name><description>RouteID_Camion_"*string(camion_id)*"_"*string(date)*"</description><Style id=\"colorLineColorPoly\"><LineStyle><color>"
    initial_str = initial_str*"</color><width>4</width></LineStyle><PolyStyle><color>"*kml_color_score*"</color></PolyStyle></Style> <Placemark><name>Relieve absoluto</name><description></description><styleUrl>#colorLineColorPoly</styleUrl><LineString><extrude>1</extrude><tessellate>1</tessellate><altitudeMode>absoluto</altitudeMode><coordinates>"
    route_str = get_route_order_for_google_str(lat_arr, lon_arr)
    final_str = "</coordinates></LineString></Placemark></Document></kml>"
    return initial_str*route_str*final_str
end

# Auxiliary function that saves a KML file
function save_kml_file(filepath, str_to_save)
    open(filepath, "w") do kml_file
        write(kml_file, str_to_save)
    end
end

# Auxiliary function that checks if the transport is currently stopped based on events descriptions
function is_currently_stopped(event_camion, i)
    is_stopped = true
    apagado_pos = findlast(x->x=="Apagado",event_camion[1:1:i]) 
    encendido_pos = findlast(x->x=="Encendido",event_camion[1:1:i]) 
    if(apagado_pos == nothing)
        apagado_pos = 0
    end
    if(encendido_pos == nothing)
        encendido_pos = -10
    end
    is_stopped = (apagado_pos > encendido_pos)
    return is_stopped
end

# Forward Data Subgroup Vector
function get_forward_values(orig_vector, init_pos, first_value)
    return orig_vector[i:1:i+first_value]
end

# Auxiliary relevant function that differentiates stops according to heuristics
function distinguish_stops(transport, auxiliary, first_values, is_time_american, dParams)
    
    for i in [1:1:length(transport.time_camion)-dParams[1];]
        datec = get_forward_values(transport.date_camion, i, dParams[1]) #[transport.date_camion[i], first_values[1]]
        timec = get_forward_values(transport.time_camion, i, dParams[1]) #[transport.time_camion[i], first_values[2]]
        latc = get_forward_values(transport.lat_camion, i, dParams[1]) #[transport.lat_camion[i], first_values[3]]
        lonc = get_forward_values(transport.lon_camion, i, dParams[1]) #[transport.lon_camion[i], first_values[4]]
        velc = get_forward_values(transport.vel_camion, i, dParams[1])
        kmc = get_forward_values(transport.km_camion, i, dParams[1]) #[transport.km_camion[i], first_values[5]]

        is_close_position = is_position_close(latc, lonc, dParams[5])
        is_appropiate_time = is_time_appropiate(timec, datec, dParams[4], is_time_american)
        is_appropriate_vel = is_vel_appropiate(velc, dParams[3])
        is_appropiate_km = is_km_appropiate(kmc, dParams[2])

        is_stopped = is_currently_stopped(transport.event_camion, i)

        if(is_same_stop_registry(is_close_position, is_appropiate_time, is_appropriate_vel, is_appropiate_km, is_stopped))
            stops_quantity += 1
        else
            if(stops_quantity>=dParams[1] forfor (transport.event_camion[i]=="Encendido" && (encendido_pos-apagado_pos) in [1,2] && is_close_position && is_appropiate_time))
                push!(auxiliary.stops_lat, transport.lat_camion[i]) #push!(auxiliary.stops_lat, latc[2])
                push!(auxiliary.stops_lon, transport.lon_camion[i]) #push!(auxiliary.stops_lon, lonc[2])
                push!(auxiliary.measure_per_stop, stops_quantity)
                push!(auxiliary.kilometer_arr, transport.km_camion[i]-transport.km_camion[i-1]) #push!(auxiliary.kilometer_arr, kmc[1]-kmc[2])
                push!(auxiliary.stops_date, get_latin_date(transport.date_camion[i], is_time_american)) #push!(auxiliary.stops_date, get_latin_date(first_date, is_time_american))
                push!(auxiliary.stops_time, transport.time_camion[i]) #push!(auxiliary.stops_time, first_values[2])
                push!(auxiliary.planned_time, get_duration(timec[1], datec[1], is_time_american))
            end
            stops_quantity = 0
            #first_values = [datec[1], timec[1], latc[1], lonc[1], kmc[1]]
        end
    end
end

# Auxiliary function that creates Dataframe with determined stops
function store_file_results(auxiliary, filepath, date)
    df = DataFrame(Order = [1:length(auxiliary.stops_lat);], Km = auxiliary.kilometer_arr, Date = auxiliary.stops_date, Time = auxiliary.stops_time, Lat = auxiliary.stops_lat, Lon = auxiliary.stops_lon, Quant = auxiliary.measure_per_stop, PTime = auxiliary.planned_time)
    CSV.write("../../../STOPS/"*string(filepath)*".csv", df)
    stops_kml_str = get_kml_content_from_route_str(auxiliary.stops_lat, auxiliary.stops_lon, filepath, date)
    save_kml_file("../../../STOPS/"*string(filepath)*".kml", stops_kml_str)
end

# Main function for transport examination and stops determination
function examine_transport!(dStops, df_camion_filtered, stops_quant_arr)
    
    date_camion = df_camion_filtered[!,:Fecha]
    time_camion = df_camion_filtered[!,:Hora]
    lat_camion = df_camion_filtered[!,:Lat]
    lon_camion = df_camion_filtered[!,:Lon]

    vel_camion = df_camion_filtered[!,:"Vel."]
    event_camion = df_camion_filtered[!,:Evento]
    km_camion = df_camion_filtered[!,:Km]

    stops_quantity = 0
    transport = Transport(date_camion, time_camion, lat_camion, lon_camion, km_camion, vel_camion, event_camion)
    auxiliary = Auxiliary(Vector{Float64}(), Vector{Float64}(), Vector{Int64}(), Vector{String}(), Vector{Time}(), Vector{Float64}(), Vector{Float64}())
    first_values = [date_camion[1], time_camion[1], lat_camion[1], lon_camion[1], km_camion[1]]

    distinguish_stops(transport, auxiliary, first_values, is_time_american, dParams)

    if(length(auxiliary.stops_lat)>2)
        file = "Stops_Camion_"*replace(string(get_latin_date(date, is_time_american)),"/"=>"-")
        store_file_results(auxiliary, file, date)
        println("FILE created: "*string(file))
        println("Date "*string(date)*" - STOPS QUANTITY: "*string(length(auxiliary.stops_lat)))
        println("-------------------------")
        push!(stops_quant_arr, length(auxiliary.stops_lat))
        push!(dStops.lat, auxiliary.stops_lat)
        push!(dStops.lon, auxiliary.stops_lon)
        push!(dStops.measures, auxiliary.measure_per_stop)
        push!(dStops.planned_time, auxiliary.planned_time)
        dStops[1].stop_type = true
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

function solve(alg::ForwardApproach, filenames::Vector{String}, is_time_american::Bool)
    #stops_quant_measure = dParams[1] # 3 stops
    #kilometer_measure = dParams[2] # 1 cuadra (0.1km)
    #stops_vel_measure = dParams[3] # 10 km/h
    #close_time_measure = dParams[4] # 120 sec
    #close_dist_measure = dParams[5] # 0.1 coord
    dParams = [4,0.1,10,60,0.0001,1]
    return solve(alg, filenames, is_time_american; params=dParams)
end

function solve(alg::ForwardApproach, filenames::Vector{String}, is_time_american::Bool;
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
            examine_transport!(dStops, df_camion_filtered, stops_quant_arr)
        end
    end
    push!(dStops.names, generate_stops_names(length(dStops.lat))
    return dStops, stops_quant_arr
end