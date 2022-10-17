"""

    rstudio_preprocessing

Algorithm that prepares the zone identification data model using the zone centers 
coordinates and zone identifiers from specified training routes for posterior 
Principal Component Analysis processing in R Studio.

Inputs: 
* routes: List of Routes (which contain determined lists of stops and zones).
* filename: Complete Address of CSV File with relevant parameters for Identification.

[Optional]
* params: Relevant parameters to configure Identification process.

Output:
* File with .csv format compatible for R Studio containing data model specification.

--------------------------------------------------------------------------------------

    rstudio_postprocessing

Algorithm that further processes the zones identification output collected from R Studio 
Principal Component Analysis, in order to obtain and define the ZonesVar structure for
application with different approaches.

Inputs: 
* pca_filename: Complete Address of CSV File where R Studio stored the PCA results.
* zones_quantity: Amount of processed zones. [Default: length(zoneids)]

Outputs:
* zvar: ZoneVar structure needed for Zone Identification procedures.

"""

using CSV, DataFrames, Dates

struct ZonesVar
    v_alpha_pos::Vector{Float64}
    v_alpha_neg::Vector{Float64}
    v_beta_pos::Vector{Float64}
    v_beta_neg::Vector{Float64}
    v_gamma_pos::Vector{Float64}
    v_gamma_neg::Vector{Float64}
    v_delta_pos::Vector{Float64}
    v_delta_neg::Vector{Float64}
    v_count::Vector{Float64}
end

ZonesVar() = ZonesVar(zeros(2,1), zeros(2,1), zeros(2,1), zeros(2,1), zeros(2,1), zeros(2,1), zeros(2,1), zeros(2,1), zeros(8,1))

function normalize_vector(zv::ZonesVar)
    dims = 1:2
    for i in dims
        zv.v_alpha_pos[i] = zv.v_alpha_pos[i]/zv.v_count[1]
        zv.v_alpha_neg[i] = zv.v_alpha_neg[i]/zv.v_count[2]
        zv.v_beta_pos[i] = zv.v_beta_pos[i]/zv.v_count[3]
        zv.v_beta_neg[i] = zv.v_beta_neg[i]/zv.v_count[4]
        zv.v_gamma_pos[i] = zv.v_gamma_pos[i]/zv.v_count[5]
        zv.v_gamma_neg[i] = zv.v_gamma_neg[i]/zv.v_count[6]
        zv.v_delta_pos[i] = zv.v_delta_pos[i]/zv.v_count[7]
        zv.v_delta_neg[i] = zv.v_delta_neg[i]/zv.v_count[8]
    end
    return zv
end

function value_in_between(value, min_val, max_val)
    return value >= min_val && value <= max_val
end

function stops_in_region(stops, lat_range, lng_range)
    return !isnothing(findfirst(x-> value_in_between(x.lat, lat_range[1], lat_range[2]) && value_in_between(x.lng, lng_range[1], lng_range[2]), stops))
end

function filter_routes(routes, city_selector)
    lat_range = [0.0, 0.0] # zeros(2,1)
    lng_range = [0.0, 0.0] # zeros(2,1)
    if(city_selector == 1)      # Austin
        lat_range = [30.134, 30.517]
        lng_range = [-97.951, -97.531]
    elseif(city_selector == 2)  # Boston
        lat_range = [42.235, 42.397]
        lng_range = [-71.199, -70.967]
    elseif(city_selector == 3)  # Chicago
        lat_range = [41.653, 42.015]
        lng_range = [-87.866, -87.487]
    elseif(city_selector == 4)  # Los Angeles
        lat_range = [33.748, 34.312]
        lng_range = [-118.639, -118.152]
    elseif(city_selector == 5)  # Seattle
        lat_range = [47.502, 47.732]
        lng_range = [-122.435, -122.244]
    else                        # Montevideo
        lat_range = [-34.922, -34.705]
        lng_range = [-56.433, -56.014]
    end
    return routes[findall(x->stops_in_region(x.stops, lat_range, lng_range), routes)]
end

function get_mean_dist(dist_lng, dist_lat)
    dist_lng_mean = 0
    dist_lat_mean = 0
    for i in 1:length(dist_lat)
        dist_lng_mean = dist_lng_mean + dist_lng[i]
        dist_lat_mean = dist_lat_mean + dist_lat[i]
    end
    return dist_lng_mean/length(dist_lng), dist_lat_mean/length(dist_lat)
end

function get_mid_points(curr_stops)
    zone_lng_max = maximum([curr_stops.lng for i in 1:length(curr_stops)])
    zone_lng_min = minimum([curr_stops.lng for i in 1:length(curr_stops)]) 
    zone_lat_max = maximum([curr_stops.lat for i in 1:length(curr_stops)])
    zone_lat_min = minimum([curr_stops.lat for i in 1:length(curr_stops)])
    return [(zone_lng_max+zone_lng_min)/2, (zone_lat_max+zone_lat_min)/2]
end

function add_center_zone_coords(route, zone_centers, zoneids, centerid)
    route_lng_max = maximum([route.stops[i].lng for i in 1:length(route.stops)])
    route_lng_min = minimum([route.stops[i].lng for i in 1:length(route.stops)]) 
    route_length = route_lng_max - route_lng_min
    route_lat_max = maximum([route.stops[i].lat for i in 1:length(route.stops)])
    route_lat_min = minimum([route.stops[i].lat for i in 1:length(route.stops)]) 
    route_width = route_lat_max - route_lat_min
    for Z in unique(zones(route))
        if(Z!='0')
            curr_stops = route.stops[findall(x->zone(x)==Z,route.stops)]
            if(centerid == 1)
                curr_center = get_mean_dist([curr_stops[i].lng for i in 1:length(curr_stops)], [curr_stops[i].lat for i in 1:length(curr_stops)])
            elseif(centerid == 2)
                curr_center = get_mid_points(curr_stops) 
            elseif(centerid == 3)
                curr_center = (route_width/route_length)*get_mean_dist([curr_stops[i].lng for i in 1:length(curr_stops)], [curr_stops[i].lat for i in 1:length(curr_stops)])
            else
                curr_center = (route_width/route_length)*get_mid_points(curr_stops) 
            end
            push!(zones_centers, curr_center)
            push!(zoneids, zone_format([Z[0:1],Z[2:3],Z[4:5],Z[5:6]]))
        end
    end
    return zones_centers, zoneids
end

function zone_format(zone_id_arr)
    return string(zone_id_arr[1])*"-"*string(zone_id_arr[2])*"."*string(zone_id_arr[3])*"*"*string(zone_id_arr[4])
end

# Model Content obtention and Data recollection
function get_model_content(zones_centers, zoneids)
    num_zones = length(zoneids)
    model_str = "centerLat,centerLng,zone\n"
    for i in 1:num_zones
        model_str = model_str * zones_centers[i,1] * "," * zones_centers[i,2] * "," * zoneids[i]
    end
    return model_str
end

# Auxiliary function that saves a CSV file
function save_csv_file(filepath, str_to_save)
    open(filepath, "w") do csv_file
        write(csv_file, str_to_save)
    end
end

function rstudio_preprocessing(routes::Vector{Route}, filename::String)
    #zone_center_calculation = dParams[1] # {1, 2, 3, 4} = {Baricenter, MidPoint, Baricenter Scaled, MidPoint Scaled}
    #selected_city = dParams[2] # {0, 1, 2, 3, 4, 5, 6} = {All, Austin, Boston, Chicago, Los Angeles, Seattle, Montevideo}
    #routes_quantity = dParams[3] # {10, 20, 50}
    dParams = [1, 6, 50]
    filename = "./zonesData.csv"
    return solve(routes, filename; params=dParams)
end

function rstudio_preprocessing(routes::Vector{Route}, filename::String; params=dParams)
    route_num = 1
    params[2] != 0 && routes = filter_routes(routes, params[2])
    for R in routes
        center_coords, zone_identifiers = add_center_zone_coords(R, center_coords, zone_identifiers, params[1])
        route_num = route_num + 1
        route_num > params[3] && break
    end
    model_str = get_model_content(center_coords, zone_identifiers)
    save_csv_file(filename, model_str)
end

function rstudio_postprocessing(pca_filename::String, zones_quantity::Int64)
    zvar = ZonesVar()
    var_matrix = CSV.File(pca_filename) |> Tables.matrix
    zvar.v_alpha_pos = [var_matrix[1,1], var_matrix[2,1]]
    zvar.v_alpha_neg = [var_matrix[1,2], var_matrix[2,2]]
    zvar.v_beta_pos = [var_matrix[1,3], var_matrix[2,3]]
    zvar.v_beta_neg = [var_matrix[1,4], var_matrix[2,4]]
    zvar.v_gamma_pos = [var_matrix[1,5], var_matrix[2,5]]
    zvar.v_gamma_neg = [var_matrix[1,6], var_matrix[2,6]]
    zvar.v_delta_pos = [var_matrix[1,7], var_matrix[2,7]]
    zvar.v_delta_neg = [var_matrix[1,8], var_matrix[2,8]]
    zvar.count = zones_quantity
    return zvar
end