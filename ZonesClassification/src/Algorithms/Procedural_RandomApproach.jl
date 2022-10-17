# 1. Procedural - OK

"""
    ProRandomApproach <: AbstractClassificationAlgorithm

    Algorithm that clasifies the already defined stops according to documented
    rules considering already segmented and identified zones per route while using
    customized procedural metholodogies and geographical information from previously
    established zones relative to the training routes. Random route selection for 
    training and validation/testing datasets.
    
    Inputs: 
    * routes: List of Training and Validation Transports Routes (unsequenced).
    
    or 
    * routes: List of Training and Validation/Testing Transports Routes (unsequenced).
    * zclust: List of Zones Clusterings from Training Routes Processing.
    
    [Optional]
    * params: Relevant parameters to configure Classification process.
    
    Outputs: trainroutes, valroutes, zclassif, zscores, fig_arr
    * trainroutes: List of Training routes (selected through random algorithm).
    * valroutes: List of updated Validation routes (selected through random algorithm).
    * zclassif: List of configured Classification attributes (e.g., predicted zones).
    * zscores: Scores pero Zone Layer result of current Classification procedure.
    * fig_arr: List of Figures of Zones and Associated Stops for Plotting.

"""

using CSV, DataFrames, Dates, RouteSequencing
using LazySets, Plots, Colors, ReverseGeocode, Random
using StatsBase, StaticArrays, NamedArrays, LinearAlgebra

struct ProRandomApproach <: AbstractClassificationAlgorithm end

struct ZonesAnalysis
    zone_pred::String
    polygon_count::Float64 
    min_dist_sides::Float64 
    min_dist_vertex::Float64 
    internal_count::Float64
end

ZonesAnalysis() = ZonesAnalysis("", 0.0, 1000.0, 1000.0, 0.0)

struct ZonesCluster{ST}
    id::String
    city::String
    polygon_arr::Vector{ST}
end

ZonesCluster(id::String,city::String,poly_arr) = ZonesCluster(id,city,poly_arr)

function is_same_zonecluster(zc1::ZonesCluster, zc2::ZonesCluster)
    is_equal_id = zc1.id == zc2.id
    is_equal_city = zc1.city == zc2.city
    return is_equal_id && is_equal_city
end

cities(ZC::ZonesCluster) = ZC.city

function is_equal_zone_clust(zc1, zc2)
    is_same_id = zc1.id == zc2.id
    is_same_city = zc1.city == zc2.city
    return is_same_id && is_same_city
end

struct ZonesDetails{ST}
    poly_count_arr::Vector{Float64}
    min_dist_sides_arr::Vector{Float64}
    min_dist_vertex_arr::Vector{Float64}
    zone_cluster_list::Vector{ZonesCluster}
end

ZonesDetails() = ZonesDetails(Vector{Float64}(),Vector{Float64}(),Vector{Float64}(),Vector{ZonesCluster}())

struct ZonesClassif{ST}
    zones_actual::Vector{String}
    zones_predicted::Vector{String}
    zones_unregistered::Vector{Stop}
    stops_used::Vector{Stop}
    zones_details::Vector{ZonesDetails}
end

ZonesClassif() = ZonesClassif(Vector{String}(),Vector{String}(),Vector{Stop}(),Vector{Stop}(),Vector{ZonesDetails}())

struct ZonesScores{ST}
    score_alpha::Int64
    score_beta::Int64
    score_gamma::Int64
    score_delta::Int64
    zscores.score_zones::Int64
    total_zones::Int64
end

ZonesScores(zones_quantity) = ZonesScores(0,0,0,0,0,zones_quantity)

function stops_2(R::AbstractRoute; remove_station=false)
    S = R.stops
    if remove_station
        return deleteat!(S, findall(x->x.id==station(R).id,S))
    else
        return S
    end
end

function split_by_zone_2(R::AbstractRoute; remove_station=true)

    S = stops_2(R, remove_station=remove_station)
    deleteat!(S, findall(x->zone(x)=="0",S))
    num_stops = length(S)

    zone_names = unique(zone.(S))
    num_zones = length(zone_names)

    # preallocate output vectors
    zone_idx = [Vector{Int}() for _ in 1:num_zones]
    stop_names = [Vector{String}() for _ in 1:num_zones]

    for i in 1:num_stops
        idx = findfirst(x -> x == zone(S[i]), zone_names)
        push!(zone_idx[idx], i)
        push!(stop_names[idx], name(S[i]))
    end
    return zone_idx, zone_names, stop_names
end

function split_by_zone_3(R::AbstractRoute; remove_station=true)

    S = stops_2(R, remove_station=remove_station)
    #deleteat!(S, findall(x->zone(x)=="0",S))
    num_stops = length(S)
    zone_names = unique(zone.(S))
    num_zones = length(zone_names)

    # preallocate output vectors
    zone_idx = [Vector{Int}() for _ in 1:num_zones]
    stop_names = [Vector{String}() for _ in 1:num_zones]

    for i in 1:num_stops
        idx = findfirst(x -> x == zone(S[i]), zone_names)
        push!(zone_idx[idx], i)
        push!(stop_names[idx], name(S[i]))
    end
    return zone_idx, zone_names, stop_names
end

function fill_zones_cluster(data_route_build)
    zones_cluster = Vector{ZonesCluster}()
    data_route_build_c = deepcopy(data_route_build)

    for R in data_route_build_c
        #R = interpolate_zones(R)
        #out = split_by_zone_3(R, remove_station=true) # OJO remove_station=true esta BIEN AHORA
        out = split_by_zone_2(R, remove_station=false)
        P = out[1]
        latitudes = [latitude.(R.stops[P[j]]) for j in 1:length(P)]
        longitudes = [longitude.(R.stops[P[j]]) for j in 1:length(P)]
        
        zones = Vector{String}()
        stop_ref = Vector{Stop}()
        for px in P
            push!(stop_ref, R.stops[px[1]])
            push!(zones, zone(R.stops[px[1]]))
        end
        deleteat!(stop_ref, findall(x->x=="0",zones))
        deleteat!(zones, findall(x->x=="0",zones))
        
        poly = [[[x, y] for (x, y) in zip(longitudes[i], latitudes[i])] for i in 1:length(longitudes)];
        poly = VPolygon.(poly);
        deleteat!(poly, length(poly)) # sacar la station
        
        println("LEN ZONES: "*string(length(zones)))
        println("LEN POLY: "*string(length(poly)))
        println("------------------------------")
        
        for i in 1:length(poly)
            id = zones[i]
            city = decode(gc, SA[latitude(stop_ref[i]), longitude(stop_ref[i])]).city
            poly_arr = Vector{VPolygon}()
            push!(poly_arr, poly[i])
            zc = ZonesCluster(id,city,poly_arr)
            
            repeated_pos = findall(x->is_same_zonecluster(x,zc),zones_cluster)
            if(length(repeated_pos)>0)
                #println("REPEATED_POS: "*string(repeated_pos))
                println("LENGTH ZC: "*string(length(zones_cluster)))
                println(typeof(zones_cluster[repeated_pos[1]]))
                push!(zones_cluster[repeated_pos[1]].polygon_arr, poly[i])
            else
                push!(zones_cluster, zc)
            end
        end
    end
    return zones_cluster
end

function get_euclidian_distance(p1, p2)
    dist = sqrt((p1[1]-p2[1])^2 + (p1[2]-p2[2])^2)
    #if(dist==0)
    #    dist = 100
    #end
    return dist
end

function get_manhattan_distance(p1, p2)
    dist = abs(p1[1]-p2[1]) + abs(p1[2]-p2[2])
    #if(dist==0)
    #    dist = 100
    #end
    return dist
end

function get_baricenter(polygon_arr)
    #vertices_count = length(polygon_vertices)
    v_count = 0
    vertices_lng = 0
    vertices_lat = 0
    for P in polygon_arr
        for V in P.vertices
            vertices_lng = vertices_lng + V[1]
            vertices_lat = vertices_lat + V[2]
            v_count = v_count + 1
        end
    end
    return vertices_lng/v_count, vertices_lat/v_count
end

function get_zones_per_city(zones_cluster)
    unique_cities = unique(cities.(zones_cluster))
    zones_cities = zeros(length(unique_cities))
    for ZC in zones_cluster
        zones_cities_pos = findfirst(x->x==ZC.city, unique_cities)
        zones_cities[zones_cities_pos] = zones_cities[zones_cities_pos]+1
    end
    return zones_cities
end

function get_mean_point(polygon_arr)
    #vertices_count = length(polygon_vertices)
    v_count = 0
    vertices_lng = 0
    vertices_lat = 0
    for P in polygon_arr
        for V in P.vertices
            vertices_lng = vertices_lng + V[1]
            vertices_lat = vertices_lat + V[2]
            v_count = v_count + 1
        end
    end
    return vertices_lng/v_count, vertices_lat/v_count
end

# Zid: alpha-beta.gammadelta
# alpha: Macromacro
# beta: Macro
# gamma: Micro
# delta: Minimicro

function get_mean_dist(dist_lng, dist_lat)
    dist_lng_mean = 0
    dist_lat_mean = 0
    for i in 1:length(dist_lat)
        dist_lng_mean = dist_lng_mean + dist_lng[i]
        dist_lat_mean = dist_lat_mean + dist_lat[i]
    end
    return dist_lng_mean/length(dist_lng), dist_lat_mean/length(dist_lat)
end

function is_stop_in_polygon(stop, vertex)
    measure_arr = Vector{Float64}()
    stop_lat = latitude(stop)
    stop_lon = longitude(stop)
    for i in 1:length(vertex)
        v1 = vertex[i]
        if(i==length(vertex))
            v2 = vertex[1]
        else
            v2 = vertex[i+1]
        end
        measure = (stop_lat-v1[2])*(v2[1]-v1[1])-(stop_lon-v1[1])*(v2[2]-v1[2])
        push!(measure_arr, measure)
    end
    return ((all(>=(0), measure_arr) || all(<=(0), measure_arr)) && length(measure_arr)>1), minimum(abs.(measure_arr))
end

function get_measure(stop, v2, v1, tolerance)
    stop_lat = latitude(stop)
    stop_lon = longitude(stop)
    new_v1 = v1 - ones(length(v1))*tolerance
    new_v2 = v2 - ones(length(v2))*tolerance
    return (stop_lat-new_v1[2])*(new_v2[1]-new_v1[1])-(stop_lon-new_v1[1])*(new_v2[2]-new_v1[2])
end

function is_stop_in_polygon_with_tolerance2(stop, vertex, tolerance)
    measure_arr_low = Vector{Float64}()
    measure_arr = Vector{Float64}()
    measure_arr_high = Vector{Float64}()
    for i in 1:length(vertex)
        v1 = vertex[i]
        if(i==length(vertex))
            v2 = vertex[1]
        else
            v2 = vertex[i+1]
        end
        push!(measure_arr_high, get_measure(stop, v2, v1, tolerance))
        push!(measure_arr, get_measure(stop, v2, v1, 0))
        push!(measure_arr_low, get_measure(stop, v2, v1, -tolerance))
    end
    same_sign_arr_high = (all(>=(0), measure_arr_high) || all(<=(0), measure_arr_high)) && length(measure_arr_high)>1
    same_sign_arr = (all(>=(0), measure_arr) || all(<=(0), measure_arr)) && length(measure_arr)>1
    same_sign_arr_low = (all(>=(0), measure_arr_low) || all(<=(0), measure_arr_low)) && length(measure_arr_low)>1
    length_high = length(measure_arr_high)>1
    length_normal = length(measure_arr)>1
    length_low = length(measure_arr_low)>1
    return same_sign_arr_high+same_sign_arr+same_sign_arr_low, length_high+length_normal+length_low, minimum(abs.(measure_arr))
end

function get_min_dist_vertex(stop, poly_vertices, distance_measure)
    dist_arr = Vector{Float64}()
    for vertex in poly_vertices
        stop_coords = [longitude(stop), latitude(stop)]
        if(distance_measure == 0)
            dist = get_manhattan_distance(vertex, stop_coords)
        else
            dist = get_euclidian_distance(vertex, stop_coords)
        end
        push!(dist_arr, dist)
    end
    return minimum(dist_arr)
end

function get_proj_mods(zv, dist_stop)
    proj_mods = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # zeros(8,1)
    proj_mods[1] = get_proj(dist_stop,zv.v_alpha_pos)
    proj_mods[2] = get_proj(dist_stop,zv.v_alpha_neg)
    proj_mods[3] = get_proj(dist_stop,zv.v_beta_pos)
    proj_mods[4] = get_proj(dist_stop,zv.v_beta_neg)
    proj_mods[5] = get_proj(dist_stop,zv.v_gamma_pos)
    proj_mods[6] = get_proj(dist_stop,zv.v_gamma_neg)
    proj_mods[7] = get_proj(dist_stop,zv.v_delta_pos)
    proj_mods[8] = get_proj(dist_stop,zv.v_delta_neg)
    return proj_mods
end

function predict_zone_with_tolerance2(stop, zone_cluster, params)
    zanalysis = ZonesAnalysis()
    zone_length_cant = Vector{String}()
    min_dist_length_cant = Vector{Float64}()
    stop_internal_arr = zeros(length(zone_cluster.polygon_arr))
    zone_polygon_bool = zeros(length(zone_cluster.polygon_arr))
    zone_polygon_arr = ones(length(zone_cluster.polygon_arr))*1000
    zone_vertices_arr = ones(length(zone_cluster.polygon_arr))*1000
    for i in 1:length(zone_cluster.polygon_arr)
        poly = zone_cluster.polygon_arr[i]
        poly_vertices = poly.vertices
        stop_internal_cant, length_cant, min_dist_side = is_stop_in_polygon_with_tolerance2(stop, poly_vertices, params[6])
        min_dist_vertex = get_min_dist_vertex(stop, poly_vertices, params[9])
        push!(zone_length_cant, zone_cluster.id)
        push!(min_dist_length_cant, min_dist_vertex)
        zone_polygon_arr[i] = min_dist_side
        zone_vertices_arr[i] = min_dist_vertex
        if(stop_internal_cant>0) #|| min_dist_side<=params[6] || min_dist_vertex<=params[6])
            zanalisis.zone_pred = zone_cluster.id
            stop_internal_arr[i] = stop_internal_cant
            zone_polygon_bool[i] = 1
            zone_polygon_arr[i] = min_dist_side # 1
            zone_vertices_arr[i] = min_dist_vertex
            #break
            #else
            #min_dist_vertex<=0.001 && zanalisis.zone_pred = zone_cluster.id
            #end
        end
    end
    zanalisis.polygon_count = length(findall(x->x==1,zone_polygon_bool))
    zanalisis.min_dist_sides = minimum(abs.(zone_polygon_arr))
    zanalisis.min_dist_vertex = minimum(abs.(zone_vertices_arr))
    zanalisis.internal_count = maximum(stop_internal_arr)
    #if(zanalisis.polygon_count==0 && zanalisis.min_dist_vertex<=0.001)
    #    zanalisis.zone_pred = zone_cluster.id
    #end
    #if(zanalisis.min_dist_vertex in min_dist_length_cant)
    #    zanalisis.zone_pred = zone_length_cant[findfirst(x->x==min_dist_vertex, min_dist_length_cant)] 
    #end
    return zanalysis #zone_pred, polygon_count, min_dist_sides, min_dist_vertex, internal_count 
end

function get_predicted_actual(routes, zones_cluster, params)
    route_count = 1
    zclassif = ZonesClassif()
    #zone_limit = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0] # No se usa
    for R in routes
        println("ROUTE: "*string(route_count))
        for S in R.stops
            zdetails = ZonesDetails()
            if(zone(S)!="0")
                city_exists = false
                stop_polygon_max = 0
                min_dist_sides_min = params[4]
                min_dist_vertex_min = params[5]
                zone_pred_max = ""
                push!(zclassif.stops_used, S)
                push!(zclassif.zones_actual, zone(S))
                stop_city = decode(gc, SA[latitude(S), longitude(S)]).city
                for ZC in zones_cluster
                    zone_pred = ""
                    if(stop_city==ZC.city)
                        city_exists = true
                        #zone_pred, polygon_count, min_dist_sides, min_dist_vertex, internal_count = predict_zone_with_tolerance2(S, ZC, params)
                        zanalisis = predict_zone_with_tolerance2(S, ZC, params)
                        push!(zdetails.zone_cluster_list, ZC)
                        push!(zdetails.min_dist_sides_arr, min_dist_sides)
                        push!(zdetails.min_dist_vertex_arr, zanalisis.min_dist_vertex)
                        stop_polygon_count = 0
                        zanalisis.polygon_count >= params[3] && stop_polygon_count = zanalisis.polygon_count
                        push!(zdetails.poly_count_arr, stop_polygon_count)
                        if( min_dist_vertex_min>=zanalisis.min_dist_vertex) #zanalisis.zone_pred!="" &&
                            min_dist_vertex_min = zanalisis.min_dist_vertex
                        end
                        if(zanalisis.zone_pred!="" && zanalisis.min_dist_vertex==0) # ==0
                            zone_pred_max = zanalisis.zone_pred
                            min_dist_vertex_min = zanalisis.min_dist_vertex
                            break
                        end
                        if(min_dist_sides_min>=zanalisis.min_dist_sides)
                            min_dist_sides_min = zanalisis.min_dist_sides
                        end
                        if(zanalisis.zone_pred!="" && zanalisis.min_dist_sides==0) # ==0
                            zone_pred_max = zanalisis.zone_pred
                            min_dist_sides_min = zanalisis.min_dist_sides
                            break
                        end
                        if(zanalisis.zone_pred!="" && stop_polygon_count>=stop_polygon_max)
                            zone_pred_max = zanalisis.zone_pred
                            stop_polygon_max = stop_polygon_count
                        end
                    end
                end
                if(isnothing(findfirst(x->x.id==zone(S)&&x.city==stop_city,zones_cluster)))
                    push!(zclassif.zones_unregistered, S)
                end
                #!city_exists && println("CIUDAD NO EXISTE")
                if(zone_pred_max=="")
                    #push!(zclassif.zones_predicted, "NaN")
                    #zone_pred_max = zone(S)
                    #zone_pred_max = zdetails.zone_cluster_list[argmin(zdetails.min_dist_vertex_arr)].id
                    closest_zone_cluster = zdetails.zone_cluster_list[argmin(zdetails.min_dist_vertex_arr)]
                    closest_zc_bari = get_baricenter(closest_zone_cluster.polygon_arr)
                    dist_stop = [S.lng - closest_zc_bari[1], S.lat - closest_zc_bari[2]]
                    proj_mods = get_proj_mods(dist_stop)
                    dif_mods = [proj_mods[5]-proj_mods[6], proj_mods[7]-proj_mods[8]]
                    rel_dif = dif_mods[1]/dif_mods[2]
                    closest_zone = split(closest_zone_cluster.id,".")[2]
                    closest_zone_mod = closest_zone
                    if(abs(rel_dif)>params[7]) # Change GAMMA
                        change_mod = (-1)^(dif_mods[1]>0)
                        closest_zone_mod = Char(Int(only(closest_zone[1:1]))+change_mod)*string(closest_zone[2:end])
                    end
                    if(abs(rel_dif)<params[8]) # Change DELTA
                        change_mod = (-1)^(dif_mods[2]>0)
                        closest_zone_mod = string(closest_zone[1:1])*Char(Int(only(closest_zone[2:2]))+change_mod)
                    end
                    zone_pred_max = split(closest_zone_cluster.id,".")[1]*"."*string(closest_zone_mod)
                end
                push!(zclassif.zones_predicted, zone_pred_max)
                if(zone_pred_max!=zone(S))
                    poly_count_act_pos = findfirst(x->x.id==zone(S)&&x.city==stop_city,zdetails.zone_cluster_list)
                    poly_count_act = "" #zdetails.poly_count_arr[]
                    min_dist_sides_act = params[4]
                    min_dist_vertex_act = params[5]
                    if(isnothing(poly_count_act_pos))
                        sum = 1
                        #println("ZONE_CLUSTER_LIST: "*string(zdetails.zone_cluster_list))
                    else
                        poly_count_act = zdetails.poly_count_arr[poly_count_act_pos]
                        min_dist_sides_act = zdetails.min_dist_sides_arr[poly_count_act_pos]
                        min_dist_vertex_act = zdetails.min_dist_vertex_arr[poly_count_act_pos]
                    end
                end
            end
            S.zone_id = Zone(zone_pred)
            push!(zclassif.zones_details, zdetails)
        end
        route_count = route_count+1
    end
    return routes, zclassif
end

function get_plot_route_zones(zones_cluster, chosen_city)
    zones_cluster_city = zones_cluster[findall(x->x.city==chosen_city, zones_cluster)]
    zones_used = Vector{String}()
    sort!(zones_cluster_city, by=x->x.id)
    cols = distinguishable_colors(length(zones_cluster_city), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    fig_zones = plot(fmt = :png)
    limit = 500
    zones_limit = (limit<=length(zones_cluster_city) ? limit : length(zones_cluster_city))
    for i in 1:zones_limit
        #push!(zones_used, zones_cluster_city[i].id)
        zones_city = zones_cluster_city[i]
        first_vertices = zones_city.polygon_arr[1].vertices[1]
        #first_vertices_lng = zones_city.polygon_arr[1].vertices[1][1]
        #first_vertices_lat = zones_city.polygon_arr[1].vertices[1][2]
        first_vertices_lng, first_vertices_lat = get_mean_point(zones_city.polygon_arr)
        #plot!(fig_zones2, zones_city.polygon_arr[1], xlims=x_lims, ylims=y_lims, color=cols[i], annotation=(first_vertices_lng,first_vertices_lat, Plots.text(zones_city.id,12)), labels=zones_city.id, legends=(0.1,0.9))
        #[plot!(fig_zones2, zones_city.polygon_arr[j], xlims=x_lims, ylims=y_lims, color=cols[i]) for j in 2:length(zones_city.polygon_arr)]
        plot!(fig_zones, zones_city.polygon_arr[1], color=cols[i], annotation=(first_vertices_lng,first_vertices_lat, Plots.text(zones_city.id,8)))
        [plot!(fig_zones, zones_city.polygon_arr[j], color=cols[i]) for j in 2:length(zones_city.polygon_arr)]
    end
    return fig_zones
end

function get_plot_route_zones_zoom(zones_cluster, chosen_city, x_lims, y_lims)
    zones_cluster_city = zones_cluster[findall(x->x.city==chosen_city, zones_cluster)]
    zones_used = Vector{String}()
    sort!(zones_cluster_city, by=x->x.id)
    cols = distinguishable_colors(length(zones_cluster_city), [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    fig_zones2 = plot(fmt = :png)
    limit = 5000
    zones_limit = (limit<=length(zones_cluster_city) ? limit : length(zones_cluster_city))
    macromacrozone = "D"
    for i in 1:zones_limit
        #push!(zones_used, zones_cluster_city[i].id)
        zones_city = zones_cluster_city[i]
        println("ZONES_CITY: "*string(zones_city.id))
        if(zones_city.id[1:1] == macromacrozone)
            first_vertices = zones_city.polygon_arr[1].vertices[1]
            #first_vertices_lng = zones_city.polygon_arr[1].vertices[1][1]
            #first_vertices_lat = zones_city.polygon_arr[1].vertices[1][2]
            first_vertices_lng, first_vertices_lat = get_mean_point(zones_city.polygon_arr)
            #plot!(fig_zones2, zones_city.polygon_arr[1], xlims=x_lims, ylims=y_lims, color=cols[i], annotation=(first_vertices_lng,first_vertices_lat, Plots.text(zones_city.id,12)), labels=zones_city.id, legends=(0.1,0.9))
            #[plot!(fig_zones2, zones_city.polygon_arr[j], xlims=x_lims, ylims=y_lims, color=cols[i]) for j in 2:length(zones_city.polygon_arr)]
            plot!(fig_zones2, zones_city.polygon_arr[1], xlims=x_lims, ylims=y_lims, color=cols[i], annotation=(first_vertices_lng,first_vertices_lat, Plots.text(zones_city.id,8)))
            [plot!(fig_zones2, zones_city.polygon_arr[j], xlims=x_lims, ylims=y_lims, color=cols[i]) for j in 2:length(zones_city.polygon_arr)]
        end
    end
    return fig_zones2
end

function is_same_zone(stops_arr1, stops_arr2)
    lat_arr1 = sort!([stops_arr1[i].lat for i in 1:length(stops_arr1)])
    lat_arr2 = sort!([stops_arr2[i].lat for i in 1:length(stops_arr2)])
    lng_arr1 = sort!([stops_arr1[i].lng for i in 1:length(stops_arr1)])
    lng_arr2 = sort!([stops_arr2[i].lng for i in 1:length(stops_arr2)])
    min_dim = minimum([length(stops_arr1), length(stops_arr2)])
    return all([lat_arr1[i]==lat_arr2[i] for i in 1:min_dim]) && all([lng_arr1[i]==lng_arr2[i] for i in 1:min_dim])
end

function get_closest_zone(zones_centers, ref_point, params)
    if(params[13] == 1)
        _, closest_zone = findmin([get_manhattan_distance(zones_centers[i],ref_point) for i in 1:length(zones_centers)])
    else
        _, closest_zone = findmin([get_euclidian_distance(zones_centers[i],ref_point) for i in 1:length(zones_centers)])
    end
    return closest_zone
end

function zone_format(zone_id_arr)
    return string(zone_id_arr[1])*"-"*string(zone_id_arr[2])*"."*string(zone_id_arr[3])*string(zone_id_arr[4])
end

function update_valroutes!(route, old_zone_id, new_zone_id)
    curr_stops = route.stops[findall(x->zone(x)==old_zone_id,route.stops)]
    for S in curr_stops
        S.zone_id = Zone(new_zone_id)
    end
end

function get_score_prediction(zones_classif, zones_cluster, params)
    zscores = ZonesScores(params[2])
    fig_arr = Vector{Plots.Plot}()
    cols = distinguishable_colors(5, [RGB(1,1,1), RGB(0,0,0)], dropseed=true)
    error_count = 0
    for i in 1:params[2]
        zone_pred = zones_classif.zones_predicted[i]
        zone_act = zones_classif.zones_actual[i]
        if(zone_pred=="NaN" || zone_act=="NaN")
            if(zone_pred==zone_act)
                zscores.score_beta += 1
                zscores.score_gamma += 1
                zscores.score_zones += 1
            end
        else
            macrozone_pred = split(zone_pred,'.')[1]
            microzone_pred = split(zone_pred,'.')[2]
            macrozone_act = split(zone_act,'.')[1]
            microzone_act = split(zone_act,'.')[2]
            if(macrozone_pred[1] == macrozone_act[1])
                zscores.score_alpha = zscores.score_alpha+1
            end
            if(macrozone_pred == macrozone_act)
                zscores.score_beta = zscores.score_beta+1
            end
            if(microzone_pred == microzone_act)
                zscores.score_gamma = zscores.score_gamma+1
            end
            if(microzone_pred[1] == microzone_act[1])
                zscores.score_delta = zscores.score_delta+1
            end
            if(zone_pred == zone_act)
                zscores.score_zones = zscores.score_zones+1
            else
                error_count = error_count + 1
                fig_aux = plot(fmt = :png)
                stop_city = decode(gc, SA[latitude(zones_classif.stops_used[i]), longitude(zones_classif.stops_used[i])]).city
                zone_pred_pos = findfirst(x->x.id==zone_pred && x.city==stop_city, zones_cluster)
                if(!isnothing(zone_pred_pos))
                    zone_clust = zones_cluster[zone_pred_pos]
                    [plot!(fig_aux, zone_clust.polygon_arr[j], color=cols[1]) for j in 1:length(zone_clust.polygon_arr)]
                end
                #zone_act_polygons = zones_cluster[findfirst(x->x.id==zone_act, zones_cluster)].polygon_arr
                #[plot!(fig, zone_act_polygons[j], color=cols[2]) for j in 1:length(zone_act_polygons)]
                zone_act_pos = findfirst(x->x.id==zone_act && x.city==stop_city, zones_cluster)
                if(!isnothing(zone_act_pos))
                    zone_act_draw = zones_cluster[zone_act_pos]
                    [plot!(fig_aux, zone_act_draw.polygon_arr[j], color=cols[2]) for j in 1:length(zone_act_draw.polygon_arr)]
                else
                    println("ZONE NOT REGISTERED")
                end
                #plot!(fig, zones_cluster[i].polygon_arr[1], color=cols[1])
                #poly = [[[x, y] for (x, y) in zip(zones_classif.stops_used[i].lng, zones_classif.stops_used[i].lat)] for j in [1:2]];
                poly = VPolygon([[zones_classif.stops_used[i].lng, zones_classif.stops_used[i].lat]])
                #deleteat!(zones_cluster[i].polygon_arr, findfirst(x->x==poly, zones_cluster[i].polygon_arr))
                #poly = VPolygon.(poly);
                plot!(fig_aux, poly, color=RGB(0,0,0))
                push!(fig_arr, fig_aux)
            end
        end
    end
    return zscores, fig_arr
end

function get_random_routes_selection(routes, params)
    routes_quantity = length(routes)
    Random.seed!(params[12])
    rng = MersenneTwister(params[12])
    if(params[13] == 1)
        perm_seg = randcycle(rng, routes_quantity)
        train_seg = perm_seg[1:params[11]]
        val_seg =  perm_seg[params[11]+1:routes_quantity]
    elseif(params[13] == 2)
        perm_seg = randperm(rng, routes_quantity)
        train_seg = perm_seg[1:params[11]]
        val_seg =  perm_seg[params[11]+1:routes_quantity]
    elseif(params[13] == 3)
        train_seg = randsubseq(rng, 1:routes_quantity, 0.7)
        val_seg = findall(x->!(x in train_seg), 1:routes_quantity)
    else
        perm_seg = shuffle(rng, Vector(1:routes_quantity))
        train_seg = perm_seg[1:params[11]]
        val_seg =  perm_seg[params[11]+1:routes_quantity]
    end
    return routes[train_seg], routes[val_seg]
end

function solve(alg::ProRandomApproach, routes::Vector{Route})
    #route_quantity = dParams[1] # {0, 10, 20, 50} = {All, N Routes}
    #max_zones = dParams[2] # {100, 150, 200}
    #min_polygon_count = dParams[3] # {1, 2, 3}
    #min_dist_sides = dParams[4] # {0.0002, 0.0004, 0.0005}
    #min_dist_vertex = dParams[5] # {0.0007, 0.0010, 0.00014}
    #tolerance = dParams[6] # {0.00001, 0.00003, 0.00005}
    #dif_limit_up = dParams[7] # {0.004, 0.005, 0.006}
    #dif_limit_down = dParams[8] # # {0.001, 0.002, 0.003}
    #distance_measure = dParams[9] # {1, 2} = {Manhattan, Euclidian}
    #update_valroutes = dParams[10] # {0, 1} = {True, False}
    #train_routes_quantity = dParams[11] # {N < Total}
    #random_seed = dParams[12] # {1000, 2500, 4000}
    #random_distribution = dParams[13] # {1, 2, 3, 4} = {RandCycle, RandPerm, RandSubSeq, Shuffle}
    dParams = [20, 150, 1, 0.0004, 0.0007, 0.00005, 0.004, 0.001, 2, 1, 6000, 4000, 4]
    return solve(alg, routes; params=dParams)
end

function solve(alg::ProRandomApproach, routes::Vector{Route}; params)

    trainroutes, valroutes = get_random_routes_selection(routes, params)
    zclust = fill_zones_cluster(trainroutes)

    if(params[1]==0 && params[10]==0])
        _, zones_classif = get_predicted_actual(valroutes, zones_cluster, params)
    elseif(params[1]!=0 && params[10]==0)
        _, zones_classif = get_predicted_actual(valroutes[1:params[1]], zones_cluster, params)
    elseif(params[1]==0 && params[10]!=0)
        routes, zones_classif = get_predicted_actual(valroutes, zones_cluster, params)
    else
        routes, zones_classif = get_predicted_actual(valroutes[1:params[1]], zones_cluster, params)
    end
    zscores, fig_arr = get_score_prediction(zones_classif, zones_cluster, params)

    return trainroutes, valroutes, zclassif, zscores, fig_arr
end