# 1. Procedural - OK

"""
    AdaptativeApproach <: AbstractIdentificationAlgorithm

Algorithm that identifies the already segmented zones per route while using
customized identifiers distributions and variational analysis based on mean
differentiation of training routes' Zid and direct application of documented
adaptative identification rules dynamically configured per city distintion.

Inputs: 
* trainroutes: List of Training Transports Routes (unsequenced).
* valroutes: List of Validation/Testing Transports Routes (unsequenced).

or

* valroutes: List of Validation/Testing Transports Routes (unsequenced).
* zvar: Zone Identifiers Variation on Training Routes (struct ZonesVar)

Outputs:
* routes: List of updated routes (i.e. stops with identified zones).
* zones_arr: List of zones (i.e. Segmented Stops with respective Zone Identifiers)
* zones_ref_arr: List of Reference Zones of processed Validation/Testing Routes.

"""

using CSV, DataFrames, Dates, RouteSequencing
using LazySets, Plots, Colors, ReverseGeocode
using StatsBase, StaticArrays, NamedArrays, LinearAlgebra

struct AdaptativeApproach <: AbstractIdentificationAlgorithm end

struct DZones
    ids::Vector{String}
    center::Vector{Float64}
    stops::Vector{Stop}
end

struct ZonesDistrib
    major_prev::String
    zone_prev::String
    major_next::String
    zone_next::String
    city::String
    dist_lat::Vector{Float64}
    dist_lng::Vector{Float64}
end

ZonesDistrib(m_prev::String,z_prev::String,m_next::String,z_next::String,city::String,lat,lng) = ZonesDist(m_prev,z_prev,m_next,z_next,city,lat,lng)

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

function is_equal_zone_distrib(zd1, zd2)
    is_same_major_prev = zd1.major_prev == zd2.major_prev
    is_same_zone_prev = zd1.zone_prev == zd2.zone_prev
    is_same_major_next = zd1.major_next == zd2.major_next
    is_same_zone_next = zd1.zone_next == zd2.zone_next
    is_same_city = zd1.city == zd2.city
    return is_same_major_prev && is_same_zone_prev && is_same_major_next && is_same_zone_next && is_same_city
end

function is_equal_zone_clust(zc1, zc2)
    is_same_id = zc1.id == zc2.id
    is_same_city = zc1.city == zc2.city
    return is_same_id && is_same_city
end

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
    # sum(dist_lat)/length(dist_lat)
    return dist_lng_mean/length(dist_lng), dist_lat_mean/length(dist_lat)
end

function fill_zones_distrib(zones_cluster, closeness_measure)
    zones_distrib2 = Vector{ZonesDistrib}()
    zone_count = 1
    #closeness_measure = 0.01

    for ZC in zones_cluster
        zc_major = split(ZC.id,".")[1]
        zc_id = split(ZC.id,".")[2]
        zc_bari = get_baricenter(ZC.polygon_arr)
        zones_close = zones_cluster[findall(x->get_euclidian_distance(zc_bari, get_baricenter(x.polygon_arr))<=closeness_measure,zones_cluster)]
        if(zone_count%100==0)
            println(string(zone_count)*" - ZC: "*string(ZC.id))
            println("LENGTH ZONES_CLOSE: "*string(length(zones_close)))
        end
        for Z_close in zones_close
            z_close_split = split(Z_close.id,".")
            if(Z_close.id != ZC.id && Z_close.city == ZC.city && length(z_close_split)>1) # zc_major == z_close_major &&
                z_close_major = z_close_split[1]
                z_close_id = z_close_split[2]
                z_close_bari = get_baricenter(Z_close.polygon_arr)
                dist_lng = z_close_bari[1] - zc_bari[1]
                dist_lat = z_close_bari[2] - zc_bari[2]

                dist_lng_vec = Vector{Float64}()
                push!(dist_lng_vec, dist_lng)
                dist_lat_vec = Vector{Float64}()
                push!(dist_lat_vec, dist_lat)

                new_zd = ZonesDistrib(zc_major, zc_id, z_close_major, z_close_id, ZC.city, dist_lat_vec, dist_lng_vec)
                new_zd_pos = findfirst(x->is_equal_zone_distrib(new_zd, x), zones_distrib2)
                if(new_zd_pos==nothing)
                    push!(zones_distrib2, new_zd)
                else
                    push!(zones_distrib2[new_zd_pos].dist_lng, new_zd.dist_lng[1])
                    push!(zones_distrib2[new_zd_pos].dist_lat, new_zd.dist_lat[1])
                end
            end
        end
        #end
        zone_count = zone_count + 1
    end
    return zones_distrib2
end

function process_zones_distrib_variations(zones_distrib2)
    zvar = ZonesVar()
    #v_alpha_pos = [0.0, 0.0] #v_alpha_neg = [0.0, 0.0] #v_beta_pos = [0.0, 0.0] #v_beta_neg = [0.0, 0.0] #v_gamma_pos = [0.0, 0.0] #v_gamma_neg = [0.0, 0.0] #v_delta_pos = [0.0, 0.0] #v_delta_neg = [0.0, 0.0] #v_count = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]

    for ZD in zones_distrib2
        macromacro_dist = Int(only(ZD.major_next[1:1]))-Int(only(ZD.major_prev[1:1]))
        macro_dist = parse(Int64,ZD.major_next[2:end])-parse(Int64,ZD.major_prev[2:end])
        micro_dist = parse(Int64,ZD.zone_next[1:(end-1)])-parse(Int64,ZD.zone_prev[1:(end-1)])
        #println("NEXT: "*string(ZD.zone_next[(end):end]))
        #println("PREV: "*string(ZD.zone_prev[(end):end]))
        minimicro_dist = Int(only(ZD.zone_next[end:end]))-Int(only(ZD.zone_prev[end:end]))
        dist_lng_mean, dist_lat_mean = get_mean_dist(ZD.dist_lng, ZD.dist_lat)
        
        if(macromacro_dist>0)
            dist_lng_mean_mod = abs(macromacro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(macromacro_dist)*dist_lat_mean
            zvar.v_alpha_pos[1] = zvar.v_alpha_pos[1] + dist_lng_mean_mod #(v_alpha_pos[1] + dist_lng_mean)/2
            zvar.v_alpha_pos[2] = zvar.v_alpha_pos[2] + dist_lat_mean_mod #(v_alpha_pos[2] + dist_lat_mean)/2
            zvar.v_count[1] = zvar.v_count[1]+1
        end
        if(macromacro_dist<0)
            dist_lng_mean_mod = abs(macromacro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(macromacro_dist)*dist_lat_mean
            zvar.v_alpha_neg[1] = zvar.v_alpha_neg[1] + dist_lng_mean_mod #(v_alpha_neg[1] + dist_lng_mean)/2
            zvar.v_alpha_neg[2] = zvar.v_alpha_neg[2] + dist_lat_mean_mod #(v_alpha_neg[2] + dist_lat_mean)/2
            zvar.v_count[2] = zvar.v_count[2]+1
        end
        if(macro_dist>0)
            dist_lng_mean_mod = abs(macro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(macro_dist)*dist_lat_mean
            zvar.v_beta_pos[1] = zvar.v_beta_pos[1] + dist_lng_mean_mod #(v_beta_pos[1] + dist_lng_mean)/2
            zvar.v_beta_pos[2] = zvar.v_beta_pos[2] + dist_lat_mean_mod #(v_beta_pos[2] + dist_lat_mean)/2
            zvar.v_count[3] = zvar.v_count[3]+1
        end
        if(macro_dist<0)
            dist_lng_mean_mod = abs(macro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(macro_dist)*dist_lat_mean
            zvar.v_beta_neg[1] = zvar.v_beta_neg[1] + dist_lng_mean_mod #(v_beta_neg[1] + dist_lng_mean)/2
            zvar.v_beta_neg[2] = zvar.v_beta_neg[2] + dist_lat_mean_mod #(v_beta_neg[2] + dist_lat_mean)/2
            zvar.v_count[4] = zvar.v_count[4]+1
        end
        if(micro_dist>0)
            dist_lng_mean_mod = abs(micro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(micro_dist)*dist_lat_mean
            zvar.v_gamma_pos[1] = zvar.v_gamma_pos[1] + dist_lng_mean_mod #(v_gamma_pos[1] + dist_lng_mean)/2
            zvar.v_gamma_pos[2] = zvar.v_gamma_pos[2] + dist_lat_mean_mod #(v_gamma_pos[2] + dist_lat_mean)/2
            zvar.v_count[5] = zvar.v_count[5]+1
        end
        if(micro_dist<0)
            dist_lng_mean_mod = abs(micro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(micro_dist)*dist_lat_mean
            zvar.v_gamma_neg[1] = zvar.v_gamma_neg[1] + dist_lng_mean_mod #(v_gamma_neg[1] + dist_lng_mean)/2
            zvar.v_gamma_neg[2] = zvar.v_gamma_neg[2] + dist_lat_mean_mod #(v_gamma_neg[2] + dist_lat_mean)/2
            zvar.v_count[6] = zvar.v_count[6]+1
        end
        if(minimicro_dist>0)
            dist_lng_mean_mod = abs(minimicro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(minimicro_dist)*dist_lat_mean
            zvar.v_delta_pos[1] = zvar.v_delta_pos[1] + dist_lng_mean_mod #(v_delta_pos[1] + dist_lng_mean)/2
            zvar.v_delta_pos[2] = zvar.v_delta_pos[2] + dist_lat_mean_mod #(v_delta_pos[2] + dist_lat_mean)/2
            zvar.v_count[7] = zvar.v_count[7]+1
        end
        if(minimicro_dist<0)
            dist_lng_mean_mod = abs(minimicro_dist)*dist_lng_mean
            dist_lat_mean_mod = abs(minimicro_dist)*dist_lat_mean
            zvar.v_delta_neg[1] = zvar.v_delta_neg[1] + dist_lng_mean_mod #(v_delta_neg[1] + dist_lng_mean)/2
            zvar.v_delta_neg[2] = zvar.v_delta_neg[2] + dist_lat_mean_mod #(v_delta_neg[2] + dist_lat_mean)/2
            zvar.v_count[8] = zvar.v_count[8]+1
        end
        
    end

    return normalize_vector(zvar)
    
end

function fill_zones_distrib_dict()
    zones_distrib_dict = Dict() 
    #Dict{String,Tuple{String,String,String,String,String,Vector{Float64},Vector{Float64}}}()
    for i in 1:length(zones_distrib2)
        ZD = zones_distrib2[i]
        dict_to_add = Dict("m_prev"=>ZD.major_prev, "z_prev"=>ZD.zone_prev, 
            "m_next"=>ZD.major_next, "z_next"=>ZD.zone_next, "city"=>ZD.city, "d_lng"=>ZD.dist_lng, "d_lat"=>ZD.dist_lat)
        zones_distrib_dict["Z_D_"*string(i)] = dict_to_add
    end
    return zones_distrib_dict
end

function detect_repeated_values(zones_distrib)
    repeated_zd = Vector{ZonesDistrib}()
    for ZD in zones_distrib
        all_elements = findall(x->is_equal_zone_distrib(x,ZD),zones_distrib)
        if(length(all_elements)>1)
            push!(repeated_zd, ZD)
        end
    end
    return repeated_zd
end

function get_proj(a, b)
    axb = dot(a,b)
    bxb = dot(b,b)
    bproj = (axb/bxb)*b
    bproj_mod = sqrt(bproj[1]^2 + bproj[2]^2)
    return bproj_mod
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
        #println("ZONES_CITY: "*string(zones_city.id))
        #println("FIRST LETTER: "*string(zones_cluster_city[i].id[1:1]))
        #println("FIRST_LETTER?: "*string(zones_city.id[1:1] == macromacrozone))
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
        #println("FIRST LETTER: "*string(zones_cluster_city[i].id[1:1]))
        #println("FIRST_LETTER?: "*string(zones_city.id[1:1] == macromacrozone))
        if(zones_city.id[1:1] == macromacrozone)
            println("ENTRO IF")
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

function is_same_zone(stops_arr1, stops_arr2)
    lat_arr1 = sort!([stops_arr1[i].lat for i in 1:length(stops_arr1)])
    lat_arr2 = sort!([stops_arr2[i].lat for i in 1:length(stops_arr2)])
    lng_arr1 = sort!([stops_arr1[i].lng for i in 1:length(stops_arr1)])
    lng_arr2 = sort!([stops_arr2[i].lng for i in 1:length(stops_arr2)])
    min_dim = minimum([length(stops_arr1), length(stops_arr2)])
    return all([lat_arr1[i]==lat_arr2[i] for i in 1:min_dim]) && all([lng_arr1[i]==lng_arr2[i] for i in 1:min_dim])
end

function get_mid_points(curr_stops)
    zone_lng_max = maximum([curr_stops.lng for i in 1:length(curr_stops)])
    zone_lng_min = minimum([curr_stops.lng for i in 1:length(curr_stops)]) 
    zone_lat_max = maximum([curr_stops.lat for i in 1:length(curr_stops)])
    zone_lat_min = minimum([curr_stops.lat for i in 1:length(curr_stops)])
    return [(zone_lng_max+zone_lng_min)/2, (zone_lat_max+zone_lat_min)/2]
end

function get_zones_centers(route, route_length, route_width, params)
    zones_centers = Vector{Vector{Float64}}()
    for Z in zones(route)
        curr_stops = route.stops[findall(x->zone(x)==Z,route.stops)]
        if(params[12] == 1)
            curr_center = get_mean_dist([curr_stops[i].lng for i in 1:length(curr_stops)], [curr_stops[i].lat for i in 1:length(curr_stops)])
        elseif(params[12] == 2)
            curr_center = get_mid_points(curr_stops) 
        elseif(params[12] == 3)
            curr_center = (route_width/route_length)*get_mean_dist([curr_stops[i].lng for i in 1:length(curr_stops)], [curr_stops[i].lat for i in 1:length(curr_stops)])
        else
            curr_center = (route_width/route_length)*get_mid_points(curr_stops) 
        end
        push!(zones_centers, curr_center)
    end
    return zones_centers
end

function get_closest_zone(zones_centers, ref_point, params)
    if(params[13] == 1)
        _, closest_zone = findmin([get_manhattan_distance(zones_centers[i],ref_point) for i in 1:length(zones_centers)])
    else
        _, closest_zone = findmin([get_euclidian_distance(zones_centers[i],ref_point) for i in 1:length(zones_centers)])
    end
    return closest_zone
end

function get_ref_zone(route, params) # Get Closest Point To
    route_lng_max = maximum([route.stops[i].lng for i in 1:length(route.stops)])
    route_lng_min = minimum([route.stops[i].lng for i in 1:length(route.stops)]) 
    route_length = route_lng_max - route_lng_min
    route_lat_max = maximum([route.stops[i].lat for i in 1:length(route.stops)])
    route_lat_min = minimum([route.stops[i].lat for i in 1:length(route.stops)]) 
    route_width = route_lat_max - route_lat_min
    zones_centers = get_zones_centers(route, route_length, route_width, params)
    if(params[9]==1)
        ref_point = [route_lng_min, route_lat_max]
    elseif(params[9]==2)
        ref_point = [route_lng_min + route_length/2, route_lat_max]
    elseif(params[9]==3)
        ref_point = [route_lng_min + route_length, route_lat_max]
    elseif(params[9]==4)
        ref_point = [route_lng_min, route_lat_max - route_width/2]
    elseif(params[9]==5)
        ref_point = [route_lng_min + route_length/2, route_lat_max - route_width/2]
    elseif(params[9]==6)
        ref_point = [route_lng_min + route_length, route_lat_max - route_width/2]
    elseif(params[9]==7)
        ref_point = [route_lng_min, route_lat_max - route_width]
    elseif(params[9]==8)
        ref_point = [route_lng_min + route_length/2, route_lat_max - route_width]
    else
        ref_point = [route_lng_min + route_length, route_lat_max - route_width]
    end
    
    closest_zone_num = get_closest_zone(zones_centers, ref_point, params)
    ref_zone_center = zones_centers[closest_zone_num]
    ref_zone_stops = route.stops[findall(x->zone(x)==zones(route)[closest_zone_num],route.stops)]
    ref_zone = DZones(Vector{String}(),ref_zone_center,ref_zone_stops)
    return ref_zone, zones_centers
end

function get_ref_identifier(route, ref_center, params)
    ref_id = ['B', '1', '1', 'D']
    route_lng_max = maximum([route.stops[i].lng for i in 1:length(route.stops)])
    route_lng_min = minimum([route.stops[i].lng for i in 1:length(route.stops)]) 
    route_length = route_lng_max - route_lng_min
    route_lat_max = maximum([route.stops[i].lat for i in 1:length(route.stops)])
    route_lat_min = minimum([route.stops[i].lat for i in 1:length(route.stops)]) 
    route_width = route_lat_max - route_lat_min
    ref_id[1] = ref_id[1] + fld((route_lng_min + route_length/2) - ref_center[1], 1.0) #ALPHA
    ref_id[2] = ref_id[2] + fld((route_lng_min + route_length/3) - ref_center[2], 1.0) #BETA
    ref_id[3] = ref_id[3] + fld(abs(ref_center[1]-params[10]*ref_center[2]), params[11]) # GAMMA
    ref_id[4] = ref_id[4] + fld(abs(-ref_center[1]-params[10]*ref_center[2]), params[11]) # DELTA
    return ref_id
end

function manage_identifiers(zone_proj_mods, ref_zone, params)
    zone_layer_ids = ref_zone.ids
    comb_num = 1
    for i in 1:2:7;
        zone_comb = params[comb_num]*(zone_proj_mods[i] - zone_proj_mods[i+1])
        if(abs(zone_comb) >= params[comb_num+4])
            zone_layer_ids[comb_num] = zone_layer_ids[comb_num] + fld(zone_comb-params[comb_num+4], 1.0)
        end
        comb_num = comb_num + 1
    end
    return zone_layer_ids
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

function value_in_between(value, min_val, max_val)
    return value >= min_val && value <= max_val
end

function stops_in_region(stops, lat_range, lng_range)
    return !isnothing(findfirst(x-> value_in_between(x.lat, lat_range[1], lat_range[2]) && value_in_between(x.lng, lng_range[1], lng_range[2]), stops))
end

function get_adaptative_parameters(route)
    params = Vector{Float64}()
    route_lng_max = maximum([route.stops[i].lng for i in 1:length(route.stops)])
    route_lng_min = minimum([route.stops[i].lng for i in 1:length(route.stops)]) 
    route_length = abs(route_lng_max - route_lng_min)
    route_lat_max = maximum([route.stops[i].lat for i in 1:length(route.stops)])
    route_lat_min = minimum([route.stops[i].lat for i in 1:length(route.stops)]) 
    route_width = abs(route_lat_max - route_lat_min)
    rel_route = route_length/route_width
    push!(params, round(1.2*rel_route,digits=3)) # comb_alpha
    push!(params, round(0.9*rel_route,digits=3)) #comb_beta
    push!(params, round(0.6*rel_route,digits=3)) #comb_gamma
    push!(params, round(0.5*rel_route,digits=3)) #comb_delta
    if(stops_in_region(route.stops,[30.134, 30.517],[-97.951, -97.531]))      # Austin
        push!(params, fld(route_length, 1.4)) #compareVal_alpha
        push!(params, fld(route_width, 2.1)) #compareVal_beta
        push!(params, fld(route_length, 2.5)) #compareVal_gamma
        push!(params, fld(route_width, 3.2)) #compareVal_delta
        push!(params, 4) # initial_zone_location
        push!(params, -3) # initial_identifier_slope
        push!(params, 3) # initial_identifier_scale
        push!(params, 2) # zone_center_calculation
        push!(params, 2) # distance_measure
    elseif(stops_in_region(route.stops,[42.235, 42.397],[-71.199, -70.967]))  # Boston
        push!(params, fld(route_length, 1.6)) #compareVal_alpha
        push!(params, fld(route_width, 1.9)) #compareVal_beta
        push!(params, fld(route_length, 2.2)) #compareVal_gamma
        push!(params, fld(route_width, 2.7)) #compareVal_delta
        push!(params, 2) # initial_zone_location
        push!(params, -1) # initial_identifier_slope
        push!(params, 3) # initial_identifier_scale
        push!(params, 1) # zone_center_calculation
        push!(params, 2) # distance_measure
    elseif(stops_in_region(route.stops,[41.653, 42.015],[-87.866, -87.487]))  # Chicago
        push!(params, fld(route_length, 1.4)) #compareVal_alpha
        push!(params, fld(route_width, 2.3)) #compareVal_beta
        push!(params, fld(route_length, 2.6)) #compareVal_gamma
        push!(params, fld(route_width, 3.1)) #compareVal_delta
        push!(params, 6) # initial_zone_location
        push!(params, -4) # initial_identifier_slope
        push!(params, 2) # initial_identifier_scale
        push!(params, 4) # zone_center_calculation
        push!(params, 1) # distance_measure
    elseif(stops_in_region(route.stops,[33.748, 34.312],[-118.639, -118.152]))  # Los Angeles
        push!(params, fld(route_length, 1.3)) #compareVal_alpha
        push!(params, fld(route_width, 2.4)) #compareVal_beta
        push!(params, fld(route_length, 2.7)) #compareVal_gamma
        push!(params, fld(route_width, 2.9)) #compareVal_delta
        push!(params, 5) # initial_zone_location
        push!(params, -4) # initial_identifier_slope
        push!(params, 4) # initial_identifier_scale
        push!(params, 1) # zone_center_calculation
        push!(params, 2) # distance_measure
    elseif(stops_in_region(route.stops,[47.502, 47.732],[-122.435, -122.244]))  # Seattle
        push!(params, fld(route_length, 1.4)) #compareVal_alpha
        push!(params, fld(route_width, 2.3)) #compareVal_beta
        push!(params, fld(route_length, 2.8)) #compareVal_gamma
        push!(params, fld(route_width, 3.3)) #compareVal_delta
        push!(params, 7) # initial_zone_location
        push!(params, -3) # initial_identifier_slope
        push!(params, 3) # initial_identifier_scale
        push!(params, 3) # zone_center_calculation
        push!(params, 2) # distance_measure
    else                                                                        # Montevideo
        push!(params, fld(route_length, 1.5)) #compareVal_alpha
        push!(params, fld(route_width, 2.3)) #compareVal_beta
        push!(params, fld(route_length, 2.7)) #compareVal_gamma
        push!(params, fld(route_width, 3.1)) #compareVal_delta
        push!(params, 8) # initial_zone_location
        push!(params, -1) # initial_identifier_slope
        push!(params, 4) # initial_identifier_scale
        push!(params, 3) # zone_center_calculation
        push!(params, 2) # distance_measure
    end
    return params
end

function solve(alg::AdaptativeApproach, trainroutes::Vector{Route}, valroutes::Vector{Route})
    #comb_alpha= dParams[1] # {1.8, 2.2, 2.4}
    #comb_beta = dParams[2] # {1.5, 1.7, 1.8}
    #comb_gamma = dParams[3] # {1.3, 1.4, 1.7}
    #comb_delta = dParams[4] # {1.1, 1.3, 1.5}
    #compareVal_alpha = dParams[5] # {0.467, 0.614, 1.105}
    #compareVal_beta = dParams[6] # {0.231, 0.378, 0.681}
    #compareVal_gamma = dParams[7] # {0.112, 0.175, 0.245}
    #compareVal_delta = dParams[8] # {0.065, 0.089, 0.134}
    zclust = fill_zones_cluster(trainroutes)
    zdistr = fill_zones_distrib(zclust, 0.01)
    zvar = process_zones_distrib_variations(zdistr) # Can be filled with adaptability for cities
    return solve(alg, valroutes, zvar)
end

function solve(alg::AdaptativeApproach, routes::Vector{Route}, zvar::ZonesVar)

    vector_variations = 1:4
    zones_arr = Vector{DZones}()
    zones_ref_arr = Vector{DZones}()
    
    for R in routes
        params = get_adaptative_parameters(R)
        route_num = findfirst(x->name(x)=name(R), routes)
        zone_quantity = length(unique(zones(R)))
        ref_zone, zones_centers = get_ref_zone(R, params)
        ref_zone.ids = get_ref_identifier(R, ref_zone.center, params)
        push!(zones_arr, ref_zone)
        push!(zones_ref_arr, ref_zone)

        for Z in zones(R)
            curr_zone = DZones(Vector{String}(),zeros(2,1),Vector{Stop}())
            curr_zone.center = zones_centers[findfirst(x->x==Z, zones(R))]
            zone_proj_mods = get_proj_mods(zvar, curr_zone.center)
            curr_zone.ids = manage_identifiers(zone_proj_mods, ref_zone, params)
            #curr_zone.ids = [ manage_identifiers(zone_proj_mods, v, zvar, ref_zone) for v in vector_variations ]
            curr_zone.stops = R.stops[findall(x->x.zone==Z,R.stops)]
            push!(zones_arr, curr_zone)
            update_valroutes!(R, Z, zone_format(curr_zone.ids))
        end
    end
    return routes, zones_arr, zones_ref_arr
end