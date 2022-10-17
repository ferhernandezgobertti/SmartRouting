## Matriz de tiempos por zona

```julia
using RouteSequencing

# cargar la primer ruta de las 13 de prueba
R = parse_route_data()[1];

@show RouteSequencing.name(R)
# cargar las secuencias actuales/exactas
S = parse_sequence()[1];

# reordenar R de acuerdo a S
apply_sequence!(R, S);

# mostrar las 5 primeras zonas recorridas (sin incluir la station)
@show zones(R, remove_station=true)[1:5]

# indices en el vector de rutas de acuerdo a la zona
zoneidx, _, _ = split_by_zone(R);

@show zoneidx[1]
 
# tomar los nombres de las paradas 2 a 7 (primera zona)
@show stopnames = [RouteSequencing.name(s) for s in R.stops[2:7]]

# cargar la matriz de tiempos de la primer ruta
T = parse_travel_times()[1].travel_times;

# seleccionar solamente la submatriz que consiste en las paradas 2 a 7 en la ruta ordenada
M = view(T, stopnames, stopnames)
```

Salida:

```julia
RouteSequencing.name(R) = "bcc07fea-86d2-41e4-9a58-cfc78956dcc7"

(zones(R, remove_station = true))[1:5] = ["H24.2C", "H24.1C", "H24.1B", "H24.2B", "H24.3B"]

zoneidx[1] = [1, 2, 3, 4, 5, 6]

stopnames = [RouteSequencing.name(s) for s = R.stops[2:7]] = ["LZ", "ZT", "ZI", "YC", "WZ", "HS"]

6×6 Named SubArray{Float64, 2, Matrix{Float64}, Tuple{Vector{Int64}, Vector{Int64}}, false}
A ╲ B │    LZ     ZT     ZI     YC     WZ     HS
──────┼─────────────────────────────────────────
LZ    │   0.0   29.1  220.8   49.7   74.7  119.3
ZT    │  38.9    0.0   66.3   20.5   45.5   90.1
ZI    │ 223.6   75.5    0.0   68.0   93.0  137.6
YC    │  55.4   19.1   77.4    0.0   24.9   69.5
WZ    │  80.5   44.2  102.5   25.0    0.0   44.5
HS    │ 137.4  101.1  159.4   81.9   56.8    0.0
```

## Optimizacion (minimum length Hamiltonian path)

```julia
using GLPK
using RouteSequencing: _min_length_path_mip

_min_length_path_mip(T, R.stops[2:7], optimizer=GLPK.Optimizer)
```

## Docker en Fedora

```bash
 sudo systemctl start docker
```
