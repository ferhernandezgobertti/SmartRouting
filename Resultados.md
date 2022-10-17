## Tabla de resultados

|Algorithm|Minimum|Maximum|Mean|
|-----|-----|----|-----|
|MinTime - LocalTime |0.0388|0.224|0.12|
|MinTimeZ - LocalTime_Zones (*)|2.41|9.09|5.03|
|MinTimeZ - LocalTime_Zones (**)|0.0013|0.0036|0.011|
|MinTimeFinal - GlobalTime_IF_LI_MIP (*)|0.03548|2.3857|0.6045|
|MinTimeFinal - GlobalTime_IF_LI_MIP (**)|0.0031|0.8389|0.2173|
|MinTimeFinal2 - GlobalTime_IF_FO_MIP (*)|0.0353|2.3835|0.6166|
|MinTimeFinal2 - GlobalTime_IF_FO_MIP (**)|0.0033|0.8160|0.2194|
|MinTimeFinal1B - GlobalTime_I_LI_MIP (*)|0.0495|2.8726|0.7232|
|MinTimeFinal1B - GlobalTime_I_LI_MIP (**)|0.0053|1.2328|0.3018|
|MinTimeFinal3 - GlobalTime_IF_LI (*)|0.0336|2.3814|0.6032|
|MinTimeFinal3 - GlobalTime_IF_LI (**)|0.0028|0.8228|0.2077|

- MinTime -- Este algoritmo no decide en base a zonas sino que usa la tabla de MinTime local.
 
- MinTimeZ (*) -- Zones are not sorted, just pick the next one in the given route.

- MinTimeZ (**) -- Zones are sorted correctly according to the actual sequence.

- MinTimeFinal (*) -- Zones sorted by Hausdorf distance, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP.

- MinTimeFinal (**) -- Zones sorted correctly according to actual sequence, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP.

- MinTimeFinal2 (*) -- Zones sorted by Hausdorf distance, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == First_Outer_Stop and Next_Stop decided base on Min(Local_Travel_Time), using MIP.

- MinTimeFinal2 (**) -- Zones sorted correctly according to actual sequence, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == First_Outer_Stop and Next_Stop decided base on Min(Local_Travel_Time), using MIP.

- MinTimeFinal1B (*) -- Zones sorted by Hausdorf distance, Stops ordered per Zone as [Min(Global_Travel_Time), First_Stop_Fixed, Last_Stop_Iterated] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP.

- MinTimeFinal1B (**) -- Zones sorted correctly according to actual sequence, Stops ordered per Zone as [Min(Global_Travel_Time), First_Stop_Fixed, Last_Stop_Iterated] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP.

- MinTimeFinal3 (*) -- Zones sorted by Hausdorf distance, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP and Exhaustive depending on quantity of Stops per Zone.

- MinTimeFinal3 (**) -- Zones sorted correctly according to actual sequence, Stops ordered per Zone as [Min(Global_Travel_Time), First&Last_Stop_Fixed] considering 1st_Subsequence_Stop == Last_Inner_Stop, using MIP and Exhaustive depending on quantity of Stops per Zone.


## Evaluation script


```julia 
# ------------------------------------------
# Cargar datos
# ------------------------------------------
using RouteSequencing

# cargar datos de travel times
data_times = parse_travel_times();

# cargar datos de rutas
data_route = parse_route_data();

# ------------------------------------------
# Calcular secuencias
# ------------------------------------------

# -- algoritmo de tiempo minimo ---
alg = MinTime()

# -- algoritmo de tiempo minimo por zonas ---
# alg = MinTimeZ()
# data_sequence = parse_sequence()
# apply_sequence!.(data_route, data_sequence)

# calcular la secuencia
sequences = [solve(alg, data_route[i], data_times[i]) for i in 1:length(data_route)];

# ------------------------------------------
# Evaluar score y calcular valores extremos
# ------------------------------------------
using Statistics

result = get_score(sequence)

@show mean(result)
@show extrema(result)
```
