int: nreg;
int: nparams;
int: minstops;
int: maxstops;
set of int: REGISTRIES = 1..nreg;
set of int: PARAMS = 1..nparams;
array[PARAMS] of float: dparams; 

% Distance
set of float: LATITUDE = -90.0..90.0;
array[REGISTRIES] of LATITUDE: lat;
set of float: LONGITUDE = -180.0..180.0;
array[REGISTRIES] of LONGITUDE: lon; 

% Time Duration
set of float: TIME = 1.0..86400.0; % Min 1 second, Max 24*60*60 seconds
array[REGISTRIES] of TIME: time; 

% Velocity
set of float: VELOCITY = 0.0..60.0;
array[REGISTRIES] of VELOCITY: vel; 

% Kilometer
set of float: KILOMETER = 0.0..1000.0;
array[REGISTRIES] of KILOMETER: km; 

% Events (Encendido, Apagado, Otro/Other)
enum EVENT = { E, A, O };
array[REGISTRIES] of EVENT: events;

set of int: BINARY = {0, 1};
array[REGISTRIES] of BINARY: stopsdeterm;
int: stopsquantity;

% General Constraints
constraint stopsquantity = sum(i in REGISTRIES)(stopsdeterm[i]);
constraint sum(i in REGISTRIES)(stopsdeterm[i]) >= minstops;
constraint sum(i in REGISTRIES)(stopsdeterm[i]) <= maxstops;

% Particular Constraints
set of int: RQM = 1..round(dparams[1]);
set of int: FORWARD = 1..round(nreg-dparams[1]);
% Kilometer Constraint
constraint forall(i in FORWARD)(sum(j in RQM)(km[i+j-1]*stopsdeterm[i+j-1]-km[i+j]*stopsdeterm[i+j]) < dparams[2]); 
% Velocity Constraint
constraint forall(i in FORWARD)(sum(j in RQM)(vel[i+j-1]*stopsdeterm[i+j-1]-vel[i+j]*stopsdeterm[i+j]) < dparams[3]);
% Temporal Constraint
constraint forall(i in FORWARD)(sum(j in RQM)(time[i+j-1]*stopsdeterm[i+j-1]-time[i+j]*stopsdeterm[i+j]) < dparams[4]);
% Distance Constraint
constraint forall(i in FORWARD)(sum(j in RQM)(sqrt((lat[i+j-1]*stopsdeterm[i+j-1]-lat[i+j]*stopsdeterm[i+j])^2+(lon[i+j-1]*stopsdeterm[i+j-1]-lon[i+j]*stopsdeterm[i+j])^2)) < dparams[5]);
% Events Constraint
constraint forall(i in FORWARD)(count(events[1..i], E)*stopsdeterm[i]-count(events[1..i], A)*stopsdeterm[i] == dparams[6]);

%ALTERNATIVES

%constraint forall(i in FORWARD)(exists(j in RQM)(km[i+j-1]*stopsdeterm[i+j-1]-km[i+j]*stopsdeterm[i+j]) < dparams[2]); 
%constraint forall(i in FORWARD)(exists(j in RQM)(vel[i+j-1]*stopsdeterm[i+j-1]-vel[i+j]*stopsdeterm[i+j]) < dparams[3]);
%constraint forall(i in FORWARD)(exists(j in RQM)(time[i+j-1]*stopsdeterm[i+j-1]-time[i+j]*stopsdeterm[i+j]) < dparams[4]);
%constraint forall(i in FORWARD)(exists(j in RQM)(sqrt((lat[i+j-1]*stopsdeterm[i+j-1]-lat[i+j]*stopsdeterm[i+j])^2+(lon[i+j-1]*stopsdeterm[i+j-1]-lon[i+j]*stopsdeterm[i+j])^2)) < dparams[5]);


%include "global_cardinality_fn";
%array[REGISTRIES] of var int: occurrences = global_cardinality(events, [i | i in REGISTRIES]);
%var int: num_diff = count(o in occurrences) (o > 0);


solve minimize sum(i in REGISTRIES)(stopsdeterm[i]);

output
      ["STOPS: "] ++ [ " " ++ show_int(2,stopsdeterm[i]) | i in REGISTRIES ]
      ++ ["\n"] ++
      [if stopsquantity > 0 then " Quantity: " ++ show_int(2,stopsquantity) else " ." endif]
      ++ ["\n"] ++
      [ if stopsdeterm[i] = 1 then " POS: " ++ show_int(1,i) ++ ": " ++ show_int(2,stopsdeterm[i]) else "" endif
      | i in REGISTRIES ]
      ++ ["\n"];
