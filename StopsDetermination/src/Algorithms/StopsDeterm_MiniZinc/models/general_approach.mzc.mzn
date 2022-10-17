% SOLUTION FORMAT CHECKER: DO NOT EDIT!!!

% --- DATA INPUT ---
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

% --- DATA OUTPUT ---
set of int: BINARY = {0, 1};
array[REGISTRIES] of BINARY: stopsdeterm;
int: stopsquantity;
int: _objective;

% --- SUCCESS MESSAGE ---
output["SUCCESS: Your model contains all variables required..."];