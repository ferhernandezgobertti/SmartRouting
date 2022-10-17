% SOLUTION FORMAT CHECKER: DO NOT EDIT!!!

% --- DATA INPUT ---
int: nstops;
int: nparams;
set of int: STOPS = 1..nstops;
set of int: PARAMS = 1..nparams;
%array[PARAMS] of float: dparams; % Deployed on following integers

int: minnum;
int: maxnum;
int: minsize;
int: maxsize;
int: minstops;
int: maxstops;
int: minsplit;
int: maxsplit;
int: radiuseps;
int: radiusquant;
int: tolerance;

% Coordinates
set of float: LATITUDE = -90.0..90.0;
array[STOPS] of LATITUDE: lat;
set of float: LONGITUDE = -180.0..180.0;
array[STOPS] of LONGITUDE: lon; 

% Cost Matrices
array[STOPS,STOPS] of float: costmatrix;


% --- DATA OUTPUT ---
array[STOPS] of int: clustass;
int: zonesquantity;
int: _objective;

% --- SUCCESS MESSAGE ---
output["SUCCESS: Your model contains all variables required..."];