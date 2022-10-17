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

% Outputs and Auxiliaries
array[STOPS] of int: clustass;
int: zonesquantity;
array[STOPS] of var bool: zonesoccur;

% General Constraints
constraint clustass[1] = 1; % Symmetry Reduction
constraint forall(i,j in STOPS) (clustass[i] = j -> zonesoccur[j]);
constraint zonesquantity = count(zonesoccur);
%ALTERNATIVE
%include "global_cardinality_fn";
%array[STOPS] of var int: zonesoccur = global_cardinality(clustass, [i | i in STOPS]);
%var int: zonesquantity = count(o in occurrences) (o > 0);

% Cluster Number Constraints
constraint forall(i in 1..zonesquantity)(sum(s in STOPS)(zonesoccur[s]) >= minnum);
constraint forall(i in 1..zonesquantity)(sum(s in STOPS)(zonesoccur[s]) <= maxnum);
%ALTERNATIVE
%constraint zonesquantity >= minnum;
%constraint zonesquantity <= maxnum;

% Cluster Sizes Constraints
array[1..maxnum] of var int: latmax;
constraint forall(i in 1..zonesquantity)(latmax[i] = max([lat[j] | j in STOPS where clustass[j] == i])); 
array[1..maxnum] of var int: latmin;
constraint forall(i in 1..zonesquantity)(latmin[i] = min([lat[j] | j in STOPS where clustass[j] == i])); 
array[1..maxnum] of var int: lonmax;
constraint forall(i in 1..zonesquantity)(lonmax[i] = max([lon[j] | j in STOPS where clustass[j] == i])); 
array[1..maxnum] of var int: lonmin;
constraint forall(i in 1..zonesquantity)(lonmin[i] = min([lon[j] | j in STOPS where clustass[j] == i])); 
constraint forall(i in 1..zonesquantity)(((latmax[i]-latmin[i])+(lonmax[i]-lonmin[i]))/2 >= minsize);
constraint forall(i in 1..zonesquantity)(((latmax[i]-latmin[i])+(lonmax[i]-lonmin[i]))/2 <= maxsize);

% Cluster Stops Constraints
include "at_least.mzn";
include "at_most.mzn";
constraint forall(i in 1..zonesquantity)(at_least(minstops, clustass, i));
constraint forall(i in 1..zonesquantity)(at_most(maxstops, clustass, i));
%ALTERNATIVE
%include "count.mzn";
%constraint forall(i in 1..zonesquantity)(count(clustass, i) >= minnum);
%constraint forall(i in 1..zonesquantity)(count(clustass, i) <= maxnum);

% Delta, Must/Cannot Link Constraints
array[1..maxnum,STOPS] of var int: distances;
constraint forall(i, j in STOPS where clustass[i]!=clustass[j])(distances[clustass[i],i] = costmatrix[i,j]);
constraint forall(i in 1..zonesquantity)(min([distances[i,j] | j in STOPS]) >= minsplit);
%constraint forall(i in 1..zonesquantity)(min([distances[i,j] | j in STOPS]) <= maxsplit);

% Epsilon Constraints
constraint forall(i in STOPS)(count([costmatrix[i,j] <= radiuseps | j in STOPS where clustass[j]==clustass[i]], 1) >= radiusquant);

% OPTIMIZATION CRITERIA: Min Distance
solve minimize sum(i in 1..zonesquantity)(min([distances[i,j] | j in STOPS]));

output
      ["ZONES: "] ++ [ " " ++ show_int(2,clustass[i]) | i in STOPS ]
      ++ ["\n"] ++
      [if zonesquantity > 0 then " Quantity: " ++ show_int(2,zonesquantity) else " ." endif]
      ++ ["\n"] ++
      [ if clustass[i] > 0 then " STOP : " ++ show_int(1,i) ++ ": (" ++ show_int(5,round(lat[i])) ++ ", " ++ show_int(5,round(lon[i])) 
      ++ ") - Zone Num: " ++ show_int(2,clustass[i]) else "" endif | i in STOPS ]
      ++ ["\n"];
