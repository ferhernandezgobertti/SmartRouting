```@meta
DocTestSetup = :(using RouteSequencing)
CurrentModule = RouteSequencing
```

# Hamiltonian paths

In this section we define the minimum-length Hamiltonian path problem and consider different solution strategies, their implementation and validation in a set of simple examples.


## Problem statement

Hamiltonian paths are graph paths between two vertices of a directed or undirected graph that visits each vertex exactly once. A graph that contains a Hamiltonian path (also known as traceable paths) is called a *traceable graph*. A graph is *Hamiltonian-connected* if for every pair of vertices there is a Hamiltonian path between the two vertices. If a Hamiltonian path exists whose endpoints are adjacent, then the resulting graph cycle is called a *Hamiltonian cycle*. Any Hamiltonian cycle can be converted to a Hamiltonian path by removing one of its edges, but a Hamiltonian path can be extended to Hamiltonian cycle only if its endpoints are adjacent.
"(...) Hamiltonian paths and cycles are named after William Rowan Hamilton who invented [in 1857] the icosian game, now also known as Hamilton's puzzle, which involves finding a Hamiltonian cycle in the edge graph of the dodecahedron. Hamilton solved this problem using the icosian calculus, an algebraic structure based on roots of unity with many similarities to the quaternions (also invented by Hamilton). [However,] This solution does not generalize to arbitrary graphs. (...)" (Berge, Claude; Ghouila-Houiri, A. (1962), Programming, games and transportation networks, New York: Sons, Inc.)
These proposed strategies aim to avoid these Hamiltonian cycles and focus on finding this minimum-length path as a solution for the *NP-Complete Hamiltonian Path Problem*. As such, being a NP-Complete problem (`P!=NP`) there is no predefined standard non-exponential-time algorithm to better resolve it, so the current usual way to determine whether a given general graph has a Hamiltonian path is to undertake an exhaustive search. The number of different Hamiltonian cycles in a complete undirected graph on $n$ vertices is $\frac{(n−1)!}{2}$ and in a complete directed graph on n vertices is $(n−1)!$. These counts assume that cycles that are the same apart from their starting point are not counted separately.
The Hamiltonian Path Problem, typical on graph theory, is a special case of the *Travelling Salesman Problem* (also NP-Complete), where the objective is to get the sequence of the minimum-cost route while stopping exactly once per stop and then returning to origin. This last problem, perhaps more worked upon over the years, has known possible solutions as Greedy Algorithms or Dynamic Programming Algorithms, although in neither case the resources to solve them are less than exponential or potential time. 
Another related NP-Complete problem is the *Quadratic Assignment Problem*, fundamental in combinatorial optimization, where given the flow costs between facilites and distance weights between locations, the objective is to assign all facilities to different locations with the goal of minimizing the sum of the distances multiplied by the corresponding flows (e.g. in a custom delivery system). This problem type can be directly formulated in QUBO form and translated to Quantum Annealing Processors, having already been publicly tested for Vehicule Routing and other Logistics problems.

## Problem context
The Routing Optimization Problem proposed in this document consists on finding the minimum-length Hamiltonian path of several routes, with a previously defined best solution based on unknown practical heuristics to which to compare the nodes sequences arrived on the different approaches explicited below. 
The problem input is a selection of information on thousands of routes, splitted on a training and a validation set with similar data. Each route is defined by an identifier, a departure time, a maximum cargo transport dimensions, several unordered stops and tons of packages to deliver to these stops from a unique defined station. The objective is to give a sequence of stops (i.e. nodes) that is as much similar as possible to the optimal already provided.
Then, each stop (both station and dropoff) has also an identifier, as well as coordinates (latitude and longitude) and association to a zone, i.e. an identified non-determined group of stops. Several groups may have stops or zones with the same identifier but belonging to different cities, which makes these identifiers an untrustworthy differentiator. Each stop or zone identifier is only meaningful within the context of each route. However, their coordinates are absolute.
Continuing, each package has its own identifer and is linked to a destination stop, but also has its dimensions (i.e. volume), planned service time (i.e. how much it takes to deliver the package once this reaches destination) and a sometimes specified time window (i.e. a temporal interval in which the package must be delivered). In case the window is not provided, the package can be delivered anytime.
Finally, the cost matrix defined for every combination of stops in a route is the non-symmetric travel times matrix (#Origins, #Destionations), where it is defined the time it takes to go from one stop to another considering usual elements like usual transit, signalization delays and permitted speeds. There is no differentiation regarding dates or time of day.
The training data is composed of all this common information, plus the defined actual stops sequences, routes quality scores and packages delivery status.

## Problem analysis
Upon first procedure diagramation, the problem is divided into three global approximations, each focused on a split of main attributes and its relationships. 

* **1st Approximation**: Using zones definitions and travel times, as well as general route and stop data. The relevant parameters to consider are:
    *  Relations: Stop->Zone, Zone->Stop
    *  Stops clustering: by zone, other
    *  Travel times sequencing: Local, Global
    *  Zones sequencing: Default (as given), Baricenter Euclidian Distance, Hausdorf Distance
    *  Next Stop: Last inner stop (i.e. last stop from previous zone), First outer stop (i.e. first stop from current zone)
*  **2nd Approximation**: Using packages time considerations, these being planned service time and time windows. The relevant parameters to consider are:
    *  Current datetime: Departure Time + Travel Time + Planned Service Time
    *  Time window conditionality (i.e. hierarchical boosts on destination stops whose packages have time windows specified): Zones clustered, Global
*  **3rd Approximation**: Using packages dimensions and transport maximum capacity. The relevant parameters to consider are:
    *  Maximum capacity reached: need to resupply on station.
    *  Cargo carried: divided per zone, max-saturation.

## Local approach
After initial consideration, the first approach taken is to consider local times minimization, that is, checking only the next stop deciding by the immediately next minimum travel time. 

### Method
Iterative algorithm considering the travel times between not yet utilized stops, making decisions on the spot with respect to which is the upcoming minimum-cost stop. 
If two possible destination stops with same cost, pick the first on the list.
If sequence first stop (i.e. station), decide on every other stop besides itself. If middle-sequence stop, decide on every stop besides itself. If sequence last stop, finish iteration.

### Pseudo-code

    Find Station 
    Add Station to Sequence
    Mark Station as Utilized
    while (Sequence Length != Stops Quantity):
        while (Next Stop is Utilized):
            Get Minimum-Cost Next Stop
            Check Next Stop is Utilized
        Add Next Stop to Sequence
        Mark Next Stop as Utilized

### Examples
Suppose a route with 6 stops A, B, C, D, E, F (where C is station) and a cost matrix (travel times) with the following format:
$$\ Cost(TT) = \begin{pmatrix}0 & 4 & 9 & 7 & 3 & 2\\
9 & 0 & 3 & 5 & 1 & 2\\ 2 & 11 & 0 & 6 & 1 & 5\\ 3 & 6 & 7 & 0 & 4 & 10\\ 9 & 8 & 5 & 7 & 0 & 8\\ 4 & 3 & 7 & 5 & 4 & 0\\
\end{pmatrix}$$

The resulting stops sequencing for this route based on Local Approach is C, E, D, A, F, B.

## Brute-force approach

### Method

### Pseudo-code

### Examples

## Ford-Fulkerson approach

### Method

### Pseudo-code

### Examples

## Comparison of the approaches

## References

