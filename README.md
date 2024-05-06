<h1> Heat Flow in Inhomogeneous Materials</h1>

2024-05-06 by Fritz Crusius<br><br>

Project means to fill a gap in publications solving elliptic PDEs in
steady state without heat generation. Finite-difference examples I found 
assumed homogeneous materials, whereas the challenge in civil engineering
is to model thermal bridge effects of component inhomogeneities.<br><br>
This project serves as as compendium to a German blog article on heat bridges
https://energieoekonom.de/waermeleitung-finite-difference/.

<h2>Minimal System<h2>

minimal_system.py<br><br>
The minimal system example reflects the explanation of setting up a
finite-difference system of linear equations, and solving it for node temperatures.

<h2>Verification System<h2>

For verification of correctness and accuracy of results, a comparison of 
simulated and theoretical values across a wall of 10 cm thickness consisting
of two layers of 4 cm of solid brick and an in between insulation layer of
2 cm polyurethane.<br><br>
This project also demonstrates how heat conductivities reflected in a
finite-difference-grid should be scaled to take account of the actual 
dimensions of a modeled component. Furthermore, the project shows how 
to model the inner room and outer environment air.



