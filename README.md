# MPI
Realization of left angle (lab1Triangle) and rectangle (lab1Rectangle) schemes for transport equation.
  Theory:
There are some theory for programming with MPI;
  lab1triangle:
Two realizations: parallel.c and consistent.c.
Also there is realization "cross" scheme but only for consistent way + require one step and then add right boundary with left angle scheme (commented).
There is bash script "Start" for run programs with different quantity of process (from 1 to 10) and you should type quantity of time steps as parameter. Execution time will be outputted in File results.txt in OutPut. 
Also you can see results of process modeling in folder Results. There are two graphs and one table there.
First graph for results with quantity of processes from 1 to 5 and second one for results with quantity of processes from 6 to 10.
  lab1Rectangle:
Two realizations: parallel.c and consistent.c.
There is bash script "Start" for run programs with different quantity of process (from 1 to 10) and you should type quantity of time steps as parameter. Execution time will be outputted in File results.txt in OutPut. 
If you add in bash script +s then results of modeling will be outputed in OutPut ("solution" + the number of processes + ".csv").
I took task from this site (there is analytical solution for this transport equation): http://math.phys.msu.ru/data/374/tema5.pdf
