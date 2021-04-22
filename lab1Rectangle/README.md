## Realization of rectangle scheme for transport equation.
### Equation:
The task was the realization consistent and parallel solution for transport equation: 

![](https://github.com/zigal0/MPI/blob/main/theory/Equation.png)

### Scheme
To solve this problem I use rectangle scheme:

![](https://github.com/zigal0/MPI/blob/main/lab1Rectangle/Results/Solution.png)

### Solution
Here you can solution of this equation (timeSteps = 100, spaceSteps = 100):

![](https://github.com/zigal0/MPI/blob/main/theory/Scheme.png)

So there are two realizations: parallel.c and consistent.c. Also you can find realization schemes such as explicit left corner and cross (commented). 
You can find these schemes in folder theory "lab1.pdf".
### Bash
```c
#!/bin/bash
# Add "+s" for output solutions
cd OutPut
rm *
cd ..
for (( i = 1; i <= 10; i++ ))
do
mpirun -np $i ./lab1Rectangle $1
done > OutPut/results.txt
echo "Finished"

```
There is bash script "Start" for running modeling with different quantity of process (form 1 to 10) and you should type quantity of time steps as parameter. Execution time will be outputted in file results.txt in folder OutPut.
If you add in bash script +s then results of modeling will be outputed in OutPut ("solution" + the number of processes + ".csv").

### Research
Several simulations were performed with different number of processes and different number of time steps. 
Here you can see data obtained from these simulations:

![](https://github.com/zigal0/MPI/blob/main/lab1Rectangle/Results/Table.png)

And here there are 2 Charts:

* From 1 process to 5 processes.
![](https://github.com/zigal0/MPI/blob/main/lab1Rectangle/Results/From1To5.png)

* From 6 processes to 10 processes.
![](https://github.com/zigal0/MPI/blob/main/lab1Rectangle/Results/From6To10.png)

As you can see, from some point time of execution time starts to increase. This effect is related with increasing influence of processes communication(MPI_Send and MPI_Recv).


I took equation from this site (there is analytical solution for this transport equation): http://math.phys.msu.ru/data/374/tema5.pdf

![image](https://user-images.githubusercontent.com/43552875/115744180-0095ae80-a39b-11eb-9322-118fc35017cd.png)
