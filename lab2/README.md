## Implementation of the integration of the function sin(1 / x)
### Task
it is required to integrate the function sin(1 / x) from leftB to rightB with accuracy e = 10^(-power).

## Scheme
The scheme "3/8" was chosen for this task:

![](https://github.com/zigal0/MPI/blob/main/lab2/results/scheme.png)
 
 And accuracy of this approximation:
 
 ![](https://github.com/zigal0/MPI/blob/main/lab2/results/accuracy.png)
 
 ## Results
 There are 3 implementations of this task:
 * The interval is divided by parts like interval / 2^(NUM_THREADS - 1). This is the fastest.

 ![](https://github.com/zigal0/MPI/blob/main/lab2/results/integral1.png)
 
 * The interval is divided by the intervals between zeros with different numbers of steps. This is the most precise.

 ![](https://github.com/zigal0/MPI/blob/main/lab2/results/integral2.png)
 
 * The interval is divided by the intervals between zeros with the same numbers of steps. This is simplified version of previous implementation.

 ![](https://github.com/zigal0/MPI/blob/main/lab2/results/integral3.png)
 
 As you can see, we get different results for different implementation and I guess the reason is machine accuracy.
 
 ![image](https://user-images.githubusercontent.com/43552875/115744180-0095ae80-a39b-11eb-9322-118fc35017cd.png)
