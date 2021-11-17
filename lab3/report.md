# Laboratory work № 3.
## Files:
* consistentMain.c - consistent main task;
* parallelMain.c - contains 2 different realizatin of parallel main task (ordinary structure & changed structure);
* parallelMainOther.c - another one realization of parallel main task (master calculate as well + all in main);
* consistent2d.c - consistent 2d task;
* parallel2d.c - contains 2 different realizatin of parallel 2d task (inner & outer);

## Main task.
The task was to remake void compute_solo(double **a) from parallelMain.c into parallel version.
```c
void compute_solo(double **a) {
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = sin(0.00001 * a[i][j]);
        }
    }
}
```
There are no loop dependencies in this task, so we can parallelize the loops as we want (outer, inner, blocks).
So, I decided to implement several case:
* Ordinary - did not change structure of cycles and paralleld on outer loops
* Unioun - united 2 groups of cycles. Paralleled on outer loop as well.
* Other - united 2 groups of cycles as well. All calculation take place in main. Different distribution: by blocks, not by rows. Mater take part in calculation as well. 

### Results
Measurement table:

![](https://github.com/zigal0/MPI/blob/main/lab3/pics/t1.jpg)

SpeedUp graphs:

![](https://github.com/zigal0/MPI/blob/main/lab3/pics/g1.png)

As i said, in this task there are no loop dependencies, so it's quite easy. Let's take a look at another task 2d.

## 2d task.
The task was to remake void compute_solo(double **a) from parallel2d.c into parallel version.
```c
void compute_solo(double **a) {
    int i, j;
    for (i = 0; i < ISIZE; i++) {
        for (j = 0; j < JSIZE; j++) {
            a[i][j] = 10 * i + j;
        }
    }
    for (i = 8; i < ISIZE; i++) {
        for (j = 0; j < JSIZE - 3; j++) {
            a[i][j] = sin(0.00001 * a[i - 8][j + 3]);
        }
    }
}
```
Here we can see loop-carried dependency. Let's look at it closely and define distance vector and direction vector:

![](https://github.com/zigal0/MPI/blob/main/lab3/pics/f.jpg)

So we can parallelize inner cycle with cashing necessary data. And also we can parallelize outer cycle but only for 8 or less processes.
I implemented both parallel variant:
* Method # 1: parallelize inner loop (there are no restriction for number of processes).
* Method # 2: parallelize outer loop (restriction - not more than 8 slave-processes).

### Results
Measurement table:

![](https://github.com/zigal0/MPI/blob/main/lab3/pics/t2.jpg)

SpeedUp graphs:

![](https://github.com/zigal0/MPI/blob/main/lab3/pics/g2.png)

## Conclusion:
Taking into account how many hours I spent on paralleling these tasks and how speedup is low It is not amazing that writing parallel programs is very expensive task. Moreover, It is often not beneficialю
