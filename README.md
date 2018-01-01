# Primal Simplex Algorithm
Primal simplex algorithm implementation on C using LAPACK for linear optimization.

## How to use it

```
$ make
$ ./simplex < (input file) > (output file)
```

## Input file format

```
m n
c_1  c_2  ... c_n
b_1  b_2  ... b_m
a_11 a_12 ... a_1n
a_21 a_22 ... a_2n
...
a_m1 a_m2 ... a_mn
```

Where:
* m = number of restrictions (rows)
* n = number of variables (columns)
* c = cost array
* b = resources array
* a = restrictions coefficients matrix

![](https://latex.codecogs.com/gif.latex?%5C%5C%20A%20%5Cin%20%5Cmathbb%7BR%7D%5E%7Bm%5Ctimes%20n%7D%20%5C%5C%20b%20%5Cin%20%5Cmathbb%7BR%7D%5En%20%5C%5C%20c%2C%20x%20%5Cin%20%5Cmathbb%7BR%7D%5Em%20%5C%5C%20min%7Bf%28x%29%7D%3Dc%5ETx%20%5C%5C%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%20Ax%20%3D%20b%20%5Cend%7Bmatrix%7D%5Cright.)

Example:

The optimization problem:

![](https://latex.codecogs.com/gif.latex?min%7Bf%28x%29%7D%20%3D%20-1x_1%20-%202x_2%20&plus;%200x_3%20&plus;%200x_4%20&plus;%200x_5%20%5C%5C%20%5Cleft%5C%7B%5Cbegin%7Bmatrix%7D%201x_1&plus;1x_2&plus;1x_3&plus;0x_4&plus;0x_5%26%20%3D%20%26%206%20%5C%5C%201x_1-1x_2&plus;0x_3&plus;1x_4&plus;0x_5%26%20%3D%20%26%204%20%5C%5C%20-1x_1&plus;1x_2&plus;0x_3&plus;0x_4&plus;0x_5%26%20%3D%20%26%204%20%5Cend%7Bmatrix%7D%5Cright.)

Input file:

```
 3  5
-1 -2 0 0 0
 6  4 4
 1  1 1 0 0
 1 -1 0 1 0
-1  1 0 0 0
 ```
 
 ### Output file format:
 
 ```
 Matrix A
 1.000  1.000  1.000  0.000  0.000 
 1.000 -1.000  0.000  1.000  0.000 
-1.000  1.000  0.000  0.000  1.000 

 Cost Array
-1.000 -2.000  0.000  0.000  0.000 

 Array b
 6.000  4.000  4.000 

.: Starting Phase I :.

.: Phase I Ended in 0.000000 Seconds with 3 Iterations :.

.: Starting Phase II :.

.: Phase II Ended in 0.000000 Seconds with 1 Iterations :.

                    [ 1.000]
                    [ 5.000]
The solution x_b = [ 0.000] is optimal.
                    [ 8.000]
                    [ 0.000]

The optimal result is f(x_b) = -11.000
 ```
