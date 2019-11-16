I wrote the code in c++ and visualized my results in python. 

The following lines compile and run the code. 

./compile.sh
./program

Output from the run, which includes the overall L2 error and a look at the errors at the nodes and at the midpoints of the elements, is written to output.h.

The approximate solution uh and true solution u are evaluated on the domain and the results are stored in the results folder. In this folder, files are named results_(n value)_(order of f).csv. The errors.csv file contains the L2 error for each combination of n values and order of f. 

I then plot my results using the plot_results.py file and the plots are stored in the plots folder. 

1. As visible in the file output.h, for f(x) = c and f(x) = x, the error on the nodes is very close to zero, as we would expect, and the error on the midpoints of the elements is not as small. However for f(x) = x^2, the error on the nodes is similar in order to the error on the midpoints, after some inspection, I determined that in this case the error on the nodes is higher than we would expect because the approximation of f(x) has higher error when f(x) is higher order. When I tried implementing an exact evaluation of f(x) in python using a built in quadrature method, the error on the nodes was again close to zero as expected. 

2. As the problem stated, I did notice errors decreased as n increased.

3. The slope of the graph is about 1.97 when f(x) = x and 2.017 when f(x) = x^2. Because this is a loglog plot. If we let s be the slope, we have that error = C * (1/n) ^ s, where C is some constant that corresponds to the y intercept in the loglog plot. This means that for larger values of n the error is decreasing exponentially by the slope.
