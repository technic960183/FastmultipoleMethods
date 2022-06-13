# FastmultipoleMethods
The fast multipole methods with N-body self gravity

The bonus point we done:
- We have extend the code into 3D !

The main code: FMM.cpp

Total time_consumed v.s. Number of Threads (with N=2000)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/T_N_thr.png)

Perfomance Comparison (Without doing parallel)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/Comparison_performance.png)
From this figure, we can see that the FMM makes the line time_consumed to N become a straight line, while the direct_N method starts to grow dramatically when the N is large. 

Accuracy Check (N=5000)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/Value_compare.png)
The x-axis of the figure states for the value of potential calculated by Direct N body method and the y-axis is for the values evaluated by FMM. The red-dashed line describe the position that the values are exactly the same(i.e. Direct N = FMM). We can clearly see that the result lies close to the red line, which means the FMM made by us provided a well accuracy on potential evaluation.

Performance scaling of FMM
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/time_N.png)
The image is a normalized image of Total time consumed divided by N compares to N, which shows in FMM the wall-time varies about proportional to N. Notice that there is an abruptly increasing line between N=4000 to N=5000. The reason of such difference is highly connected to the level of the tree constructed, the level is 4 under N=5000 and is 5 for 5000 and above in our test.
