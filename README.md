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
