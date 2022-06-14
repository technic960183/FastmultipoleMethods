# FastmultipoleMethods
The fast multipole methods with N-body self gravity

Brief Description about the FMM and our work:
In order to solve the potential problem in a more efficient way (rather than Direct N body), FMM takes the advantages of tree structure (with same amazing mathematical tricks of spherical harmonic) that reduce the time complexity from $O(N^2)$ to $O(N)$. Our work mainly focuses on building up the program to implement the FMM on the problem of gravitational potential (which can easily extend to electronstatic). At final, we provided the code of FMM for potential problem with High Accuracy and about O(N) scaling which will make it become faster when the number of particles is large. We have also parallelized the code in C++ version with wonderful performance(FIG.1).

The bonus point we done:
- We have extend the code into 3D !

The main code: FMM.cpp   
The parallelized version:FMM_parallel.cpp

Total time_consumed v.s. Number of Threads (with N=2000)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/T_N_thr.png)
FIG.1   
In FIG.1, we can observe that the total time consumed is aroud twice faster when the number of threads is double. This figure shows that the parallelism of this code(FMM_parallel.cpp) works well and can really reduce the total calcuating time. The value of evaluation is confirmed to be matched with the one calculated by direct N body, you may also check it by yourself through turning on the function by making the line 10 #define COMPARE_TO_DIRECT false to be true, which will print fifty particles' data (in both FMM and direct N).   
Perfomance Comparison (Without doing parallel)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/Comparison_performance.png)
FIG.2   
From this figure, we can see that the FMM makes the line time_consumed to N become a straight line, while the direct_N method starts to grow dramatically when the N is large.    

Accuracy Check (N=5000)
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/Value_compare.png)
FIG.3   
The x-axis of the figure states for the value of potential calculated by Direct N body method and the y-axis is for the values evaluated by FMM. The red-dashed line describe the position that the values are exactly the same(i.e. Direct N = FMM). We can clearly see that the result lies close to the red line, which means the FMM made by us provided a well accuracy on potential evaluation.    

Performance scaling of FMM
![Image text](https://github.com/technic960183/FastmultipoleMethods/blob/main/Figure/time_N.png)
FIG.4   
The image is a normalized image of Total time consumed divided by N compares to N, which shows in FMM the wall-time varies about proportional to N. Notice that there is an abruptly increasing line between N=4000 to N=5000. The reason of such difference is highly connected to the level of the tree constructed, the level is 4 under N=5000 and is 5 for 5000 and above in our test.    

