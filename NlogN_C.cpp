#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Constant_Setup();

double *CarToSph(double XYZ[3]);
double *SphToCar(double SPH[3]);
double *CellCenter(int *idx, int level);
int *NeighboursRange(int *center_idx, int center_level);
int *NeighboursChildRange(int *center_idx, int center_level);

double eps;
int N, P, Level;

int *Pow2;

double *particles_loc, *particles_mass;

int *particles_idx;
int Tran_idx_pt(int n, int lv, int axis);

int *tree_os;

int *tree_idx;
int *tree_idx_os;
int *tree_idx_l;
int XYZToL(int x, int y, int z, int lv);
int XYZToL(int *idx, int lv);
int *LToXYZ(int L, int lv);

double *M_tree;
int Tran_M_tree(int L, int lv, int n, int m);

double *potential;

double factorial(int n);
double Y(int m, int n, double theta, double phi);

int main()
{
    printf("Start.\n");
    Constant_Setup();

    // produce particle list
    srand(time(NULL));

    particles_loc = (double *)malloc(N * 3 * sizeof(double));
    particles_mass = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        int particle_id = 3 * i; //(x,y,z) for each
        for (int j = 0; j < 3; j++)
            particles_loc[particle_id + j] = ((double)rand() / (RAND_MAX));
        particles_mass[i] = ((double)rand() / (RAND_MAX));
    }

    // idx_particle
    particles_idx = (int *)malloc(N * Level * 3 * sizeof(int));
    for (int i = 0; i < N; i++)
    {
        for (int k = 0; k < 3; k++)
            particles_idx[Tran_idx_pt(i, 0, k)] = 0;
        double res_x = particles_loc[i * 3];
        double res_y = particles_loc[i * 3 + 1];
        double res_z = particles_loc[i * 3 + 2];
        for (int lv = 1; lv < Level; lv++)
        {
            int split_x = (res_x > 0.5) ? 1 : 0;
            int split_y = (res_y > 0.5) ? 1 : 0;
            int split_z = (res_z > 0.5) ? 1 : 0;
            int particle_id = i * Level * 3;
            particles_idx[particle_id + lv * 3] = particles_idx[particle_id + (lv - 1) * 3] * 2 + split_x;
            particles_idx[particle_id + lv * 3 + 1] = particles_idx[particle_id + (lv - 1) * 3 + 1] * 2 + split_y;
            particles_idx[particle_id + lv * 3 + 2] = particles_idx[particle_id + (lv - 1) * 3 + 2] * 2 + split_z;
            res_x = 2 * res_x - split_x;
            res_y = 2 * res_y - split_y;
            res_z = 2 * res_z - split_z;
        }
    }

    // Offset of the tree layer
    tree_os = (int *)malloc(Level * sizeof(int));
    tree_os[0] = 0;
    for (int lv = 1; lv < Level; lv++)
        tree_os[lv] = tree_os[lv - 1] + Pow2[3 * (lv - 1)];

    // tree_idx
    tree_idx_l = (int *)malloc((Pow2[3 * Level] - 1) / (8 - 1) * sizeof(int));
    for (int i = 0; i < (Pow2[3 * Level] - 1) / (8 - 1); i++)
        tree_idx_l[i] = 0;
    for (int lv = 0; lv < Level; lv++)
        for (int id = 0; id < N; id++)
        {
            int *idx = &particles_idx[Tran_idx_pt(id, lv, 0)]; // ?? particles_idx + Tran_idx_pt(id, lv, 0) ??
            tree_idx_l[tree_os[lv] + XYZToL(idx, lv)]++;
        }

    tree_idx_os = (int *)malloc((Pow2[3 * Level] - 1) / (8 - 1) * sizeof(int));
    tree_idx_os[0] = 0;
    for (int i = 1; i < (Pow2[3 * Level] - 1) / (8 - 1); i++)
        tree_idx_os[i] = tree_idx_os[i - 1] + tree_idx_l[i - 1];

    tree_idx = (int *)malloc(N * Level * sizeof(int));
    for (int i = 0; i < N * Level; i++)
        tree_idx[i] = -1;
    for (int lv = 0; lv < Level; lv++)
        for (int id = 0; id < N; id++)
        {
            int *idx = &particles_idx[Tran_idx_pt(id, lv, 0)]; // ?? particles_idx + Tran_idx_pt(id, lv, 0) ??
            int offset_layer = tree_os[lv] + XYZToL(idx, lv);
            int offset_tree = tree_idx_os[offset_layer];
            for (int i = 0; i < tree_idx_l[offset_layer]; i++)
                if (tree_idx[offset_tree + i] == -1)
                {
                    tree_idx[offset_tree + i] = id;
                    break;
                }
        }

    // M_tree
    M_tree = (double *)malloc((Pow2[3 * Level] - 1) / (8 - 1) * (P + 1) * (P + 2) / 2 * sizeof(double));
    for (int lv = 0; lv < Level; lv++)
        for (int L = 0; L < Pow2[3 * lv]; L++)
        {
            double *center = CellCenter(LToXYZ(L, lv), lv);
            for (int i = 0; i < tree_idx_l[tree_os[lv] + L]; i++)
            {
                int id = tree_idx[tree_idx_os[tree_os[lv] + L] + i];
                double *corr = &particles_loc[3 * id];
                double *deltaX = new double[3];
                for (int k = 0; k < 3; k++)
                    deltaX[k] = corr[k] - center[k];
                double *corr_sph = CarToSph(deltaX);
                for (int n = 0; n <= P; n++)
                    for (int m = 0; m <= n; m++)
                    {
                        double y = Y(m, n, corr_sph[1], corr_sph[2]);
                        M_tree[Tran_M_tree(L, lv, n, m)] += particles_mass[id] * pow(corr_sph[0], n) * y;
                        // cos(m phi) - i sin(m phi)
                    }
            }
        }
    /*for (int lv = 0; lv < Level-1; lv++)
        for (int L = 0; L < Pow2[3 * lv]; L++)
            for (int n = 0; n <= P; n++)
                for (int m = 0; m <= n; m++)
                {
                    printf("%d,%d,%d,%d: %e\n", lv, L, n, m, M_tree[Tran_M_tree(L, lv, n, m)]);
                }*/

    potential = (double *)malloc(N * sizeof(double));
    for (int i = 0; i < N; i++)
        potential[i] = 0;
    for (int id = 0; id < N; id++)
    {
        int lv = Level - 1;
        int *incl = NeighboursRange(&particles_idx[Tran_idx_pt(id, 0, 0)], lv);
        for (int x = incl[0]; x <= incl[3]; x++)
            for (int y = incl[1]; y <= incl[4]; y++)
                for (int z = incl[2]; z <= incl[5]; z++)
                {
                    int L = XYZToL(x, y, z, lv);
                    for (int i = 0; i < tree_idx_l[tree_os[lv] + L]; i++)
                    {
                        int id_src = tree_idx[tree_idx_os[tree_os[lv] + L] + i];
                        if (id_src == id)
                            continue;
                        double *corr_src = &particles_loc[3 * id_src];
                        double *corr = &particles_loc[3 * id];
                        double r = 0;
                        for (int k = 0; k < 3; k++)
                            r += pow(corr_src[k] - corr[k], 2);
                        r = sqrt(r);
                        potential[id] += particles_mass[id_src] / r;
                    }
                }
    }
    for (int id = 0; id < N; id++)
    {
        for (int lv = 2; lv < Level; lv++)
        {
            int *incl = NeighboursChildRange(&particles_idx[Tran_idx_pt(id, 0, 0)], lv - 1);
            int *excl = NeighboursRange(&particles_idx[Tran_idx_pt(id, 0, 0)], lv);
            if(id==0)
            {
                printf("%d: %d,%d,%d,%d\n",lv-1,incl[0],incl[1],incl[3],incl[4]);
                printf("%d: %d,%d,%d,%d\n",lv,excl[0],excl[1],excl[3],excl[4]);
            }
            for (int x = incl[0]; x <= incl[3]; x++)
                for (int y = incl[1]; y <= incl[4]; y++)
                    for (int z = incl[2]; z <= incl[5]; z++)
                    {
                        if (excl[0] <= x && excl[3] >= x && excl[1] <= y && excl[4] >= y && excl[2] <= z && excl[5] >= z)
                            continue;
                        int L = XYZToL(x, y, z, lv);
                        double *center = CellCenter(LToXYZ(L, lv), lv);
                        for (int i = 0; i < tree_idx_l[tree_os[lv] + L]; i++)
                        {
                            int id_src = tree_idx[tree_idx_os[tree_os[lv] + L] + i];
                            double *corr_src = &particles_loc[3 * id_src];
                            double *deltaX = new double[3];
                            for (int k = 0; k < 3; k++)
                                deltaX[k] = corr_src[k] - center[k];
                            double *corr_sph = CarToSph(deltaX);
                            for (int n = 0; n <= P; n++)
                                for (int m = 0; m <= n; m++)
                                {
                                    double y = Y(m, n, corr_sph[1], corr_sph[2]) / pow(corr_sph[0], n + 1);
                                    // cos(m phi) + i sin(m phi)
                                    potential[id] += M_tree[Tran_M_tree(L, lv, n, m)] * y * (m == 0 ? 1 : 2);
                                }
                        }
                    }
        }
    }

    for (int id = 0; id < N; id++)
    {
        printf("%f, ",potential[id]);
    }
    



    free(particles_loc);
    free(particles_mass);
    free(particles_idx);
    free(tree_idx);
    free(tree_idx_os);
    free(tree_idx_l);
    free(tree_os);
    free(M_tree);
    free(potential);

    return EXIT_SUCCESS;
}

void Constant_Setup()
{
    N = 1000; // Number of particles
    eps = pow(10, -2);
    P = ceil(-log(eps) / log(pow(3, 0.5)));
    Level = 5; // ceil(log2(N) / 3) + 1;
    printf("Grid level: %d\n", Level - 1);

    Pow2 = new int[3 * Level + 1];
    Pow2[0] = 1;
    for (int i = 1; i < 3 * Level + 1; i++)
        Pow2[i] = Pow2[i - 1] * 2;
}

double *CarToSph(double XYZ[3])
{
    double x = XYZ[0];
    double y = XYZ[1];
    double z = XYZ[2];
    static double sph[3];
    double r = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
    double theta = acos(z / r);
    double phi = atan2(y, x);
    sph[0] = r;
    sph[1] = theta;
    sph[2] = phi;
    return sph;
}
double *SphToCar(double SPH[3])
{
    double r = SPH[0];
    double theta = SPH[1];
    double phi = SPH[2];
    static double car[3];
    double x = r * cos(phi) * sin(theta);
    double y = r * sin(phi) * sin(theta);
    double z = r * cos(theta);
    car[0] = x;
    car[1] = y;
    car[2] = z;
    return car;
}
double *CellCenter(int *idx, int level)
{
    static double center[3];
    for (int i = 0; i < 3; i++)
        center[i] = idx[i] * pow(2, -level) + pow(2, -level - 1);
    return center;
}
int *NeighboursRange(int *center_idx, int center_level)
{
    static int R[3 * 2];
    for (int k = 0; k < 3; k++)
    {
        R[k] = (0 > center_idx[3 * center_level + k] - 1) ? 0 : center_idx[3 * center_level + k] - 1;
        R[k+3] = (Pow2[center_level] - 1 < center_idx[3 * center_level + k] + 1) ? Pow2[center_level] - 1 : center_idx[3 * center_level + k] + 1;
    }
    return R;
}
int *NeighboursChildRange(int *center_idx, int center_level)
{
    static int R[3 * 2];
    for (int k = 0; k < 3; k++)
    {
        int a = 2 * (center_idx[3 * center_level + k] - 1);
        int b = 2 * (center_idx[3 * center_level + k] + 1) + 1;
        R[k] = (0 > a) ? 0 : a;
        R[k + 3] = (Pow2[center_level + 1] - 1 < b) ? Pow2[center_level + 1] - 1 : b;
    }
    return R;
}

int Tran_idx_pt(int id, int lv, int axis)
{
    return (id * Level + lv) * 3 + axis;
}

int XYZToL(int x, int y, int z, int lv)
{
    return z * Pow2[lv] * Pow2[lv] + y * Pow2[lv] + x;
}
int XYZToL(int *idx, int lv)
{
    return idx[2] * Pow2[lv] * Pow2[lv] + idx[1] * Pow2[lv] + idx[0];
}
int *LToXYZ(int L, int lv)
{
    static int XYZ[3];
    XYZ[2] = L / Pow2[2 * lv];
    XYZ[1] = (L - XYZ[2] * Pow2[2 * lv]) / Pow2[lv];
    XYZ[0] = (L - XYZ[2] * Pow2[2 * lv]) % Pow2[lv];
    return XYZ;
}

int Tran_M_tree(int L, int lv, int n, int m)
{
    return ((tree_os[lv] + L) * (P + 1) * (P + 2) / 2) + n * (n + 1) / 2 + m;
}

double factorial(int n)
{
    double result = 1.;
    for (int i = 1; i < n + 1; i++)
        result *= i;
    return result;
}
double Y(int m, int n, double theta, double phi)
{
    double output = sqrt(factorial(n - abs(m)) / factorial(n + abs(m))) * pow(-1, m);
    return output * std::assoc_legendre(n, abs(m), cos(theta));
}

/*
particle_loc	(N, 3)
particle_mass	(N)

particle_idx	(N, Level, 3)

tree_idx	size = N*level  dynamic offset
tree_idx_os (8**lv)
tree_idx_l	(8**lv)
tree_os	(Level)

M_tree		(8**lv, p+1, p+1, 2) level offset (L_idx, n, m, 2)

potential	(N)
*/
