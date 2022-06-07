#include <cmath>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void Constant_Setup();

double *CarToSph(double XYZ[3]);
double *SphToCar(double SPH[3]);

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
int Tran_tree_idx();

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

    // for (int i = 0; i < N; i++)
    //     printf("%d, %d, %d\n", particles_idx[i * 3 * Level + 6], particles_idx[i * 3 * Level + 9], particles_idx[i * 3 * Level + 12]);

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
            if (lv == 1 && id == 2)
            {
                printf("%d %d %d \n", idx[0], idx[1], idx[2]);
                printf("%d \n", offset_layer);
                printf("%d \n", offset_tree);
            }
            for (int i = 0; i < tree_idx_l[offset_layer]; i++)
                if (tree_idx[offset_tree + i] == -1)
                {
                    tree_idx[offset_tree + i] = id;
                    break;
                }
        }

    for (int i = 0; i < N * Level; i++)
    {
        if (tree_idx[i] == -1)
            printf("Q");
    }

    free(particles_loc);
    free(particles_mass);
    free(particles_idx);
    free(tree_idx);
    free(tree_idx_os);
    free(tree_idx_l);
    free(tree_os);

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
    double phi = atan(y / x);
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