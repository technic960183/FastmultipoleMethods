//To main repo
#include <cmath>
#include <iostream>
#include <stdlib.h>


double* CarToSph(double XYZ[3]){
   double x = XYZ[0];
   double y = XYZ[1];
   double z = XYZ[2];
   static double sph[3];
   double r = sqrt(pow(x,2)+pow(y,2)+pow(z,2));
   double theta = acos(z/r);
   double phi = atan(y/x);
   sph[0] = r;
   sph[1] = theta;
   sph[2] = phi;
   return  sph;
}
double* SphToCar(double SPH[3]){
    double r = SPH[0];
    double theta = SPH[1];
    double phi = SPH[2];
    static double car[3];
    double x = r*cos(phi)*sin(theta);
    double y = r*sin(phi)*sin(theta);
    double z = r*cos(theta);
    car[0] = x;
    car[1] = y;
    car[2] = z;
    return car;
}

int main(){
    int N = 1000;// Number of particles
    srand(time(NULL));
    //produce particle list
    double *particles_loc = (double *) malloc(N*3*sizeof(double));
    double *particles_mass = (double *) malloc(N*sizeof(double));
    if (particles_loc==NULL){return 1;}
    for (int i=0;i<N;i++){
        int particle_id = 3*i;//(x,y,z) for each
        for (int j=0;j<3;j++){
            particles_loc[particle_id+j] =((double) rand() / (RAND_MAX));
        }
        particles_mass[i] = ((double) rand() / (RAND_MAX));
    }
    //set p
    double eps = pow(10,-2);
    int p = ceil(-log(eps)/log(pow(3,0.5)));
    //Level
    int level = ceil(log2(N)/3)+1;
    //idx_particle
    int *idx_particles = (int *) malloc(N*level*3*sizeof(int));
    for (int i=0;i<N*level*3;i++){idx_particles[i]=-1;}
    for (int i=0;i<N;i++){
        for (int num=0;num<3;num++){
            idx_particles[i*level*3+num] = 0;
        }
        double res_x = particles_loc[i*3];
        double res_y = particles_loc[i*3+1];
        double res_z = particles_loc[i*3+2];
        for (int lv=1;lv<level;lv++){
            int split_x = int (res_x>0.5);
            int split_y = int (res_y>0.5);
            int split_z = int (res_z>0.5);
            int particle_id = i*level*3;
            idx_particles[particle_id+lv*3] = idx_particles[particle_id+(lv-1)*3]*2+split_x;
            idx_particles[particle_id+lv*3+1] = idx_particles[particle_id+(lv-1)*3+1]*2+split_y;
            idx_particles[particle_id+lv*3+2] = idx_particles[particle_id+(lv-1)*3+2]*2+split_z;
            res_x = 2*res_x-split_x;
            res_y = 2*res_y-split_y;
            res_z = 2*res_z-split_z;

        }
    }
    //XYZToL

    free (particles_loc);
    free (particles_mass);
    free (idx_particles);
    return EXIT_SUCCESS;

}