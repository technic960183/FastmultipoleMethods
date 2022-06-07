//second commit
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <complex>



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

double XYZToL(int idx[3], int level){
    int n = pow(2,level);
    return idx[2]*pow(n,2)+idx[1]*n+idx[0];
}

int* LToXYZ(int l, int level){
    int n = pow(2,level);
    static int XYZ[3];
    int z =l/(n*n);
    int y = (l-z*n*n)/n;
    int x = (l-z*n*n) % n;
    XYZ[0] = x;
    XYZ[1] = y;
    XYZ[2] = z;
    return XYZ;
}

double* Cellcenter(double idx[3], int level){
    static double center[3];
    for (int i=0;i<3;i++){
        center[i] = idx[i]*pow(2,-level)+pow(2,-level-1);
    }
    return center;
}

/*
Neighbor Range ??
*/

double factorial(int n){
    double result = 1.;
    for (int i=1;i<n+1;i++){
        result *= i;
    }
    return result;
}

std::complex<double> Y(int m,int n,double theta,double phi){
    double output = sqrt(factorial(n-abs(m))/factorial(n+abs(m)))\
                    *pow(-1,m)*std::assoc_legendre(n,abs(m),std::cos(theta))\
                    ;
    std::complex<double> y(output*cos(m*phi),output*sin(m*phi));
    return y;
}
//Need to sum over k
std::complex<double> M(int m,int n,int k,double *mass,double *xyz){
    std::complex<double> sum =(0,0);
    for (int i=0;i<k;i++){
        double pos[3] = {xyz[i*3],xyz[i*3+1],xyz[i*3+2]};
        double *sph = CarToSph(pos);
        std::complex<double> y =Y(-m,n,sph[1],sph[2]);
        sum += mass[i]*pow(sph[0],n)*y;
    }
    return sum;
}
/*
std::complex<double>phi(double *xyz, std::complex<double>M_coeff){
    double pos[3] = {xyz[0],xyz[1],xyz[2]};
    std::complex<double>tot = (0,0);

}
*/

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
    /*
    for (int i=0;i<N;i++){
        for (int j=0;j<3;j++){
            std::cout << particles_loc[i*3+j] <<',';
        }
        std::cout << '\n' << std::endl;
    }
    */
    
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


    free (particles_loc);
    free (particles_mass);
    free (idx_particles);
    return EXIT_SUCCESS;

}