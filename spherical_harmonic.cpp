#include <cmath>
#include <iostream>

const int g_k = 3;
double g_mp[g_k] = {1.,1.,2.};
double g_rho[g_k] = {1.,3.,2.};
double alpha[g_k] = {0.5,2,1.3};
double beta[g_k] = {3.1,2,1.5};


double factorial(int n){
    double result = 1.;
    for (int i=1;i<n+1;i++){
        result *= i;
    }
    return result;
}

double Y(int m,int n,double theta,double phi){
    double output = sqrt(factorial(n-abs(m))/factorial(n+abs(m)))\
                    *pow(-1,m)*std::assoc_legendre(n,abs(m),std::cos(theta))\
                    *std::cos(m*phi);
    return output;
}

double M(int m,int n){
    //k states for number of particles
    double sum = 0.;
    for (int i=0;i<g_k;i++){
        //m_p is the mass of particle
        //rho is the distance parameter in spherical coord.
        //alpha is the zenith angle 
        //beta is the azimuth
        sum += g_mp[i]*pow(g_rho[i],n)*Y(-m,n,alpha[i],beta[i]);
    }
    return sum;
}

double potential_approx(int p, double r,double theta, double phi){
    //p is number of boxes
    //r is the distance parameter of the point wanted in spherical coord.
    double pot = 0.;
    for (int n=0;n<p+1;n++){
        for (int m=-n;m<n+1;m++){
            pot+=M(m,n)*Y(m,n,theta,phi)/pow(r,n+1);
        }
    }
    return pot;
}

int main(){
    //std::cout << M(0,0) <<std::endl;
    //std::cout << Y(0,0,0.3,1.2) << std::endl;
    std::cout << potential_approx(20,5,0.5,1.2) <<std::endl;
}