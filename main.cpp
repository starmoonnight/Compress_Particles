#include <iostream>
#include <fstream>
#include <cmath>
#include "nr.h"
//This code is written based on package of Numerical Recipes, detailed information
//about the functions used can be found in the book Numerical Recipes: The Art of
//Scientific Computing

using namespace std;
using namespace NR;


//define globle variables accessed by func and dfunc
extern const int Np=1000, NMAX=3000, Ni=1000;
double rnew;
int n;
int ii[Np][Ni];


// Driver for routine frprmn
// fun calculates the function value at a given double vector p
DP func(Vec_I_DP &p)
{
    int i,j;
    double dij,u;
    double x2,xi,xj,xij;
    double y2,yi,yj,yij;
    double z2,zi,zj,zij;

    u=0.0;
    for (i=0;i<Np-1;i++)
    {
        xi=p[3*i];
        yi=p[3*i+1];
        zi=p[3*i+2];



        for (j=0;j<Np-1;j++)
        {
            if (ii[i][j]!=0)
            {
                xj=p[3*ii[i][j]];
                xij=xi-xj;
                xij=xij-round(xij); // displacement are adjusted because of periodic boundary condition
                x2=xij*xij;

                yj=p[3*ii[i][j]+1];
                yij=yi-yj;
                yij=yij-round(yij);
                y2=yij*yij;

                zj=p[3*ii[i][j]+2];
                zij=zi-zj;
                zij=zij-round(zij);
                z2=zij*zij;

                dij=sqrt(x2+y2+z2)/rnew;

                if (dij<=2.0)
                {
                    u = u+pow((1.0-dij/2),2.5);
                }

            }
        }
    }
    return u;

}
//derivative of function
void dfunc(Vec_I_DP &p, Vec_O_DP &df)
{
    int i,j;
    double dij,u;
    double x2,xi,xj,xij;
    double y2,yi,yj,yij;
    double z2,zi,zj,zij;

    double dx,dy,dz;

    double myval;



    for(i=0;i<n;i++)
         df[i]=0.0;

    for (i=0;i<Np-1;i++)
    {
        xi=p[3*i];
        yi=p[3*i+1];
        zi=p[3*i+2];

        for (j=0;j<Np-1;j++)
        {
            if (ii[i][j]!=0)
            {
                xj=p[3*ii[i][j]];
                xij=xi-xj;
                xij=xij-round(xij); // displacement are adjusted because of periodic boundary condition
                x2=xij*xij;

                yj=p[3*ii[i][j]+1];
                yij=yi-yj;
                yij=yij-round(yij);
                y2=yij*yij;

                zj=p[3*ii[i][j]+2];
                zij=zi-zj;
                zij=zij-round(zij);
                z2=zij*zij;

                dij=sqrt(x2+y2+z2)/rnew;

                if (dij<=2.0)
                {
                    dx=pow((1.0-dij/2.0),1.5)*xij/dij;
                    dy=pow((1.0-dij/2.0),1.5)*yij/dij;
                    dz=pow((1.0-dij/2.0),1.5)*zij/dij;

                    df[3*i]=df[3*i]-5.0/4.0*dx;
                    df[3*i+1]=df[3*i+1]-5.0/4.0*dy;
                    df[3*i+2]=df[3*i+2]-5.0/4.0*dz;

                    df[3*ii[i][j]]=df[3*ii[i][j]]+5.0/4.0*dx;
                    df[3*ii[i][j]+1]=df[3*ii[i][j]+1]+5.0/4.0*dy;
                    df[3*ii[i][j]+2]=df[3*ii[i][j]+2]+5.0/4.0*dz;
                }

            }

        }

    }
}


int main()
{
    int i,iter,j,k;
    double dia, dij, PI;
    double fret, fract,fractset;
    long double ab;

    const double ftol=1.0e-5;
    double r, range;
    double x2, xi, xij, xj, y2, yi, yij, yj;
    double z2, zi, zij, zj;
    double x[Np], y[Np], z[Np];

     Vec_DP p(NMAX),df(NMAX);
    //Vec_DP is class template NRVec for DP, (DP is double),This is noted in nrtypes_nr.h;
    //NRVec is a class template defined in nrutil_nr.h

    Vec_O_DP mydf(NMAX);
    DP myval;

    PI = acos(-1.0);
    n = NMAX;

    fract = 0.90;

    cout.precision(15);

    //Open input configuration file
    ifstream fpin ("tf_hz_p85");

    //Open output configuration file
    ofstream fpout;
    fpout.open("tf_hz_p90");

    //Volume fraction of hard spheres
    fpin>>fractset;
    dia = pow(6.0*fractset/double(Np)/PI,1.0/3.0);
    r = dia/2.0;

    cout << "initial r=" << r <<endl;

    for (i=0;i<Np;i++)
    {
        fpin >> x[i] >> y[i] >> z[i];
    }
    fpin.close();

    cout<<"compressing..."<<endl;
    rnew = r*pow((fract/fractset),1.0/3.0);
    range = 10.0*rnew;
    cout<<"new r="<<rnew<<endl;

    for (i=0;i<Np;i++)
    {
        p[3*i]=x[i];
        p[3*i+1]=y[i];
        p[3*i+2]=z[i];
        xi=x[i];
        yi=y[i];
        zi=z[i];
        k=0;
        for(j=i+1;j<Np;j++)
        {
            xj=x[j];
            xij=xi-xj;
            xij=xij-round(xij);
            x2=xij*xij;

            yj=y[j];
            yij=yi-yj;
            yij=yij-round(yij);
            y2=yij*yij;

            zj=z[j];
            zij=zi-zj;
            zij=zij-round(zij);
            z2=zij*zij;
            dij=sqrt(x2+y2+z2);

            if (dij<range)
            {
                ii[i][k]=j;
                k=k+1;
            }
        }
    }

    cout<<"equilibrating..."<<endl;

    frprmn(p, ftol, iter, fret, func, dfunc);

    fpout.precision(15);
    fpout<<fract<<endl;
    for (i=0;i<Np;i++)
    {
        x[i]=p[3*i];
        y[i]=p[3*i+1];
        z[i]=p[3*i+2];
        fpout<<x[i]<<endl<<y[i]<<endl<<z[i]<<endl;
    }
    fpout.close();

    return 0;
}

