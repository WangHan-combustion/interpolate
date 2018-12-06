/*
 * functions.h
 *
 *  Created on: Apr 26, 2017
 *      Author: simone
 */
#ifndef FUNCTIONS_H_
#define FUNCTIONS_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define myfloat float

#define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a >= _b ? _a : _b; })

#define idx_T(i,j,k,im,jm) \
                ({ i + (im+2)*( j + k*(jm+2) ) ; }) //for values going from 0 to max-1, 3-D

class Grid
{
public:
        int imax, jmax, kmax;
        myfloat Lx, Ly, Lz;
        myfloat *x,*y,*z,*dx,dy,dz;
        myfloat fact, *T;

        void allocGrid(int _im, int _jm, int _km, myfloat L1, myfloat L2, myfloat L3, myfloat _fact)
        {
                imax = _im; jmax = _jm; kmax = _km;
                Lx = L1; Ly = L2; Lz = L3;
                fact = _fact;
                x  = new myfloat[_im+1];
                y  = new myfloat[_jm+1];
                z  = new myfloat[_km+1];
                dx = new myfloat[_im+1];
		T  = new myfloat[(_im+2)*(_jm+2)*(_km+2)];               

                for (int i=0; i<(_im+2)*(_jm+2)*(_km+2); i++) T[i] = 0;
        }
        void init_comp()
        {
                myfloat xaux;

                dy = Ly/jmax;
                dz = Lz/kmax;
                for (int i=1; i<imax+1; i++) {
                        xaux  =  1.*(i-0.5)/imax - 0.5;
                        x[i]    =  1.0 + tanh(fact*xaux)/tanh(fact*0.5);
                }
                for (int j=1; j<jmax+1; j++)
                        y[j] = (j-0.5)*dy;
                for (int k=1; k<kmax+1; k++)
                        z[k] = (k-0.5)*dz;
        }
        void init()
        {
                myfloat dx,dxaux,xaux,xmax,xnorm;

      		dx = Lx/imax;
                dy = Ly/jmax;
                dz = Lz/kmax;
	       	
		xmax = Lx;
	 	
		myfloat ru[imax+1];
		ru[0]=0.;


     		for(int i=1; i<=imax/2; i++)
      		{
		        xaux  = 1.*i/imax;
		        dxaux = 0.5-1.45*(xaux-0.5)**2.;
		        ru[i]=ru[i-1]+dxaux;
		}        
		xnorm = ru[imax/2]
	        for(int i=1; i<=imax/2; i++)
		{
		        ru[i]=Lx/2.*ru[i]/xnorm;
		}
                for(int i=imax; i>=imax/2+1; i--)
		{
		        ru[i]=Lx-ru[imax-i];
		}

		for(int i=1; i<=imax; i++)
		{
			x[i]=0.5*(ru[i]+ru[i-1]);
		}

                for (int j=1; j<jmax+1; j++)
                        y[j] = (j-0.5)*dy;
                for (int k=1; k<kmax+1; k++)
                        z[k] = (k-0.5)*dz;
        }
        void printX(int p)
        {

	        FILE *fp;
	        if(p==0) fp = fopen("oldX.txt","w+");
	        if(p==1) fp = fopen("newX.txt","w+");

	        for(int i=1; i<imax+1; i++) {
	            fprintf(fp,"%f %f\n",(1.0*i)/(1.0*imax),x[i]); }
	        fclose(fp);
        }
        void printY(int p)
        {

	        FILE *fp;
	        if(p==0) fp = fopen("oldY.txt","w+");
	        if(p==1) fp = fopen("newY.txt","w+");

	        for(int i=1; i<jmax+1; i++) {
	            fprintf(fp,"%f %f\n",(1.0*i)/(1.0*jmax),y[i]); }
        	fclose(fp);
       	}
        void printZ(int p)
        {

        	FILE *fp;
        	if(p==0) fp = fopen("oldZ.txt","w+");
        	if(p==1) fp = fopen("newZ.txt","w+");

        	for(int i=1; i<kmax+1; i++) {
        	    fprintf(fp,"%f %f\n",(1.0*i)/(1.0*kmax),z[i]); }
        	fclose(fp);
        }
};


//forward function definitions
void interp3D(Grid *n, Grid *o);
void spline(myfloat x[], myfloat y[], int n, myfloat yp1, myfloat ypn, myfloat y2[]);
myfloat splint(myfloat xa[], myfloat ya[], myfloat y2a[], int n, myfloat x);
myfloat *vector12(long nl, long nh);
void free_vector(myfloat *v, long nl, long nh);

#endif /* FUNCTIONS_H_ */
