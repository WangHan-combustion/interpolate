#include <time.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include "functions.h"

using namespace std;


const myfloat pi =  4.0*atan(1.0);

const int p_row= 12;
const int p_col= 2;
const int iold = 168;  
const int jold = 120/p_row;
const int kold = 120/p_col;
const int inew = 240;
const int jnew = 240/p_row;
const int knew = 240/p_col;
const myfloat Lxold = 2.0 ; 
const myfloat Lyold = 2*pi/p_row;
const myfloat Lzold = 4*pi/p_col;
const myfloat Lxnew = 2.0 ;
const myfloat Lynew = 2*pi/p_row;
const myfloat Lznew = 4*pi/p_col;

extern "C" void interpolate_(double *Tin, double *Tout) //, const int iold, const int jold, const int kold, const int inew, const int jnew, const int knew)
{

	/**************************************************************/
	/********** CREATING GRID AND TEMPERATURE FIELD ***************/
	/**************************************************************/

	Grid gOld, gNew;

	gOld.allocGrid(iold,jold,kold,Lxold,Lyold,Lzold,3.0);
	gNew.allocGrid(inew,jnew,knew,Lxnew,Lynew,Lznew,3.0);

        gOld.init();
        gNew.init();

        gOld.printX(0);
        gNew.printX(1);
        gOld.printY(0);
        gNew.printY(1);
        gOld.printZ(0);
        gNew.printZ(1);
      
	for (int k = 0; k < (kold+2); k++)
	{
		for (int j = 0; j < (jold+2); j++)
		{
			for (int i = 0; i < (iold+2); i++)
			{
				gOld.T[idx_T(i,j,k,iold,jold)] = (myfloat) Tin[idx_T(i,j,k,iold,jold)];
			}
		}
	}

	// interpolating temperature and finding new concentration on coarser grid
	interp3D(&gNew, &gOld);

	for (int k = 0; k < (knew+2); k++)
	{
		for (int j = 0; j < (jnew+2); j++)
		{
			for (int i = 0; i < (inew+2); i++)
			{
				Tout[idx_T(i,j,k,inew,jnew)] = (double) gNew.T[idx_T(i,j,k,inew,jnew)];
			}
		}
	}
}

