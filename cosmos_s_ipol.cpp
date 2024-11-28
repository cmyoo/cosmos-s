/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* INITERPOLATION FUNCTIONS :: BSSN evolution Class of COSMOS_S                                          */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/* See "Numerical recipes" for the algorithm of the Lagrange and 3rd order spline interpolation          */ 
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos_s.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

//interpolation by using Lagrange formula
double Fmv0::ipol( double rr,double *xx,double *yy,int order )
{
	//rr:position, xx:grid, yy:variable, order: order of interpolation
	double *kk;
	kk = new double[order];
	double ans=0.;
	for(int i=0;i<order;i++)
	{
		kk[i]=1.;
		for(int n=0;n<order;n++)
		{
			if(n==i) continue;
			kk[i] *= (rr-xx[n])/(xx[i]-xx[n]);
		}
		ans += kk[i]*yy[i];
	}
	delete[] kk;
	return ans;
}

//interpolation for bv
double Fmv0::bv_ipol(int lc,double zc,int order,int number)
//jc,kc,lc:gridnumber of first smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;									// final answer 
	double *xx, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;								// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the pnterpolation
	for(int l=0;l<order;l++)
	{
		xx[l]  =get_z(lc-buf+l);
		//xx[l]  =get_coordZ(lc-buf+l);
		kxyc[l]=get_bv(lc-buf+l,0,0,number);
	}
	////////// preparation of values of variables to be used /////////////////////////

	ans = ipol( zc, xx, kxyc, order );
	//ans = ipol( funcf(zc), xx, kxyc, order );
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	delete[] xx;
	delete[] kxyc;
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for bv0
double Fmv0::bv0_ipol(int lc,double zc,int order,int number)
//jc,kc,lc:gridnumber of first smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;									// final answer 
	double *xx, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;								// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the pnterpolation
	for(int l=0;l<order;l++)
	{
		xx[l]  =get_z(lc-buf+l);
		//xx[l]  =get_coordZ(lc-buf+l);
		kxyc[l]=get_bv0(lc-buf+l,0,0,number);
	}
	////////// preparation of values of variables to be used /////////////////////////

	ans = ipol( zc, xx, kxyc, order );
	//ans = ipol( funcf(zc), xx, kxyc, order );
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	delete[] xx;
	delete[] kxyc;
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//interpolation for bv1
double Fmv0::bv1_ipol(int lc,double zc,int order,int number)
//jc,kc,lc:gridnumber of first smaller point
//xc,yc,zc:coordinate value of the point
//order:order of the interpolation
//number:specifies the variable for interpolation
{
	////////// variables definition and memory allocation ////////////////////////////
	double ans=0.;									// final answer 
	double *xx, *kxyc;
	xx  = new double[order];							// coordinate values of grid points
	kxyc= new double[order];							// interpolated for x and y
	int buf=order/2-1;								// number of shift to smaller grids
	////////// variables definition and memory allocation ////////////////////////////

	////////// preparation of values of variables to be used /////////////////////////
	//in the pnterpolation
	for(int l=0;l<order;l++)
	{
		xx[l]  =get_z(lc-buf+l);
		//xx[l]  =get_coordZ(lc-buf+l);
		kxyc[l]=get_bv1(lc-buf+l,0,0,number);
	}
	////////// preparation of values of variables to be used /////////////////////////

	ans = ipol( zc, xx, kxyc, order );
	//ans = ipol( funcf(zc), xx, kxyc, order );
	////////// z-coordinate interpolation ////////////////////////////////////////////

	////////// memory release ////////////////////////////////////////////////////////
	delete[] xx;
	delete[] kxyc;
	////////// memory release ////////////////////////////////////////////////////////

	return ans;
}

//See "Numerical recipes" for the algorithm of the spline 3rd order interpolation 
//Setting LU vectors for spline interpolation
void Fmv0::set_luvec()
{
	/////////  defs  ////////////////
	// triple diagonal matrix
	// alp:lower diagonal vec
	// bet:upper diagonal vec
	// gam:diagonal vec
	// alp
	
	//itype = 1

	//for variables of odd parity at the center
	set_alp_o(lli+1)=0.;			//y''=0 at the center the matrix starts from lli+1
	set_gam_o(lli+1)=2.*dz/3.;
	set_bet_o(lli+1)=0.25;
	
	for(int l=lli+2;l<lui;l++)
	{
		set_alp_o(l)=dz/6.;
		set_gam_o(l)=2.*dz/3.-get_alp_o(l)*get_bet_o(l-1);
		set_bet_o(l)=dz/(6.*get_gam_o(l));
	}
	
	//gives the value of the first derivative at the boundary
	set_alp_o(lui)=dz/6.;
	set_gam_o(lui)=dz/3.-get_alp_o(lui)*get_bet_o(lui-1);
	
	//for variables of even parity at the center
	set_alp_e(lli)=0.;			//y'=0 at the center
	set_gam_e(lli)=dz/3.;
	set_bet_e(lli)=0.5;
	
	for(int l=lli+1;l<lui;l++)
	{
		set_alp_e(l)=dz/6.;
		set_gam_e(l)=2.*dz/3.-get_alp_e(l)*get_bet_e(l-1);
		set_bet_e(l)=dz/(6.*get_gam_e(l));
	}
	
	//give the value of the first derivative at the boundary
	set_alp_e(lui)=dz/6.;
	set_gam_e(lui)=dz/3.-get_alp_e(lui)*get_bet_e(lui-1);

	return;
}

//Setting y''
void Fmv0::put_ddy(int i)
{
	int j=0;
	int k=0;
	double r;
	
	if(!evod[i])
	{
		//forward substitution
		set_ddy(lli,i)=0.;
		for(int l=lli+1;l<lui;l++)
		{
			r=dzi*(get_bv(l+1,k,j,i)-2.*get_bv(l,k,j,i)+get_bv(l-1,k,j,i));
			set_ddy(l,i)=(r-get_alp_o(l)*get_ddy(l-1,i))/get_gam_o(l);

		}
		
		//give the value of the first derivative
		r=get_f_z(lui,k,j,i)-dzi*(get_bv(lui,k,j,i)-get_bv(lui-1,k,j,i));
		set_ddy(lui,i)=(r-get_alp_o(lui)*get_ddy(lui-1,i))/get_gam_o(lui);
	
		//backsubstitution
		for(int l=lui-1;l>lli;l--)
		{
			set_ddy(l,i)-=(get_bet_o(l)*get_ddy(l+1,i));
		}
	}
	else
	{
		//forward substitution
		//the first derivative is set to zero at the center
		r=dzi*(get_bv(lli+1,k,j,i)-get_bv(lli,k,j,i));
		set_ddy(lli,i)=r/get_gam_e(lli);
		
		for(int l=lli+1;l<lui;l++)
		{
			r=dzi*(get_bv(l+1,k,j,i)-2.*get_bv(l,k,j,i)+get_bv(l-1,k,j,i));
			set_ddy(l,i)=(r-get_alp_e(l)*get_ddy(l-1,i))/get_gam_e(l);
		}
		
		//give the value of the first derivative
		r=get_f_z(lui,k,j,i)-dzi*(get_bv(lui,k,j,i)-get_bv(lui-1,k,j,i));
		set_ddy(lui,i)=(r-get_alp_e(lui)*get_ddy(lui-1,i))/get_gam_e(lui);
	
		//backsubstitution
		for(int l=lui-1;l>=lli;l--)
		{
			set_ddy(l,i)-=(get_bet_e(l)*get_ddy(l+1,i));
		}
	}
	
}

double Fmv0::bv_spline_ipol(int lc,double zc,int number)
{
	double A,B,C,D;
	int j=0;
	int k=0;
	
	A=(get_z(lc+1)-zc)*dzi;
	B=1.-A;
	C=(pow(A,3)-A)*dz*dz/6.;
	D=(pow(B,3)-B)*dz*dz/6.;
	
	double w=A*get_bv(lc,k,j,number)+B*get_bv(lc+1,k,j,number)+C*get_ddy(lc,number)+D*get_ddy(lc+1,number);
	
	return(w);
}	

