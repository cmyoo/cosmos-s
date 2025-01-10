/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* FUNCTIONS FOR FIXED MESH REFINENMENT :: BSSN evolution Class of COSMOS_S                              */
/*                                             ver. 1.00           coded by Chulmoon Yoo   2021/11/17    */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos_s.h"
#include <iomanip>
#include <cmath>
#include <fstream>
#include <cstdio>

//time evolution of this layer
void Fmv1::onestep(int btype)
{
	#pragma omp barrier
	set01();		//xxx1=xxx0 
	setv0();		//xxx0=xxxv 
	#pragma omp barrier

	set_zero_r();	//xxxr=0
	#pragma omp barrier
	
	//each step start
	BSSN(1);
	#pragma omp barrier
	
	set_boundary(btype);
	#pragma omp barrier
	
	//checking constraint 
	check_const();
	#pragma omp barrier

	BSSN(2);
	#pragma omp barrier
	set_boundary(btype);
	#pragma omp barrier
	BSSN(3);
	#pragma omp barrier
	set_boundary(btype+1);
	#pragma omp barrier
	BSSN(4);
	#pragma omp barrier
	set_boundary(btype+1);
	#pragma omp barrier
	set_dtpp(get_dtp());			//one before previous time step
	set_dtp(get_dt0());			//previous time step

	if(mrf)
	 ulay->evolve();
	#pragma omp barrier

	//time forward 
	t+=dt0;
	
	return;
}

//refine lower layer using higher layer
void Fmv1::refine_llay()
{
	int k=0;
	int j=0;
	
	int dell=(lui-lli)/2;
	int llli=llay->get_lli();
	
	#pragma omp parallel for
	for(int l=llli;l<=llli+dell-negb;l++)
	{
		for(int i=0;i<nn;i++)
		 llay->set_bv(l,k,j,i)=get_bv(lli+2*(l-llli),k,j,i);
	}
	
	return;
}

void Fmv1::evolve()
{
	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);
	
	dt0=0.5*llay->get_dt0();
	
	onestep(2);
	#pragma omp barrier
	onestep(4);
	#pragma omp barrier
	refine_llay();
	
}

void Fmv1::set_fmr_initial()
{
	//initial parameter setting start
	// cfl=llay->get_cfl();
	// etaa=llay->get_etaa();
	// etab=llay->get_etab();
	// etabb=llay->get_etabb();
	// lambda=llay->get_lambda();
	// dt0=0.5*llay->get_dt0();
	// dtp=dt0;
	// dtpp=dt0;
	// fluidw=llay->get_fluidw();
	// t=llay->get_t();
	// Hb=llay->get_Hb();
	// tini=llay->get_tini();
	// KOep=llay->get_KOep();
	// exg=2*llay->get_exg();
	// kap_MUSCL=llay->get_Mkap();
	// b_minmod=llay->get_b();
	negb=3;
	mrf=false;
	//initial parameter setting end
	
	set_zero_all();
	#pragma omp barrier
	
	//initial setting
	set_luvec();
	#pragma omp barrier
	set_flat();
	#pragma omp barrier
	set_refz();
	#pragma omp barrier
	
	//setting the variables by interpolating the lower layer values
	#pragma omp parallel for
	for(int l=lli;l<=lui+tab;l++)
	{
		int k=0;
		int j=0;
		
		for(int i=0;i<nn;i++)
		{
			if(l%2==0)
			set_bv(l,k,j,i)=llay->get_bv(l/2,k,j,i);
			else
			{
				//set_bv(l,k,j,i)=llay->bv_spline_ipol(int(l/2),get_z(l),i);
				set_bv(l,k,j,i)=llay->bv_ipol(int(l/2),get_z(l),6,i);
			}
		}
	}

	#pragma omp barrier
	set_boundary(5);
	
	#pragma omp barrier
	setv0();
	
	#pragma omp barrier
	dyntoprim();
	
	return;
}

void Fmv1::set_boundary(int itype)
{
	
	double delt=llay->get_dt();
	double deltp=llay->get_dtpp();
	
	double aa,bb,cc;
	
	switch(itype){
	case 1:
		aa=0.;
		bb=1.;
		cc=0.;
		break;
	case 2:
		aa=(-3*pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/16.;
		bb=(3*(delt + 4*deltp)*pow(deltp,-1))/16.;
		cc=((delt + 4*deltp)*pow(delt + deltp,-1))/16.;
		break;
	case 3:
		aa=-(pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/4.;
		bb=((delt + 2*deltp)*pow(deltp,-1))/4.;
		cc=((delt + 2*deltp)*pow(delt + deltp,-1))/4.;
		break;
	case 4:
		aa=(-3*pow(delt,2)*pow(deltp,-1)*pow(delt + deltp,-1))/16.;
		bb=0.25 + (3*delt*pow(deltp,-1))/16.;
		cc=(3*(3*delt + 4*deltp)*pow(delt + deltp,-1))/16.;
		break;
	case 5:
		aa=0.;
		bb=0.;
		cc=1.;
		break;
	default:
		aa=0.;
		bb=1.;
		cc=0.;
	}
	
	#pragma omp parallel for
	for(int l=lui+1-negb*2;l<=lui+tab;l++)
	{
		int k=0;
		int j=0;
		
		for(int i=0;i<nn;i++)
		{
			if(l%2==0)
			set_bv(l,k,j,i)=aa*llay->get_bv1(l/2,k,j,i)+bb*llay->get_bv0(l/2,k,j,i)+cc*llay->get_bv(l/2,k,j,i);
			else
			set_bv(l,k,j,i)=aa*(llay->bv1_ipol(int(l/2),get_z(l),6,i))
					+bb*(llay->bv0_ipol(int(l/2),get_z(l),6,i))
					+cc*(llay->bv_ipol(int(l/2),get_z(l),6,i));
		}
	}
	
	#pragma omp barrier
	boundary_fmr();
	
	return;
}

void Fmv1::boundary_fmr()
{
	short int *refsign;

	// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
	// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
	// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
	
	refsign = new short int[nn];

	for(int i=0;i<nn;i++)
		refsign[i]=1; 

	refsign[3]=-1;
	refsign[6]=-1;
	refsign[23]=-1;

	if(fluidevo)
	refsign[27]=-1;

	for(int l=1;l<=tab;l++)
	{
		for(int i=0;i<nn;i++)
		{
			set_bv(ll[l],0,0,i)=refsign[i]*get_bv(lli+l,0,0,i);
		}
	}
		
	delete[] refsign;
		
	short int ntmpv[ntmpipol];

	ntmpv[0]=13;
	ntmpv[1]=0;
	ntmpv[2]=20;
	ntmpv[3]=24;

	if(scalarevo && fluidevo)
	{
		ntmpv[4]=28;
		ntmpv[5]=29;
		ntmpv[6]=30;
		ntmpv[7]=3;
		ntmpv[8]=6;
		ntmpv[9]=23;
		ntmpv[10]=27;
		ntmpv[11]=7;
		ntmpv[12]=9;
		ntmpv[13]=14;
		ntmpv[14]=16;
	}
	else if(fluidevo)
	{
		ntmpv[4]=28;
		ntmpv[5]=3;
		ntmpv[6]=6;
		ntmpv[7]=23;
		ntmpv[8]=27;
		ntmpv[9]=7;
		ntmpv[10]=9;
		ntmpv[11]=14;
		ntmpv[12]=16;
	}
	else if(scalarevo)
	{
		ntmpv[4]=25;
		ntmpv[5]=3;
		ntmpv[6]=6;
		ntmpv[7]=23;
		ntmpv[8]=7;
		ntmpv[9]=9;
		ntmpv[10]=14;
		ntmpv[11]=16;
	}

	// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
	// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
	// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
	
	//second derivative setting for spline interpolation
	#pragma omp parallel for
	for(int i=0;i<ntmpipol;i++)
		put_ddy(ntmpv[i]);
	//second derivative setting for spline interpolation

	#pragma omp parallel for
	for(int i=0;i<ntmpipol;i++)
	{	
		for(int l=lmin;l<=lui;l++)
		{
			for(int k=kmin;k<=kmax;k++)
			{
				for(int j=jmin;j<=jmax;j++)
				{
					if(k==kli && j==jli)
					continue;
					
					double refz=get_refz(l,k,j);
					
					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
					
					int lc=int(refz*dzi);
					
					if(lc>=lui)
					{
						//cout << "higher " << endl;
						if(lc+ipolbufu>lu[tab])
						lc=lu[tab]-ipolbufu;
						
						set_tmpipol(l,k,j,i)=bv_ipol(lc,refz,ipolo,ntmpv[i]);
					}
					else
					set_tmpipol(l,k,j,i)=bv_spline_ipol(lc,refz,ntmpv[i]);
					
				}
			}
		}
	}
	
	if(fluidevo && scalarevo)
	{
		#pragma omp parallel for
		for(int l=lmin;l<=lui;l++)
		{
			double z=get_z(l);
			double Z=funcf(z);
			double fzz=get_flat_df2z(l);
			double dfzi=1./sqrt(fzz);
			
			for(int k=kmin;k<=kmax;k++)
			{
				double Y=get_y(k);
				
				for(int j=jmin;j<=jmax;j++)
				{
					if(k==kli && j==jli)
					 continue;
					
					double X=get_x(j);
					
					double rh=sqrt(X*X+Y*Y);
					double R=sqrt(rh*rh+Z*Z);
					
					double sith=rh/R;
					double coth=Z/R;
					
					double siph=Y/rh;
					double coph=X/rh;
					
					double refz=get_refz(l,k,j);
					
					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
					
					double dfzref=df(refz);
					
					double psi=get_tmpipol(l,k,j,0);
					double sqgamref=dfzref*exp(6.*psi);
					double sqgami=dfzi*exp(-6.*psi);
					
					//scalar quantities
					set_bv(l,k,j,0)=get_tmpipol(l,k,j,1);
					set_bv(l,k,j,13)=psi;
					set_bv(l,k,j,20)=get_tmpipol(l,k,j,2);
					set_bv(l,k,j,24)=get_tmpipol(l,k,j,3)/(sqgamref*sqgami);
					set_bv(l,k,j,28)=get_tmpipol(l,k,j,4)/(sqgamref*sqgami);
					set_bv(l,k,j,29)=get_tmpipol(l,k,j,5);
					set_bv(l,k,j,30)=get_tmpipol(l,k,j,6);
					
					//vector quantities
					double bZref=get_tmpipol(l,k,j,7)*dfzref;
					set_bv(l,k,j,1)=bZref*sith*coph;
					set_bv(l,k,j,2)=bZref*sith*siph;
					set_bv(l,k,j,3)=bZref*coth*dfzi;

					double bbZref=get_tmpipol(l,k,j,8)*dfzref;
					set_bv(l,k,j,4)=bbZref*sith*coph;
					set_bv(l,k,j,5)=bbZref*sith*siph;
					set_bv(l,k,j,6)=bbZref*coth*dfzi;
					
					double GamZref=get_tmpipol(l,k,j,9)*dfzref;
					set_bv(l,k,j,21)=GamZref*sith*coph;
					set_bv(l,k,j,22)=GamZref*sith*siph;
					set_bv(l,k,j,23)=GamZref*coth*dfzi;
					
					double p_Zref=get_tmpipol(l,k,j,10)/dfzref/sqgamref;
					set_bv(l,k,j,25)=p_Zref*sith*coph/sqgami;
					set_bv(l,k,j,26)=p_Zref*sith*siph/sqgami;
					set_bv(l,k,j,27)=p_Zref*coth/dfzi/sqgami;

					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25

					//tensor quantities
					double gXX=get_tmpipol(l,k,j,11)+1.;
					double gZZ=(get_tmpipol(l,k,j,12))/(dfzref*dfzref)+1.;
					
					set_bv(l,k,j,7)=gZZ*pow(sith*coph,2)+gXX*(pow(coth*coph,2)+siph*siph)-1.;
					set_bv(l,k,j,8)=gZZ*pow(sith*siph,2)+gXX*(pow(coth*siph,2)+coph*coph)-1.;
					set_bv(l,k,j,9)=(gZZ*pow(coth,2)+gXX*pow(sith,2)-1.)*fzz;
					
					set_bv(l,k,j,10)=gZZ*sith*sith*siph*coph+gXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,11)=(gZZ*sith*coth*coph-gXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,12)=(gZZ*sith*coth*siph-gXX*sith*coth*siph)/dfzi;

					double aXX=get_tmpipol(l,k,j,13);
					double aZZ=get_tmpipol(l,k,j,14)/(dfzref*dfzref);
					
					set_bv(l,k,j,14)=aZZ*pow(sith*coph,2)+aXX*(pow(coth*coph,2)+siph*siph);
					set_bv(l,k,j,15)=aZZ*pow(sith*siph,2)+aXX*(pow(coth*siph,2)+coph*coph);
					set_bv(l,k,j,16)=(aZZ*pow(coth,2)+aXX*pow(sith,2))*fzz;
					
					set_bv(l,k,j,17)=aZZ*sith*sith*siph*coph+aXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,18)=(aZZ*sith*coth*coph-aXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,19)=(aZZ*sith*coth*siph-aXX*sith*coth*siph)/dfzi;
				}
			}
		}
	}
	else if(fluidevo)
	{
		#pragma omp parallel for
		for(int l=lmin;l<=lui;l++)
		{
			double z=get_z(l);
			double Z=funcf(z);
			double fzz=get_flat_df2z(l);
			double dfzi=1./sqrt(fzz);
				
			for(int k=kmin;k<=kmax;k++)
			{
				double Y=get_y(k);
				
				for(int j=jmin;j<=jmax;j++)
				{
					if(k==kli && j==jli)
					 continue;
					
					double X=get_x(j);
					
					double rh=sqrt(X*X+Y*Y);
					double R=sqrt(rh*rh+Z*Z);
					
					double sith=rh/R;
					double coth=Z/R;
					
					double siph=Y/rh;
					double coph=X/rh;
					
					double refz=get_refz(l,k,j);
					
					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
					
					double dfzref=df(refz);
					
					double psi=get_tmpipol(l,k,j,0);
					double sqgamref=dfzref*exp(6.*psi);
					double sqgami=dfzi*exp(-6.*psi);
					
					
					//scalar quantities
					set_bv(l,k,j,0)=get_tmpipol(l,k,j,1);
					set_bv(l,k,j,13)=psi;
					set_bv(l,k,j,20)=get_tmpipol(l,k,j,2);
					set_bv(l,k,j,24)=get_tmpipol(l,k,j,3)/(sqgamref*sqgami);
					set_bv(l,k,j,28)=get_tmpipol(l,k,j,4)/(sqgamref*sqgami);
					
					//vector quantities
					double bZref=get_tmpipol(l,k,j,5)*dfzref;
					set_bv(l,k,j,1)=bZref*sith*coph;
					set_bv(l,k,j,2)=bZref*sith*siph;
					set_bv(l,k,j,3)=bZref*coth*dfzi;

					double bbZref=get_tmpipol(l,k,j,6)*dfzref;
					set_bv(l,k,j,4)=bbZref*sith*coph;
					set_bv(l,k,j,5)=bbZref*sith*siph;
					set_bv(l,k,j,6)=bbZref*coth*dfzi;
					
					double GamZref=get_tmpipol(l,k,j,7)*dfzref;
					set_bv(l,k,j,21)=GamZref*sith*coph;
					set_bv(l,k,j,22)=GamZref*sith*siph;
					set_bv(l,k,j,23)=GamZref*coth*dfzi;
					
					double p_Zref=get_tmpipol(l,k,j,8)/dfzref/sqgamref;
					set_bv(l,k,j,25)=p_Zref*sith*coph/sqgami;
					set_bv(l,k,j,26)=p_Zref*sith*siph/sqgami;
					set_bv(l,k,j,27)=p_Zref*coth/dfzi/sqgami;

					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25

					//tensor quantities
					double gXX=get_tmpipol(l,k,j,9)+1.;
					double gZZ=get_tmpipol(l,k,j,10)/(dfzref*dfzref)+1.;
					
					set_bv(l,k,j,7)=gZZ*pow(sith*coph,2)+gXX*(pow(coth*coph,2)+siph*siph)-1.;
					set_bv(l,k,j,8)=gZZ*pow(sith*siph,2)+gXX*(pow(coth*siph,2)+coph*coph)-1.;
					set_bv(l,k,j,9)=(gZZ*pow(coth,2)+gXX*pow(sith,2)-1.)*fzz;
					
					set_bv(l,k,j,10)=gZZ*sith*sith*siph*coph+gXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,11)=(gZZ*sith*coth*coph-gXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,12)=(gZZ*sith*coth*siph-gXX*sith*coth*siph)/dfzi;

					double aXX=get_tmpipol(l,k,j,11);
					double aZZ=get_tmpipol(l,k,j,12)/(dfzref*dfzref);
					
					set_bv(l,k,j,14)=aZZ*pow(sith*coph,2)+aXX*(pow(coth*coph,2)+siph*siph);
					set_bv(l,k,j,15)=aZZ*pow(sith*siph,2)+aXX*(pow(coth*siph,2)+coph*coph);
					set_bv(l,k,j,16)=(aZZ*pow(coth,2)+aXX*pow(sith,2))*fzz;
					
					set_bv(l,k,j,17)=aZZ*sith*sith*siph*coph+aXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,18)=(aZZ*sith*coth*coph-aXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,19)=(aZZ*sith*coth*siph-aXX*sith*coth*siph)/dfzi;
				}
			}
		}
	}
	else if(scalarevo)
	{
		#pragma omp parallel for
		for(int l=lmin;l<=lui;l++)
		{
			double z=get_z(l);
			double Z=funcf(z);
			double fzz=get_flat_df2z(l);
			double dfzi=1./df(z);
			
			for(int k=kmin;k<=kmax;k++)
			{
				double Y=get_y(k);
				
				for(int j=jmin;j<=jmax;j++)
				{
					if(k==0 && j==0)
					 continue;
					
					double X=get_x(j);
					
					double rh=sqrt(X*X+Y*Y);
					double R=sqrt(rh*rh+Z*Z);
					
					double sith=rh/R;
					double coth=Z/R;
					
					double siph=Y/rh;
					double coph=X/rh;
					
					double refz=get_refz(l,k,j);
					
					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(24):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
					
					double dfzref=df(refz);
					
					double psi=get_tmpipol(l,k,j,0);
					
					//scalar quantities
					set_bv(l,k,j,0)=get_tmpipol(l,k,j,1);
					set_bv(l,k,j,13)=psi;
					set_bv(l,k,j,20)=get_tmpipol(l,k,j,2);
					set_bv(l,k,j,24)=get_tmpipol(l,k,j,3);
					set_bv(l,k,j,25)=get_tmpipol(l,k,j,4);
					
					//vector quantities
					double bZref=get_tmpipol(l,k,j,5)*dfzref;
					set_bv(l,k,j,1)=bZref*sith*coph;
					set_bv(l,k,j,2)=bZref*sith*siph;
					set_bv(l,k,j,3)=bZref*coth*dfzi;

					double bbZref=get_tmpipol(l,k,j,6)*dfzref;
					set_bv(l,k,j,4)=bbZref*sith*coph;
					set_bv(l,k,j,5)=bbZref*sith*siph;
					set_bv(l,k,j,6)=bbZref*coth*dfzi;
					
					double GamZref=get_tmpipol(l,k,j,7)*dfzref;
					set_bv(l,k,j,21)=GamZref*sith*coph;
					set_bv(l,k,j,22)=GamZref*sith*siph;
					set_bv(l,k,j,23)=GamZref*coth*dfzi;

					// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
					//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
					//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
					//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
					//  13:wa , -> psi 
					//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
					//  20:ek  , -> trK
					// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
					// 24:E, 25:p_x(S_x), 26:p_y(S_y), 27:p_z(S_z), 28:N_B(D), 
					// 29(24):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25

					//tensor quantities
					double gXX=get_tmpipol(l,k,j,8)+1.;
					double gZZ=(get_tmpipol(l,k,j,9))/(dfzref*dfzref)+1.;
					
					set_bv(l,k,j,7)=gZZ*pow(sith*coph,2)+gXX*(pow(coth*coph,2)+siph*siph)-1.;
					set_bv(l,k,j,8)=gZZ*pow(sith*siph,2)+gXX*(pow(coth*siph,2)+coph*coph)-1.;
					set_bv(l,k,j,9)=(gZZ*pow(coth,2)+gXX*pow(sith,2)-1.)*fzz;
					
					set_bv(l,k,j,10)=gZZ*sith*sith*siph*coph+gXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,11)=(gZZ*sith*coth*coph-gXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,12)=(gZZ*sith*coth*siph-gXX*sith*coth*siph)/dfzi;

					double aXX=get_tmpipol(l,k,j,10);
					double aZZ=get_tmpipol(l,k,j,11)/(dfzref*dfzref);
					
					set_bv(l,k,j,14)=aZZ*pow(sith*coph,2)+aXX*(pow(coth*coph,2)+siph*siph);
					set_bv(l,k,j,15)=aZZ*pow(sith*siph,2)+aXX*(pow(coth*siph,2)+coph*coph);
					set_bv(l,k,j,16)=(aZZ*pow(coth,2)+aXX*pow(sith,2))*fzz;
					
					set_bv(l,k,j,17)=aZZ*sith*sith*siph*coph+aXX*(coth*coth*siph*coph-siph*coph);
					set_bv(l,k,j,18)=(aZZ*sith*coth*coph-aXX*sith*coth*coph)/dfzi;
					set_bv(l,k,j,19)=(aZZ*sith*coth*siph-aXX*sith*coth*siph)/dfzi;
				}
			}
		}
	}
	return;
}

