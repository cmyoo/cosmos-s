/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* BOUNDARY CONDITIONS :: BSSN evolution Class of COSMOS_S                                               */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos_s.h"

//Assumptions:
//fluid + scalar system
//Cartesian around the boundary
//scalar field is initially localized -> 0 in the background

//not accurate enough and the constraint violation propagates inward
//the propagation can be easily checked by the momentum constraint violation

void Fmv::asymcond(int l,int i,double bgv)
{
	//xx, yy, rr are assumed to be X,Y,Z in Cartesian coord.
	double R=funcf(get_z(l));

	double dfz=sqrt(get_flat_df2z(l));
	double delr=dz/dfz;
	
	double zc=get_z(l)-delr;
	
	int lc=int(zc*dzi);
	
	if(lc+ipolbufu>lu[tab]-3)
	lc=lu[tab]-ipolbufu-3;

	double Qv=bv_ipol(lc,zc,ipolo,i)-bgv;

	set_bv(l,0,0,i)=bgv+(1.-delr/R)*Qv;
	return;
	
}


void Fmv::asymcond(int l,int i,double bgv1,double bgv,double dtime,int itype)
{
	//xx, yy, rr are assumed to be X,Y,Z in Cartesian coord.
	double R=funcf(get_z(l));

	double rate=dt0/dtp;
	double fac0=0.5*rate+0.5;
	double fac1=-0.5*rate+0.5;

	if((itype == 1) || (itype ==2))
	{
		fac0=0.25*rate+0.5;
		fac1=-0.25*rate+0.5;
	}

	double alpha=fac0*get_bv0(l,0,0,0)+fac1*get_bv1(l,0,0,0);
	double wa=fac0*get_bv0(l,0,0,13)+fac1*get_bv1(l,0,0,13);
	double dfz=sqrt(get_flat_df2z(l));
	double delR=dtime*alpha*exp(-2.*wa);
	double delr=delR/dfz;
	
	double zc=get_z(l)-delr;
	
	int lc=int(zc*dzi);
	
	if(lc+ipolbufu>lu[tab])
	lc=lu[tab]-ipolbufu;

	double Qv=bv1_ipol(lc,zc,ipolo,i)-bgv1;

	set_bv(l,0,0,i)=bgv+(1.-delR/R)*Qv;

	return;
	
}

void Fmv::boundary_asym(int itype)
{
	double dtime;
	double bgv1[nn],bgv[nn];
	
	if(itype==0)
	{
		bgv[0]=1.;

		for(int i=1;i<nn;i++)
		{
			bgv[i]=0.;
		}

		if(fluidevo)
		{
			double psibg,trkbg,rhobg,sqgam;

			psibg=1./(3.*(1.+fluidw))*log(t/tini);
			trkbg=-2./(1.+fluidw)/t;
			sqgam=exp(6.*psibg)*(1.+amp/(1.+amp));
			rhobg=1./(24.*M_PI)*trkbg*trkbg*sqgam;

			bgv[13]=psibg;
			bgv[20]=trkbg;
			bgv[24]=rhobg;
		}
	}
	else
	{	
		if(itype==1 || itype==2)
			dtime=0.5*dt0+dtp;
		else
			dtime=dt0+dtp;

		bgv1[0]=1.;
		bgv[0]=1.;

		for(int i=1;i<nn;i++)
		{
			bgv1[i]=0.;
			bgv[i]=0.;
		}

		if(fluidevo)
		{
			double tt1=t-dtp,tt;

			if(itype==1 || itype==2)
			{
				tt=t+0.5*dt0;
			}
			else
			{
				tt=t+dt0;
			}

			double psibg1,trkbg1,rhobg1,sqgam1;
			double psibg,trkbg,rhobg,sqgam;

			psibg1=1./(3.*(1.+fluidw))*log(tt1/tini);
			trkbg1=-2./(1.+fluidw)/tt1;
			sqgam1=exp(6.*psibg1)*(1.+amp/(1.+amp));
			rhobg1=1./(24.*M_PI)*trkbg1*trkbg1*sqgam1;

			psibg=1./(3.*(1.+fluidw))*log(tt/tini);
			trkbg=-2./(1.+fluidw)/tt;
			sqgam=exp(6.*psibg)*(1.+amp/(1.+amp));
			rhobg=1./(24.*M_PI)*trkbg*trkbg*sqgam;

			bgv1[13]=psibg1;
			bgv[13]=psibg;
			bgv1[20]=trkbg1;
			bgv[20]=trkbg;
			bgv1[24]=rhobg1;
			bgv[24]=rhobg;
		}
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

	if(itype!=0)
	{
		for(int l=1;l<=tab;l++)
		{
			for(int i=0;i<nn;i++)
			{
				asymcond(lu[l],i,bgv1[i],bgv[i],dtime,itype);
				set_bv(ll[l],0,0,i)=refsign[i]*get_bv(lli+l,0,0,i);
			}
		}
	}
	else
	{
		for(int l=1;l<=tab;l++)
		{
			for(int i=0;i<nn;i++)
			{
				asymcond(lu[l],i,bgv[i]);
				set_bv(ll[l],0,0,i)=refsign[i]*get_bv(lli+l,0,0,i);
			}
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
	// 29(24):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
	
	//second derivative setting for spline interpolation
	#pragma omp parallel for
	for(int i=0;i<ntmpipol;i++)
	put_ddy(ntmpv[i]);
	//second derivative setting for spline interpolation

	#pragma omp parallel for
	for(int i=0;i<ntmpipol;i++)
	{	
		for(int l=lmin;l<=lmax;l++)
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
		for(int l=lmin;l<=lmax;l++)
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
		for(int l=lmin;l<=lmax;l++)
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
		for(int l=lmin;l<=lmax;l++)
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

