#ifndef _COSMOS_S_H_
#define _COSMOS_S_H_
/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* HEADER FILE :: BSSN evolution Class of COSMOS_S                                                       */
/*                                             ver. 1.00          coded by Chulmoon Yoo                  */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#ifndef M_PI 
#define M_PI 3.14159265358979323846
#endif
#ifndef M_E 
#define M_E 2.71828182845904523536
#endif

using namespace std;

class Fmv0{

protected:
	///////////////  general parameters and variables ////////////////////////
	bool fluidevo,scalarevo;							// on/off for fluid, scalar evolution
	bool curveval;										// on/off for curvature evaluation
	bool hform;											// horizon formation
	bool exc;											// excision
	bool mrf;											// mesh refinenment flag

	int tab;											// tab for refinement boundary for z
	int tabx;											// tab for refinement boundary for x,y
	int jmax,jmin,kmax,kmin,lmax,lmin;					// maximum/minimum number of grid points
	int jui,jli,kui,kli,lui,lli;						// range to evolve
	int nx,ny,nz,nn,nc,nv,ncoef;						// number of grid points, variables and outputs 
	int nt;												// number of energy momentum tensor
	int npr;											// number of  primitive variables of fluids
	int nsc,nscp;										// label number for scalar field variables
	int ntmpipol;										// number of interpolated variables
	int lhm;											// grid point of max constraint violation
	int lmm;											// grid point of max constraint violation
	int lkm;											// grid point of max Kre. inv.
	int lwm;											// grid point of max Weyl inv.
	int ipolo,ipolbufu,ipolbufl;						// order of interpolation
	int exg;											// excision grid number
	int hl;												// value of l just outside the horizon
	int hn;												// number of horizons so far
	int nneck;											// number of necks so far
	int hnmax;

	int layn;											// layer number
	int negb;											// neglecting boundary grids
	int llmin;											// minimum l for constraint check

	double xu,xl, yu,yl, zu,zl;							// maximum/minimum of coordinate length
	double dx,dy,dz,dxi,dyi,dzi,dxi2,dyi2,dzi2,dvol;	// grid intervals and inverse of those
	double dxi4,dyi4,dzi4;								// dx/12,dx/24 for convinience
	double dxi12,dyi12,dzi12,dxi24,dyi24,dzi24;			// dx/12,dx/24 for convinience
	double t,dt,dt0,dtp,dtpp;							// time intervals dtp:previous dtpp:one before previous
	double tmax,cfl,etaa,etab,etabb;					// CFL condition and gauge parameters
	double KOep;										// dissipation term
	double amp;											// amplitude of inhomogeneous grid
	double lambda;										// cosmological constant
	double pi2,pi4,pi8,pi16,pi32;						// 4*pi,8*pi,16*pi
	///////////////  general parameters and variables ////////////////////////
	
	///////////////  constraints and curvature invariants ////////////////////
	double ham,hammax,mom,mommax;						// values of constraint violation
	double dGam;										// deviation of Gamma: constraint
	double dGammax;										// deviation of Gamma: constraint
	double Kremax;										// value of max Kre. inv.
	double Weylmax;										// value of max Weyl inv.
	///////////////  constraints and curvature invariants ////////////////////
	
	///////////////  fluid parameters ////////////////////////////////////////
	double kap_MUSCL,b_minmod;							// parameters for MUSCL interpolation
	double fluidw;										// parameter for fluid EOS
	///////////////  fluid parameters ////////////////////////////////////////

		///////////////  scalar field parameters /////////////////////////////////
	double scalarm;										// scalar filed mass parameter
	///////////////  scalar field parameters /////////////////////////////////

	///////////////  parameters for initial data and settings ////////////////
	double Hb;											// initial Hubble
	double tini;										// initial time
	double tk2;											// initial trK^2
	double boxL;										// boxsize
	///////////////  parameters for initial data and settings ////////////////
	
	///////////////  coordinates and buffers /////////////////////////////////
	int *ju,*jl,*ku,*kl,*lu,*ll;						// from innermost grid to outermost buffer grid
	double *x;											// coordinate x
	double *y;											// coordinate y
	double *z;											// coordinate z
	///////////////  coordinates and buffers /////////////////////////////////
	
	///////////////  boundary flags //////////////////////////////////////////
	int ***bflag;										// boundary flag
	int ***hflag;										// horizon flag
	///////////////  boundary flags //////////////////////////////////////////

	///////////////  dynamical variables /////////////////////////////////////
	double ****bv;										// variables (BSSN,GAUGE,fluid,[scalar])
	double ****dbv;										// variables (BSSN,GAUGE,fluid,[scalar])
	double ****bv0;										// variables (previous step)
	double ****bv1;										// variables (two steps before)
	double ****bvr;										// variables (sum for Rungekutta)
	///////////////  dynamical variables /////////////////////////////////////

	///////////////  temporary strages for ipol //////////////////////////////
	double ****tmpipol;									// strages for interpolated values
	///////////////  temporary strages for ipol //////////////////////////////
	
	///////////////  fluid fluxes and primitive vars /////////////////////////
	double ****flux_x;									// deposit for x-flux at i-1/2
	double ****flux_y;									// deposit for y-flux at i-1/2
	double ****flux_z;									// deposit for z-flux at i-1/2
	double ****primv;									// primitive variables of fluid
	///////////////  fluid fluxes and primitive vars /////////////////////////
	
	///////////////  reference coordinate values of z ////////////////////////
	double ***refz;										// reference z values for all grid points
	///////////////  reference coordinate values of z ////////////////////////
	
	///////////////  geometrical variables for inhomogeneous coordinate //////
//	double ***flat_det;									// flat metric determinant
	double *coordZ;										// assigned Z foordinate value 
	double *flat_df2z;									// flat metric variables
	double *flat_Gamz;									// bar Gamma_uz_zz
	double *flat_dGamz;									// bar delG_z_uz_zz
	///////////////  geometrical variables for inhomogeneous coordinate //////
	
	///////////////  constraints and output storages /////////////////////////
	double ****con;										// constraints, etc...
	double ****outv;									// constraints, etc...
	///////////////  constraints and output storages /////////////////////////

	///////////////  vectors for spline interpolation ////////////////////////
	double *alp_e;										// vector for non diag of L
	double *bet_e;										// vector for non diag of U 
	double *gam_e;										// vector for diag of L
	double *alp_o;										// vector for non diag of L
	double *bet_o;										// vector for non diag of U 
	double *gam_o;										// vector for diag of L
	double **ddy;										// second derivative of variables
	///////////////  vectors for spline interpolation ////////////////////////
	
	///////////////  supplements /////////////////////////////////////////////
	bool *nonzero;										// flag for nonzero v on z-axis
	bool *evod;											// even or odd variables at the center
	double **horis;										// horizon radiuses
	double *neck;										// neck radiuses
	///////////////  supplements /////////////////////////////////////////////
	
public:
	Fmv0(int tabsz,int tabsx,int lupper,
	double xupper,double zupper,double am,bool fld, bool scl, bool cuev)
	{
		///////////////  coordinate settings /////////////////////////////////////////////
		tab=tabsz;																		// 
		tabx=tabsx;
		ju = new int[tabx+1];
		jl = new int[tabx+1];
		ku = new int[tabx+1];
		kl = new int[tabx+1];
		lu = new int[tab+1];
		ll = new int[tab+1];
		for(int n=0;n<=tab;n++)
		{
			lu[n] = lupper + n;
			ll[n] = 0 - n;
		}
		for(int n=0;n<=tabx;n++)
		{
			ju[n] = 0 + n;
			jl[n] = 0 - n;
			ku[n] = 0 + n;
			kl[n] = 0 - n;
		}
		
		jmax=ju[tabx];
		jmin=jl[tabx];
		kmax=ku[tabx];
		kmin=kl[tabx];
		lmax=lu[tab];
		lmin=ll[tab];

		jui=ju[0];
		jli=jl[0];
		kui=ku[0];
		kli=kl[0];
		lui=lu[0];
		lli=ll[0];

		nx=jmax-jmin+1;									// 2*tabx+1
		ny=kmax-kmin+1;									// 2*tabx+1
		nz=lmax-lmin+1;									//
		
		amp=am;
		
		xu=xupper;
		xl=0.;
		yu=xu;
		yl=0.;
		zu=zupper;
		zl=0.;																			//
		///////////////  coordinate settings /////////////////////////////////////////////

		///////////////  equations on/off ////////////////////////////////////////////////
		fluidevo=fld;									// for fluid
		scalarevo=scl;									// for sclar
		curveval=cuev;									// for curvareu evaluation
		hform=false;									// initially no horizon
		exc=false;										// initially no horizon
		///////////////  equations on/off ////////////////////////////////////////////////

		///////////////  number of variables /////////////////////////////////////////////
		nn=24;											// for geometry
		ntmpipol=10;									// interpolated variables for geometry
		
		if(fld)
		{
			nn+=5;										// for fluid
			ntmpipol+=3;
		}
		
		if(scl)
		{
			nsc=nn;
			nscp=nsc+1;
			
			nn+=2;										// for scalar
			ntmpipol+=2;
		}
		else
		{
			nsc=nscp=0;
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

		npr=5;											// for primitive variables of fluid
		// @@@@@@@@@@@@   primitive variables of fluid   @@@@@@@@@@@
		//  0:rho , 1:V^x  ,  2:V^y  ,  3:V^z   ,  4:epsilon

		nc=5;											// for constraints
		// @@@@@@@@@@@@   constraints   @@@@@@@@@@@
		//  0:normalized ham, 1:ham, 2:normalized momz, 3:momz, 4:Gamma^z-D_gamma^z
		
		nv=12;											// for output variables
		///////////////  number of variables /////////////////////////////////////////////

		///////////////  parameters and temporary variables //////////////////////////////
		ipolo=6;										// order of interpolation
		ipolbufu=ipolo-(ipolo/2);						// buffer number of grids for interpolation
		ipolbufl=ipolo/2-1;								// buffer number of grids for interpolation

		pi2=2.*M_PI;									// 2pi
		pi4=4.*M_PI;									// 4pi
		pi8=8.*M_PI;									// 8pi
		pi16=16.*M_PI;									//16pi
		pi32=32.*M_PI;									//16pi
		
		tk2=1.;											// temporary initial value of trK^2

		t=dt=dt0=dtp=dtpp=tmax=0.;						// temporary initial values wrt time
		cfl=0.;											// temporary initial values for cfl
		
		boxL=1.;										// box size
		negb=0;											// neglecting boundary grids
		llmin=lli;										// minimum l for constraint check
		hl=0;											// value of l just outside the horizon
		hn=0;											// number of horizons
		nneck=0;
		hnmax=10;										// max number of horizons
		///////////////  parameters and temporary variables //////////////////////////////

		///////////////  spatial coordinate intervals ////////////////////////////////////
		dx=(xu-xl)/(double(ju[tabx]-jl[0])); 			// ju[0] :: inner boundary		//
		dy=(yu-yl)/(double(ku[tabx]-kl[0]));
		dz=(zu-zl)/(double(lu[0]-ll[0]));
		
		//////////////////////////////////////////////////
		//boundary positions
		//     0 1 2 3 4 5 6 7 8 9 
		//     B B B I I I I B B B --> inner and buffer
		//           |     |       --> boundary position
		//////////////////////////////////////////////////

		dxi=1./dx;
		dyi=1./dy;
		dzi=1./dz;

		dxi2=0.5*dxi;
		dyi2=0.5*dyi;
		dzi2=0.5*dzi;

		dxi4=0.5*dxi2;
		dyi4=0.5*dyi2;
		dzi4=0.5*dzi2;
		dxi12=dxi/12.;
		dyi12=dyi/12.;
		dzi12=dzi/12.;
		dxi24=dxi/24.;
		dyi24=dyi/24.;
		dzi24=dzi/24.;

		dvol=dx*dy*dz;									// volume						//
		///////////////  spatial coordinate intervals ////////////////////////////////////

		///////////////  coordinate strages //////////////////////////////////////////////
		x = new double[nx];																//
		for(int j=jmin;j<=jmax;j++){
			x[j-jmin]=xl +dx*(double(j-jl[0]));
		}
		y = new double[ny];
		for(int k=kmin;k<=kmax;k++){
			y[k-kmin]=yl +dy*(double(k-kl[0]));
		}
		z = new double[nz];
		for(int l=lmin;l<=lmax;l++){
			z[l-lmin]=zl +dz*(double(l-ll[0]));
		}																				//
		///////////////  coordinates strages /////////////////////////////////////////////
		
		///////////////  variable strages ////////////////////////////////////////////////
		bv = new double***[nn];															//
		dbv = new double***[nn];
		bv0 = new double***[nn];
		bv1 = new double***[nn];
		bvr = new double***[nn];
		
		for(int l=0;l<nn;l++){
			bv[l] = new double**[nz];
			dbv[l] = new double**[nz];
			bv0[l] = new double**[nz];
			bv1[l] = new double**[nz];
			bvr[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				bv[l][k] = new double*[ny];
				dbv[l][k] = new double*[ny];
				bv0[l][k] = new double*[ny];
				bv1[l][k] = new double*[ny];
				bvr[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					bv[l][k][j]= new double[nx];
					dbv[l][k][j]= new double[nx];
					bv0[l][k][j]= new double[nx];
					bv1[l][k][j]= new double[nx];
					bvr[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  variable strages ////////////////////////////////////////////////

		///////////////  strages for interpolated values /////////////////////////////////
		tmpipol = new double***[ntmpipol];												//
		
		for(int l=0;l<ntmpipol;l++){
			tmpipol[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				tmpipol[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					tmpipol[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  strages for interpolated values /////////////////////////////////
		
		///////////////  boundary flags //////////////////////////////////////////////////
		bflag = new int**[nz];															//
		hflag = new int**[nz];
		for(int l=0;l<nz;l++){
			bflag[l] = new int*[ny];
			hflag[l] = new int*[ny];
			for(int k=0;k<ny;k++){
				bflag[l][k] = new int[nx];
				hflag[l][k] = new int[nx];
			}
		}																				//
		///////////////  boundary flags //////////////////////////////////////////////////

		///////////////  refz variables //////////////////////////////////////////////////
		refz = new double**[nz];														//
		for(int l=0;l<nz;l++)
		{
			refz[l] = new double*[ny];
			for(int k=0;k<ny;k++)
			{
				refz[l][k]= new double[nx];
			}
		}																				//
		///////////////  refz variables //////////////////////////////////////////////////

		///////////////  fluid primitive and fluxes //////////////////////////////////////
		primv = new double***[npr];														//
		for(int l=0;l<npr;l++){
			primv[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				primv[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					primv[l][k][j]= new double[nx];
				}
			}
		}
		
		flux_x = new double***[npr-2];
		flux_y = new double***[npr-2];
		flux_z = new double***[npr-2];
		for(int l=0;l<npr-2;l++){
			flux_x[l] = new double**[nz];
			flux_y[l] = new double**[nz];
			flux_z[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				flux_x[l][k] = new double*[ny];
				flux_y[l][k] = new double*[ny];
				flux_z[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					flux_x[l][k][j]= new double[nx];
					flux_y[l][k][j]= new double[nx];
					flux_z[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  fluid primitive and fluxes //////////////////////////////////////

		///////////////  constraints and output variables ////////////////////////////////
		con = new double***[nc];														//
		for(int l=0;l<nc;l++){
			con[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				con[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					con[l][k][j]= new double[nx];
				}
			}
		}

		outv = new double***[nv];
		for(int l=0;l<nv;l++){
			outv[l] = new double**[nz];
			for(int k=0;k<nz;k++){
				outv[l][k] = new double*[ny];
				for(int j=0;j<ny;j++){
					outv[l][k][j]= new double[nx];
				}
			}
		}																				//
		///////////////  constraints and output variables ////////////////////////////////


		///////////////  geometrical variables for inhomogeneous coordinate //////////////
		flat_df2z = new double[nz];														//
		coordZ = new double[nz];
		flat_Gamz = new double[nz];
		flat_dGamz = new double[nz];													//
		///////////////  geometrical variables for inhomogeneous coordinate //////////////
		
		///////////////  vectors for spline interpolation ////////////////////////////////
		alp_e = new double[nz];															//
		bet_e = new double[nz];
		gam_e = new double[nz];
		alp_o = new double[nz];
		bet_o = new double[nz];
		gam_o = new double[nz];
	
		ddy = new double*[nn];
		for(int l=0;l<nn;l++)
		{
			ddy[l] = new double[nz];
		}																				//
		///////////////  vectors for spline interpolation ////////////////////////////////

		///////////////  horizon radius strage ///////////////////////////////////////////
		horis = new double*[hnmax];
		neck = new double[hnmax];
		for(int i=0;i<hnmax;i++)
		{
			horis[i] = new double[2];
		}																				//
		///////////////  horizon radius strage ///////////////////////////////////////////

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

		///////////////  bool flags //////////////////////////////////////////////////////
		nonzero = new bool[nn];															//
		evod = new bool[nn];
		
		nonzero[0]=true;	//0
		evod[0]=true;
		nonzero[1]=false;	//1
		evod[1]=true;
		nonzero[2]=false;	//2
		evod[2]=true;
		nonzero[3]=true;	//3
		evod[3]=false;
		nonzero[4]=false;	//4
		evod[4]=true;
		nonzero[5]=false;	//5
		evod[5]=true;
		nonzero[6]=true;	//6
		evod[6]=false;
		nonzero[7]=true;	//7
		evod[7]=true;
		nonzero[8]=true;	//8
		evod[8]=true;
		nonzero[9]=true;	//9
		evod[9]=true;
		nonzero[10]=false;	//10
		evod[10]=true;
		nonzero[11]=false;	//11
		evod[11]=false;
		nonzero[12]=false;	//12
		evod[12]=false;
		nonzero[13]=true;	//13
		evod[13]=true;
		nonzero[14]=true;	//14
		evod[14]=true;
		nonzero[15]=true;	//15
		evod[15]=true;
		nonzero[16]=true;	//16
		evod[16]=true;
		nonzero[17]=false;	//17
		evod[17]=true;
		nonzero[18]=false;	//18
		evod[18]=false;
		nonzero[19]=false;	//19
		evod[19]=false;
		nonzero[20]=true;	//20
		evod[20]=true;
		nonzero[21]=false;	//21
		evod[21]=true;
		nonzero[22]=false;	//22
		evod[22]=true;
		nonzero[23]=true;	//23
		evod[23]=false;
		
		if(fld)
		{
			nonzero[24]=true;
			evod[24]=true;
			nonzero[25]=false;
			evod[25]=true;
			nonzero[26]=false;
			evod[26]=true;
			nonzero[27]=true;
			evod[27]=false;
			nonzero[28]=true;
			evod[28]=true;
		
			if(scl)
			{
				nonzero[29]=true;
				evod[29]=true;
				nonzero[30]=true;
				evod[30]=true;
			}
		}
		else if(scl)
		{
			nonzero[24]=true;
			evod[24]=true;
			nonzero[25]=true;
			evod[25]=true;
		}																				//
		///////////////  bool flags //////////////////////////////////////////////////////

		///////////////  horizon radius setting //////////////////////////////////////////
		for(int h=0;h<hnmax;h++)
		{
			neck[h]=10.;
			for(int i=0;i<2;i++)
			{
				horis[h][i]=10.;
			}
		}																				//
		///////////////  horizon radius setting //////////////////////////////////////////
	
		set_zero_all();

		for(int l=lmin;l<=lmax;l++){
			for(int k=kmin;k<=kmax;k++){
				for(int j=jmin;j<=jmax;j++){
					set_bflag(l,k,j)=0;
						set_hflag(l,k,j)=0;
				}
			}
		}

		cout.setf(ios_base::fixed, ios_base::floatfield);
		cout.precision(3);
		cout << "          xl : xu  / jli : jui      " << setw(9) << right << xl << " : " << setw(9) << right << xu
		<< "   /   " << setw(4) << right << jli << " : " << setw(4) << right << jui << endl
		<< "          yl : yu  / kli : kui      " << setw(9) << right << yl << " : " << setw(9) << right << yu
		<< "   /   " << setw(4) << right << kli << " : " << setw(4) << right << kui << endl
		<< "          zl : zu  / lli : lui      " << setw(9) << right << zl << " : " << setw(9) << right << zu
		<< "   /   " << setw(4) << right << lli << " : " << setw(4) << right << lui << endl
		<< endl;
	}

	virtual ~Fmv0(){

		delete[] ju;
		delete[] jl;
		delete[] ku;
		delete[] kl;
		delete[] lu;
		delete[] ll;

		delete[] x;
		delete[] y;
		delete[] z;
		
		delete[] flat_df2z;
		delete[] coordZ;
		delete[] flat_Gamz;
		delete[] flat_dGamz;
		
		delete[] alp_e;
		delete[] bet_e;
		delete[] gam_e;
		delete[] alp_o;
		delete[] bet_o;
		delete[] gam_o;
		
		delete[] neck;
	
		for(int l=0;l<hnmax;l++){
			delete[] horis[l];
			horis[l]=NULL;
		}
		delete[] horis;
		horis=NULL;
		
		for(int l=0;l<nn;l++){
			delete[] ddy[l];
			ddy[l]=NULL;
		}
		delete[] ddy;
		ddy=NULL;
		
		for(int l=0;l<nn;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] bv[l][k][j];
					bv[l][k][j]=NULL;
					delete[] dbv[l][k][j];
					dbv[l][k][j]=NULL;
					delete[] bv0[l][k][j];
					bv0[l][k][j]=NULL;
					delete[] bv1[l][k][j];
					bv1[l][k][j]=NULL;
					delete[] bvr[l][k][j];
					bvr[l][k][j]=NULL;
				}
				delete[] bv[l][k];
				bv[l][k]=NULL;
				delete[] dbv[l][k];
				dbv[l][k]=NULL;
				delete[] bv0[l][k];
				bv0[l][k]=NULL;
				delete[] bv1[l][k];
				bv1[l][k]=NULL;
				delete[] bvr[l][k];
				bvr[l][k]=NULL;
			}
			delete[] bv[l];
			bv[l]=NULL;
			delete[] dbv[l];
			dbv[l]=NULL;
			delete[] bv0[l];
			bv0[l]=NULL;
			delete[] bv1[l];
			bv1[l]=NULL;
			delete[] bvr[l];
			bvr[l]=NULL;
		}
		delete[] bv;
		bv=NULL;
		delete[] dbv;
		dbv=NULL;
		delete[] bv0;
		bv0=NULL;
		delete[] bv1;
		bv1=NULL;
		delete[] bvr;
		bvr=NULL;

		for(int l=0;l<ntmpipol;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] tmpipol[l][k][j];
					tmpipol[l][k][j]=NULL;
				}
				delete[] tmpipol[l][k];
				tmpipol[l][k]=NULL;
				}
			delete[] tmpipol[l];
			tmpipol[l]=NULL;
		}
		delete[] tmpipol;
		tmpipol=NULL;
		
		for(int l=0;l<nc;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] con[l][k][j];
					con[l][k][j]=NULL;
				}
				delete[] con[l][k];
				con[l][k]=NULL;
			}
			delete[] con[l];
			con[l]=NULL;
		}
		delete[] con;
		con=NULL;

		for(int l=0;l<nv;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] outv[l][k][j];
					outv[l][k][j]=NULL;
				}
				delete[] outv[l][k];
				outv[l][k]=NULL;
			}
			delete[] outv[l];
			outv[l]=NULL;
		}
		delete[] outv;
		outv=NULL;

		for(int k=0;k<nz;k++){
			for(int j=0;j<ny;j++){
				delete[] refz[k][j];
				refz[k][j]=NULL;
			}
			delete[] refz[k];
			refz[k]=NULL;
		}
		delete[] refz;
		refz=NULL;

		for(int l=0;l<npr;l++){
			for(int k=0;k<nz;k++){
				for(int j=0;j<ny;j++){
					delete[] primv[l][k][j];
					primv[l][k][j]=NULL;
				}
				delete[] primv[l][k];
				primv[l][k]=NULL;
			}
			delete[] primv[l];
			primv[l]=NULL;
		}
		delete[] primv;
		primv=NULL;

		for(int l=0;l<npr-2;l++){
			for(int k=0;k<nz;k++){
				delete[] flux_x[l][k];
				flux_x[l][k]=NULL;
			}
			delete[] flux_x[l];
			flux_x[l]=NULL;
		}
		delete[] flux_x;
		flux_x=NULL;
		
		for(int l=0;l<npr-2;l++){
			for(int k=0;k<nz;k++){
				delete[] flux_y[l][k];
				flux_y[l][k]=NULL;
			}
			delete[] flux_y[l];
			flux_y[l]=NULL;
		}
		delete[] flux_y;
		flux_y=NULL;

		for(int l=0;l<npr-2;l++){
			for(int k=0;k<ny;k++){
				delete[] flux_z[l][k];
				flux_z[l][k]=NULL;
			}
			delete[] flux_z[l];
			flux_z[l]=NULL;
		}
		delete[] flux_z;
		flux_z=NULL;
		
	}

	public:

	////////////////////////////////////////
	//  GET func.
	////////////////////////////////////////
	bool get_fluidevo() const{
		return fluidevo;
	}
	bool get_hform() const{
		return hform;
	}
	bool get_exc() const{
		return exc;
	}
	bool get_mrf() const{
		return mrf;
	}
	int get_jmax() const{
		return jmax;
	}
	int get_jmin() const{
		return jmin;
	}
	int get_kmax() const{
		return kmax;
	}
	int get_kmin() const{
		return kmin;
	}
	int get_lmax() const{
		return lmax;
	}
	int get_lmin() const{
		return lmin;
	}
	int get_jui() const{
		return jui;
	}
	int get_jli() const{
		return jli;
	}
	int get_kui() const{
		return kui;
	}
	int get_kli() const{
		return kli;
	}
	int get_lui() const{
		return lui;
	}
	int get_lli() const{
		return lli;
	}
	int get_nx() const{
		return nx;
	}
	int get_ny() const{
		return ny;
	}
	int get_nz() const{
		return nz;
	}
	int get_lhm() const{
		return lhm;
	}
	int get_lmm() const{
		return lmm;
	}
	int get_lkm() const{
		return lkm;
	}
	int get_lwm() const{
		return lwm;
	}
	int get_exg() const{
		return exg;
	}
	int get_hl() const{
		return hl;
	}
	
	int get_layn() const{
		return layn;
	}
	int get_llmin() const{
		return llmin;
	}

	double get_t() const{
		return t;
	}
	double get_dt() const{
		return dt;
	}
	double get_dt0() const{
		return dt0;
	}
	double get_dtp() const{
		return dtp;
	}
	double get_dtpp() const{
		return dtpp;
	}
	double get_tmax() const{
		return tmax;
	}
	double get_cfl() const{
		return cfl;
	}
	double get_fluidw() const{
		return fluidw;
	}
		double get_scalarm() const{
		return scalarm;
	}
	double get_dx() const{
		return dx;
	}
	double get_dy() const{
		return dy;
	}
	double get_dz() const{
		return dz;
	}
	double get_dvol() const{
		return dvol;
	}
	double get_dxi() const{
		return dxi;
	}
	double get_dyi() const{
		return dyi;
	}
	double get_dzi() const{
		return dzi;
	}
	double get_dxi2() const{
		return dxi2;
	}
	double get_dyi2() const{
		return dyi2;
	}
	double get_dzi2() const{
		return dzi2;
	}

	double get_dxi4() const{
		return dxi4;
	}
	double get_dyi4() const{
		return dyi4;
	}
	double get_dzi4() const{
		return dzi4;
	}

	double get_xu() const{
		return xu;
	}
	double get_xl() const{
		return xl;
	}
	double get_yu() const{
		return yu;
	}
	double get_yl() const{
		return yl;
	}
	double get_zu() const{
		return zu;
	}
	double get_zl() const{
		return zl;
	}
	double get_ham() const{
		return ham;
	}
	double get_hammax() const{
		return hammax;
	}
	double get_Kremax() const{
		return Kremax;
	}
	double get_Weylmax() const{
		return Weylmax;
	}
	double get_mom() const{
		return mom;
	}
	double get_mommax() const{
		return mommax;
	}
	double get_etaa() const{
		return etaa;
	}
	double get_etab() const{
		return etab;
	}
	double get_etabb() const{
		return etabb;
	}
	double get_lambda() const{
		return lambda;
	}
	double get_tini() const{
		return tini;
	}
	double get_KOep() const{
		return KOep;
	}
	double get_Mkap() const{
		return kap_MUSCL;
	}
	double get_b() const{
		return b_minmod;
	}
	double get_Hb() const{
		return Hb;
	}

	int get_bflag(int l,int k,int j) const{
		return bflag[l-lmin][k-kmin][j-jmin];
	}
	int get_hflag(int l,int k,int j) const{
		return hflag[l-lmin][k-kmin][j-jmin];
	}

	double get_x(int j)const{
		return x[j-jmin];
	}
	double get_y(int k)const{
		return y[k-kmin];
	}
	double get_z(int l)const{
		return z[l-lmin];
	}
	double get_ext_x(int j)const{
		return(xl +dx*(double(j-jl[0])));
	}
	double get_ext_y(int k)const{
		return(yl +dy*(double(k-kl[0])));
	}
	double get_ext_z(int l)const{
		return(zl +dz*(double(l-ll[0])));
	}	
	double get_bv(int l,int k,int j,int i) const{
		return bv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_dbv(int l,int k,int j,int i) const{
		return dbv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bv0(int l,int k,int j,int i) const{
		return bv0[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bv1(int l,int k,int j,int i) const{
		return bv1[i][l-lmin][k-kmin][j-jmin];
	}
	double get_bvr(int l,int k,int j,int i) const{
		return bvr[i][l-lmin][k-kmin][j-jmin];
	}
	double get_tmpipol(int l,int k,int j,int i) const{
		return tmpipol[i][l-lmin][k-kmin][j-jmin];
	}
	double get_con(int l,int k,int j,int i) const{
		return con[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flat_df2z(int l)const{
		return flat_df2z[l-lmin];
	}
	double get_coordZ(int l)const{
		return coordZ[l-lmin];
	}
	double get_flat_Gamz(int l)const{
		return flat_Gamz[l-lmin];
	}
	double get_flat_dGamz(int l)const{
		return flat_dGamz[l-lmin];
	}
	double get_alp_e(int l) const{
		return alp_e[l-lmin];
	}
	double get_bet_e(int l) const{
		return bet_e[l-lmin];
	}
	double get_gam_e(int l) const{
		return gam_e[l-lmin];
	}
	double get_alp_o(int l) const{
		return alp_o[l-lmin];
	}
	double get_bet_o(int l) const{
		return bet_o[l-lmin];
	}
	double get_gam_o(int l) const{
		return gam_o[l-lmin];
	}
	double get_ddy(int l,int i) const{
		return ddy[i][l-lmin];
	}
	double get_outv(int l,int k,int j,int i) const{
		return outv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_refz(int l,int k,int j) const{
		return refz[l-lmin][k-kmin][j-jmin];
	}
	double get_primv(int l,int k,int j,int i) const{
		return primv[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_x(int l,int k,int j,int i) const{
		return flux_x[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_y(int l,int k,int j,int i) const{
		return flux_y[i][l-lmin][k-kmin][j-jmin];
	}
	double get_flux_z(int l,int k,int j,int i) const{
		return flux_z[i][l-lmin][k-kmin][j-jmin];
	}
	double*** get_bv(int i) const{
		return bv[i];
	}
	double*** get_dbv(int i) const{
		return dbv[i];
	}

	////////////////////////////////////////
	//  GET derivative func.
	////////////////////////////////////////
	double get_f_x(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l,k,j+2,i)-get_bv(l,k,j-2,i))
			+8.*(get_bv(l,k,j+1,i)-get_bv(l,k,j-1,i))
			)*dxi12;
		return w;
	}
	double get_f_y(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l,k+2,j,i)-get_bv(l,k-2,j,i))
			+8.*(get_bv(l,k+1,j,i)-get_bv(l,k-1,j,i))
			)*dyi12;
		return w;
	}
	double get_f_z(int l,int k,int j,int i) const{
		double w=(   -(get_bv(l+2,k,j,i)-get_bv(l-2,k,j,i))
			+8.*(get_bv(l+1,k,j,i)-get_bv(l-1,k,j,i))
			)*dzi12;
		return w;
	}
	double get_f_xx(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l  ,k  ,j+2,i) + get_bv(l  ,k  ,j-2,i))
			+16.*( get_bv(l  ,k  ,j+1,i) + get_bv(l  ,k  ,j-1,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dxi*dxi12;
		return w;
	}
	double get_f_yy(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l  ,k+2,j  ,i) + get_bv(l  ,k-2,j  ,i))
			+16.*( get_bv(l  ,k+1,j  ,i) + get_bv(l  ,k-1,j  ,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dyi*dyi12;
		return w;
	}
	double get_f_zz(int l,int k,int j,int i) const{
		double w=(    -(     get_bv(l+2,k  ,j  ,i) + get_bv(l-2,k  ,j  ,i))
			+16.*( get_bv(l+1,k  ,j  ,i) + get_bv(l-1,k  ,j  ,i))
			-30.*  get_bv(l  ,k  ,j  ,i)
			)*dzi*dzi12;
		return w;
	}

	double get_f_xy(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l  ,k+2,j+2,i) - get_bv(l  ,k+2,j-2,i))
						+8.*(get_bv(l  ,k+2,j+1,i) - get_bv(l  ,k+2,j-1,i)))
						-(   -(get_bv(l  ,k-2,j+2,i) - get_bv(l  ,k-2,j-2,i))
						+8.*(get_bv(l  ,k-2,j+1,i) - get_bv(l  ,k-2,j-1,i))))
						+8.*(  (-(get_bv(l  ,k+1,j+2,i) - get_bv(l  ,k+1,j-2,i))
						+8.*(get_bv(l  ,k+1,j+1,i) - get_bv(l  ,k+1,j-1,i)))
						-(-(get_bv(l  ,k-1,j+2,i) - get_bv(l  ,k-1,j-2,i))
						+8.*(get_bv(l  ,k-1,j+1,i) - get_bv(l  ,k-1,j-1,i))))
						)*dxi12*dyi12;
		return w;
	}
	double get_f_xz(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l+2,k  ,j+2,i) - get_bv(l+2,k  ,j-2,i))
						+8.*(get_bv(l+2,k  ,j+1,i) - get_bv(l+2,k  ,j-1,i)))
						-(   -(get_bv(l-2,k  ,j+2,i) - get_bv(l-2,k  ,j-2,i))
						+8.*(get_bv(l-2,k  ,j+1,i) - get_bv(l-2,k  ,j-1,i))))
						+8.*(  (-(get_bv(l+1,k  ,j+2,i) - get_bv(l+1,k  ,j-2,i))
						+8.*(get_bv(l+1,k  ,j+1,i) - get_bv(l+1,k  ,j-1,i)))
						-(-(get_bv(l-1,k  ,j+2,i) - get_bv(l-1,k  ,j-2,i))
						+8.*(get_bv(l-1,k  ,j+1,i) - get_bv(l-1,k  ,j-1,i))))
						)*dxi12*dzi12;
		return w;
	}
	double get_f_yz(int l,int k,int j,int i) const{
		double w=( -( ( -(   get_bv(l+2,k+2,j  ,i) - get_bv(l+2,k-2,j  ,i))
						+8.*(get_bv(l+2,k+1,j  ,i) - get_bv(l+2,k-1,j  ,i)))
						-(   -(get_bv(l-2,k+2,j  ,i) - get_bv(l-2,k-2,j  ,i))
						+8.*(get_bv(l-2,k+1,j  ,i) - get_bv(l-2,k-1,j  ,i))))
						+8.*(  (-(get_bv(l+1,k+2,j  ,i) - get_bv(l+1,k-2,j  ,i))
						+8.*(get_bv(l+1,k+1,j  ,i) - get_bv(l+1,k-1,j  ,i)))
						-(-(get_bv(l-1,k+2,j  ,i) - get_bv(l-1,k-2,j  ,i))
						+8.*(get_bv(l-1,k+1,j  ,i) - get_bv(l-1,k-1,j  ,i))))
						)*dyi12*dzi12;
		return w;
	}
	double get_ipol_x_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l,k,j-2,i)+9.*get_bv(l,k,j-1,i)+9.*get_bv(l,k,j,i)-get_bv(l,k,j+1,i));
		return w;
	}
	double get_ipol_y_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l,k-2,j,i)+9.*get_bv(l,k-1,j,i)+9.*get_bv(l,k,j,i)-get_bv(l,k+1,j,i));
		return w;
	}
	double get_ipol_z_lower_mid(int l,int k,int j,int i) const{
		double w=0.0625*(-get_bv(l-2,k,j,i)+9.*get_bv(l-1,k,j,i)+9.*get_bv(l,k,j,i)-get_bv(l+1,k,j,i));
		return w;
	}
		
		
	////////////////////////////////////////
	//  SET func.
	////////////////////////////////////////
	void set_fluidevo(bool f){
		fluidevo=f;
		return;
	}
	void set_scalarevo(bool s){
		scalarevo=s;
		return;
	}
	void set_mrf(bool s){
		mrf=s;
		return;
	}
	void set_t(double time){
		t=time;
		return;
	}
	void set_dt(double time){
		dt=time;
		return;
	}
	void set_dt0(double time){
		dt0=time;
		return;
	}
	void set_dtp(double time){
		dtp=time;
		return;
	}
	void set_dtpp(double time){
		dtpp=time;
		return;
	}
	void set_tmax(double time){
		tmax=time;
		return;
	}
	void set_cfl(double c){
		cfl=c;
		return;
	}
	void set_etaa(double e){
		etaa=e;
		return;
	}
	void set_etab(double e){
		etab=e;
		return;
	}
	void set_etabb(double e){
		etabb=e;
		return;
	}
	void set_lambda(double l){
		lambda=l;
		return;
	}
	void set_Hb(double hb){
		Hb=hb;
		return;
	}
	void set_tini(double t){
		tini=t;
		return;
	}
	void set_KOep(double l){
		KOep=l;
		return;
	}
	void set_exg(int eg){
		exg=eg;
		return;
	}
	void set_hl(int h){
		hl=h;
		return;
	}
	void set_amp(double a){
		amp=a;
		return;
	}
	void set_fluidw(double fw){
		fluidw=fw;
		return;
	}
		void set_scalarm(double sm){
		scalarm=sm;
		return;
	}
	void set_Mkap(double k){
		kap_MUSCL=k;
		return;
	}
	void set_b(double b){
		b_minmod=b;
		return;
	}
	void set_exc(bool e){
		exc=e;
		return;
	}
	void set_llmin(int lll){
		llmin=lll;
		return;
	}
	int& set_bflag(int l,int k,int j){
		return bflag[l-lmin][k-kmin][j-jmin];
	}
	int& set_hflag(int l,int k,int j){
		return hflag[l-lmin][k-kmin][j-jmin];
	}
	double& set_bv(int l,int k,int j,int i){
		return bv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_dbv(int l,int k,int j,int i){
		return dbv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bv0(int l,int k,int j,int i){
		return bv0[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bv1(int l,int k,int j,int i){
		return bv1[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_bvr(int l,int k,int j,int i){
		return bvr[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_tmpipol(int l,int k,int j,int i){
		return tmpipol[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_refz(int l,int k,int j){
		return refz[l-lmin][k-kmin][j-jmin];
	}
	double& set_con(int l,int k,int j,int i){
		return con[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_coordZ(int l){
		return coordZ[l-lmin];
	}
	double& set_flat_df2z(int l){
		return flat_df2z[l-lmin];
	}
	double& set_flat_Gamz(int l){
		return flat_Gamz[l-lmin];
	}
	double& set_flat_dGamz(int l){
		return flat_dGamz[l-lmin];
	}
	double& set_alp_e(int l){
		return alp_e[l-lmin];
	}
	double& set_bet_e(int l){
		return bet_e[l-lmin];
	}
	double& set_gam_e(int l){
		return gam_e[l-lmin];
	}
	double& set_alp_o(int l){
		return alp_o[l-lmin];
	}
	double& set_bet_o(int l){
		return bet_o[l-lmin];
	}
	double& set_gam_o(int l){
		return gam_o[l-lmin];
	}
	double& set_ddy(int l,int i){
		return ddy[i][l-lmin];
	}
	double& set_outv(int l,int k,int j,int i){
		return outv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_primv(int l,int k,int j,int i){
		return primv[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_x(int l,int k,int j,int i){
		return flux_x[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_y(int l,int k,int j,int i){
		return flux_y[i][l-lmin][k-kmin][j-jmin];
	}
	double& set_flux_z(int l,int k,int j,int i){
		return flux_z[i][l-lmin][k-kmin][j-jmin];
	}
	double*** set_bv(int i) const{
		return bv[i];
	}
	double*** set_dbv(int i) const{
		return dbv[i];
	}
	double*** set_bv0(int i) const{
		return bv0[i];
	}
	double*** set_bv1(int i) const{
		return bv1[i];
	}
	double*** set_bvr(int i) const{
		return bvr[i];
	}
	double*** set_tmpipol(int i) const{
		return tmpipol[i];
	}

	////////////////////////////////////////
	//  SET zero or unity func.
	////////////////////////////////////////
	void set_zero_all()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv[i][l][k][j] = 0.;
						dbv[i][l][k][j] = 0.;
						bv0[i][l][k][j] = 0.;
						bvr[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<ntmpipol;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						tmpipol[i][l][k][j] = 0.;
					}
				}
			}
		}

		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nc;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						con[i][l][k][j] = 0.;
					}
				}
			}
		}

		#pragma omp parallel for
		for(int l=0;l<nz;l++){
			flat_df2z[l] = 0.;
			coordZ[l] = 0.;
			flat_Gamz[l] = 0.;
			flat_dGamz[l] = 0.;
		}

		#pragma omp parallel for
		for(int l=0;l<nz;l++){
			alp_e[l] = 0.;
			bet_e[l] = 0.;
			gam_e[l] = 0.;
			alp_o[l] = 0.;
			bet_o[l] = 0.;
			gam_o[l] = 0.;
			
			for(int i=0;i<nn;i++)
			 ddy[i][l] = 0.;
		}

		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nv;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						outv[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int k=0;k<ny;k++){
				for(int j=0;j<nx;j++){
					refz[l][k][j] = 0.;
				}
			}
		}

		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<npr;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						primv[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<npr-2;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						flux_x[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<npr-2;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						flux_y[i][l][k][j] = 0.;
					}
				}
			}
		}
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<npr-2;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						flux_z[i][l][k][j] = 0.;
					}
				}
			}
		}
		

		return;
	}
	
	void set_zero_primv()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<npr;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						primv[i][l][k][j] = 0.;
					}
				}
			}
		}
	}
	
	void set_zero()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_zero_0()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv0[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_zero_1()
	{
		#pragma omp parallel for 
		for(int i=0;i<nn;i++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv1[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}

	void set_zero_d()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				if(!nonzero[i])
				continue;
				
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						dbv[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}

	void set_zero_d_exc()
	{
		#pragma omp parallel for 
		for(int l=lli;l<=lui;l++)
		{
			for(int i=0;i<nn;i++)
			{
				if(!nonzero[i])
				continue;
				
				for(int k=kli;k<=kui;k++)
				{
					for(int j=jli;j<=jui;j++)
					{
						set_dbv(l,k,j,i) = 0.;
					}
				}
			}
		}
		return;
	}

	void set_zero_r()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bvr[i][l][k][j] = 0.;
					}
				}
			}
		}
		return;
	}
	
	void set_bflag_zero(){
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int k=0;k<ny;k++){
				for(int j=0;j<nx;j++){
					bflag[l][k][j] = 0;
				}
			}
		}
		return;
	}

	void set_hflag_zero()
	{	
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int k=0;k<ny;k++){
				for(int j=0;j<nx;j++){
					hflag[l][k][j] = 0;
				}
			}
		}
		return;
	}

	////////////////////////////////////////
	//  UPDATE func.
	////////////////////////////////////////

	void setv0()
	{
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv0[i][l][k][j] = bv[i][l][k][j];
					}
				}
			}
		}

		return;
	}
	
	void set01()
	{
		#pragma omp parallel for
		for(int i=0;i<nn;i++){
			for(int l=0;l<nz;l++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv1[i][l][k][j] = bv0[i][l][k][j];
					}
				}
			}
		}
	
		return;
	}
	
	void set0v(){
		
		#pragma omp parallel for 
		for(int l=0;l<nz;l++){
			for(int i=0;i<nn;i++){
				for(int k=0;k<ny;k++){
					for(int j=0;j<nx;j++){
						bv[i][l][k][j] = bv0[i][l][k][j];
					}
				}
			}
		}
		return;
	}

	void runge_kutta(double dt)
	{
		#pragma omp parallel for 
		for(int l=lmin;l<=lmax;l++){
			for(int i=0;i<nn;i++){
				if(!nonzero[i])
				continue;

				for(int k=kmin;k<=kmax;k++){
					for(int j=jmin;j<=jmax;j++){
						set_bvr(l,k,j,i)=get_bvr(l,k,j,i)+get_dbv(l,k,j,i)*dt;
					}
				}
			}
		}
		return;
	}
	
	void new_bv(double dt)
	{
		#pragma omp parallel for 
		for(int l=lmin;l<=lmax;l++){
			for(int i=0;i<nn;i++){
				if(!nonzero[i])
				continue;

				for(int k=kmin;k<=kmax;k++){
					for(int j=jmin;j<=jmax;j++){
						set_bv(l,k,j,i)=get_bv0(l,k,j,i)+get_dbv(l,k,j,i)*dt;
					}
				}
			}
		}
		return;
	}

	void new_bv4()
	{
		#pragma omp parallel for 
		for(int l=lmin;l<=lmax;l++){
			for(int i=0;i<nn;i++){
				if(!nonzero[i])
				continue;

				for(int k=kmin;k<=kmax;k++){
					for(int j=jmin;j<=jmax;j++){
						set_bv(l,k,j,i)=get_bv0(l,k,j,i)+get_bvr(l,k,j,i);
					}
				}
			}
		}
		return;
	}

	////////////////////////////////////////
	//  check functions
	////////////////////////////////////////

	void check_const()
	//checking constraint start
	{
		ham=0.;
		double hamtmp=0.;
		double momtmp=0.;
		hammax=0.;
		mommax=0.;
		mom=0.;
		lhm=-1;
		lmm=-1;
		
		dGam=0.;
		double dGamtmp=0.;
		dGammax=0.;
		
		int excnum=0;
			
		#pragma omp parallel 
		{
			double hammaxp=0.;
			double mommaxp=0.;
			double lhmp=-1;
			double lmmp=-1;
			double dGammaxp=0.;
			
			// @@@@@@@@@@@@   constraints   @@@@@@@@@@@
			//  0:normalized ham, 1:ham, 2:normalized momz, 3:momz
			//  4:Gamma^z-D_gamma^z

			#pragma omp for reduction(+:hamtmp,momtmp,dGamtmp,excnum) 
			for(int l=llmin;l<=lui-2*negb;l++)
			{
				int k=0;
				int j=0;

				if(get_hflag(l,k,j)!=0)
				{
					excnum++;
					continue;
				}
				
				hamtmp+= abs(get_con(l,k,j,0));
				momtmp+= abs(get_con(l,k,j,2));
				dGamtmp+= abs(get_con(l,k,j,4));
				
				if(hammaxp<abs(get_con(l,k,j,0)))
				{
					hammaxp=abs(get_con(l,k,j,0));
					lhmp=l;
				}
						
				double momcompmax=abs(get_con(l,k,j,2));
				
				if(mommaxp<momcompmax)
				{
					mommaxp=momcompmax;
					lmmp=l;
				}

				double dGamcompmax=abs(get_con(l,k,j,4));
				
				if(dGammaxp<dGamcompmax)
				{
					dGammaxp=dGamcompmax;
				}
			}
			
			if(hammaxp > hammax)
			{
				#pragma omp critical (hamcheck)
				{
					if(hammaxp > hammax) 
					{
						hammax=hammaxp;
						lhm=lhmp;
					}
				}
			}

			if(mommaxp > mommax)
			{
				#pragma omp critical (momcheck)
				{
					if(mommaxp > mommax) 
					{
						mommax=mommaxp;
						lmm=lmmp;
					}
				}
			}
			
			if(dGammaxp > dGammax)
			{
				#pragma omp critical (dGamcheck)
				{
					if(dGammaxp > dGammax) 
					{
						dGammax=dGammaxp;
					}
				}
			}

		}

		#pragma omp barrier

		double fac=lui-2*negb-llmin+1-excnum;
		ham=hamtmp/fac;
		mom=momtmp/fac;
		dGam=dGamtmp/fac;
		
	}
	
	void check_Kremax()
	//checking constraint start
	{
		Kremax=0.;
		lkm=-1;
		
		#pragma omp parallel 
		{
			double Kremaxp=0.;
			double lkmp=-1;
			
			#pragma omp for  
			for(int l=lli;l<=lui;l++)
			{
				int k=0;
				int j=0;

				if(get_hflag(l,k,j)!=0)
				continue;
				
				if(Kremaxp<abs(get_outv(l,k,j,0)))
				{
					Kremaxp=abs(get_outv(l,k,j,0));
					lkmp=l;
				}
			}
			
			if(Kremaxp > Kremax)
			{
				#pragma omp critical (Krecheck)
				{
					if(Kremaxp > Kremax) 
					{
						Kremax=Kremaxp;
						lkm=lkmp;
					}
				}
			}
		}
	}

	void check_Weylmax()
	//checking constraint start
	{
		Weylmax=0.;
		lwm=-1;
		
		#pragma omp parallel 
		{
			double Weylmaxp=0.;
			double lwmp=-1;
			
			#pragma omp for  
			for(int l=lli;l<=lui;l++)
			{
				int k=0;
				int j=0;
				if(get_hflag(l,k,j)!=0)
				continue;
				
				if(Weylmaxp<abs(get_outv(l,k,j,3)))
				{
					Weylmaxp=abs(get_outv(l,k,j,3));
					lwmp=l;
				}
			}
			
			if(Weylmaxp > Weylmax)
			{
				#pragma omp critical (Weylcheck)
				{
					if(Weylmaxp > Weylmax) 
					{
						Weylmax=Weylmaxp;
						lwm=lwmp;
					}
				}
			}
		}
	}
	
	void check_horizon(ofstream& fout)
	{
		double Comp=0.,Arad=0.;
		double temphr[hnmax][2];
		int hnc=0;

		for(int h=0;h<hnmax;h++)
		{
			temphr[h][0]=10.;
			temphr[h][1]=10.;
		}

		for(int l=llmin;l<=lui;l++)
		{
			int k=0;
			int j=0;
			
			double preArad,prerad,preComp,rad;
			double Mass;
			
			Mass=get_outv(l,k,j,9);
			preArad=Arad;
			prerad=get_z(l-1);
			Arad=get_outv(l,k,j,5);
			rad=get_z(l);
			preComp=Comp;
			
			if(l==lli)
			Comp=0.;
			else
			Comp=2.*Mass/Arad;

			if(get_hflag(l-1,k,j)==1 && exc==true)
			continue;

			if((1.-Comp)*(1.-preComp)<=0.)
			{
				if(hnc+1>=hnmax)
				{
					cout << "too many horizons" << endl;
					//exit(1);
				}
				else
				{
					temphr[hnc][0]=(prerad+(rad-prerad)*(preComp-1.)/(preComp-Comp));
					temphr[hnc][1]=(preArad+(Arad-preArad)*(preComp-1.)/(preComp-Comp));
					hnc++;
				}
			}
			
			if(hform==false)
			{
				if(get_outv(l,k,j,4)<0. && get_outv(l,k,j,6)<0.)
				{
					hform=true;

					if(exc==false && exg!=0) 
					{
						exc=true;
						set_excflags();

						for(int ll=lli;ll<l;ll++)
						{
							set_hflag(ll,k,j)=1;
						}
					}
				}
			}
		}

		double devs[hnmax];
		int hcpn[hnmax];
		int noashc=0;
		for(int h=0;h<hnmax;h++)
		{
			devs[h]=999999999.;
			hcpn[h]=hnmax;
		}
		
		for(int hc=0;hc<hnc;hc++)
		{
			int hpsqn=hnmax;	
			int hpsqn2=hnmax;	
			double dev=999999999.;
			double dev2=999999999.;
			
			for(int hp=0;hp<hn;hp++)
			{
				if(horis[hp][0]>1.)
				continue;
				
				double devtemp=abs(temphr[hc][0]-horis[hp][0]);
				if(devtemp<dev)
				{
					dev2=dev;
					hpsqn2=hpsqn;
					dev=devtemp;
					hpsqn=hp;
				}
			}

			if(dev<devs[hpsqn])
			{
				if(hcpn[hpsqn]!=hnmax)
				{
					int hcc=hcpn[hpsqn];
					double devv=999999999.;
					int hpsqnn=hnmax;
					for(int hpp=0;hpp<hn;hpp++)
					{
						if(hpp==hpsqn || horis[hpp][0]>1.)
						continue;

						double devtemp=abs(temphr[hcc][0]-horis[hpp][0]);
						if(devtemp<devv)
						{
							devv=devtemp;
							hpsqnn=hpp;
						}
					}

					if(devv<devs[hpsqnn])
					{
						if(hcpn[hpsqnn]!=hnmax)
						{
							hcpn[hn+noashc]=hcpn[hpsqnn];
							devs[hn+noashc]=0.;
							noashc++;

							if(hn+noashc>=hnmax)
							{
								cout << "too many horizons 1" << endl;
								noashc--;
							}
						}

						devs[hpsqnn]=devv;
						hcpn[hpsqnn]=hcc;
					}
					else
					{
						hcpn[hn+noashc]=hcc;
						devs[hn+noashc]=0.;
						noashc++;

						if(hn+noashc>=hnmax)
						{
							cout << "too many horizons 2" << endl;
							noashc--;
						}
					}
				}
				
				devs[hpsqn]=dev;
				hcpn[hpsqn]=hc;
			}
			else if(dev2<devs[hpsqn2])
			{
				if(hcpn[hpsqn2]!=hnmax)
				{
					hcpn[hn+noashc]=hcpn[hpsqn2];
					devs[hn+noashc]=0.;
					noashc++;

					if(hn+noashc>=hnmax)
					{
						cout << "too many horizons 3" << endl;
						noashc--;
					}
				}
				devs[hpsqn2]=dev2;
				hcpn[hpsqn2]=hc;
			}
			else
			{
				hcpn[hn+noashc]=hc;
				devs[hn+noashc]=0.;
				noashc++;

				if(hn+noashc>=hnmax)
				{
					cout << "too many horizons 4" << endl;
					noashc--;
				}
			}
		}

		if(noashc>0)
		{
			hn=hn+noashc;
			cout << "hn=" << hn << endl;
		}

		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(20) << get_t();							//1

		for(int h=0;h<hnmax;h++)
		{
			if(hcpn[h]<hnmax && temphr[hcpn[h]][0]<1.)
			{
				//cout << "horis substituted" << endl;
				horis[h][0]=temphr[hcpn[h]][0];
				horis[h][1]=temphr[hcpn[h]][1];
			}
			else
			{
				horis[h][0]=10.;
				horis[h][1]=10.;
			}

			if(horis[h][0]>1.)
			{
				fout 
				<< " "  << setw(20) << "-"
				<< " "  << setw(20) << "-";	//2
			}
			else
			{
				fout 
				<< " "  << setw(20) << horis[h][0] 
				<< " "  << setw(20) << horis[h][1];
			}
		}
		fout << endl;
		
		return;
	}


	void check_neck(ofstream& fout)
	{
		double temphr[hnmax];
		double drR=0.;
		int hnc=0;

		for(int h=0;h<hnmax;h++)
		{
			temphr[h]=10.;
		}

		for(int l=llmin;l<=lui;l++)
		{
			int k=0;
			int j=0;
			
			double prerad,predrR,rad;
			
			predrR=drR;
			prerad=get_z(l-1);
			rad=get_z(l);
			
			if(l==lli)
			{
				drR=1.0e-10;
				continue;
			}
			else
			drR=get_outv(l,k,j,11);

			if(get_hflag(l-1,k,j)==1 && exc==true)
			continue;

			if(drR*predrR<=0.)
			{
				if(hnc+1>=hnmax)
				{
					cout << "too many necks" << endl;
					//exit(1);
				}
				else
				{
					temphr[hnc]=(prerad+(rad-prerad)*predrR/(predrR-drR));
					hnc++;
				}
			}
		}

		double devs[hnmax];
		int hcpn[hnmax];
		int noashc=0;
		for(int h=0;h<hnmax;h++)
		{
			devs[h]=999999999.;
			hcpn[h]=hnmax;
		}
		
		for(int hc=0;hc<hnc;hc++)
		{
			int hpsqn=hnmax;	
			int hpsqn2=hnmax;	
			double dev=999999999.;
			double dev2=999999999.;
			
			for(int hp=0;hp<nneck;hp++)
			{
				if(neck[hp]>1.)
				continue;
				
				double devtemp=abs(temphr[hc]-neck[hp]);
				if(devtemp<dev)
				{
					dev2=dev;
					hpsqn2=hpsqn;
					dev=devtemp;
					hpsqn=hp;
				}
			}

			if(dev<devs[hpsqn])
			{
				if(hcpn[hpsqn]!=hnmax)
				{
					int hcc=hcpn[hpsqn];
					double devv=999999999.;
					int hpsqnn=hnmax;
					for(int hpp=0;hpp<nneck;hpp++)
					{
						if(hpp==hpsqn || neck[hpp]>1.)
						continue;

						double devtemp=abs(temphr[hcc]-neck[hpp]);
						if(devtemp<devv)
						{
							devv=devtemp;
							hpsqnn=hpp;
						}
					}

					if(devv<devs[hpsqnn])
					{
						if(hcpn[hpsqnn]!=hnmax)
						{
							hcpn[nneck+noashc]=hcpn[hpsqnn];
							devs[nneck+noashc]=0.;
							noashc++;

							if(nneck+noashc>=hnmax)
							{
								cout << "too many necks 1" << endl;
								noashc--;
							}
						}

						devs[hpsqnn]=devv;
						hcpn[hpsqnn]=hcc;
					}
					else
					{
						hcpn[nneck+noashc]=hcc;
						devs[nneck+noashc]=0.;
						noashc++;

						if(nneck+noashc>=hnmax)
						{
							cout << "too many neck 2" << endl;
							noashc--;
						}
					}
				}
				
				devs[hpsqn]=dev;
				hcpn[hpsqn]=hc;
			}
			else if(dev2<devs[hpsqn2])
			{
				if(hcpn[hpsqn2]!=hnmax)
				{
					hcpn[nneck+noashc]=hcpn[hpsqn2];
					devs[nneck+noashc]=0.;
					noashc++;

					if(nneck+noashc>=hnmax)
					{
						cout << "too many necks 3" << endl;
						noashc--;
					}
				}
				devs[hpsqn2]=dev2;
				hcpn[hpsqn2]=hc;
			}
			else
			{
				hcpn[nneck+noashc]=hc;
				devs[nneck+noashc]=0.;
				noashc++;

				if(nneck+noashc>=hnmax)
				{
					cout << "too many necks 4" << endl;
					noashc--;
				}
			}
		}

		if(noashc>0)
		{
			nneck=nneck+noashc;
			cout << "nneck=" << nneck << endl;
		}

		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(20) << get_t();							//1

		for(int h=0;h<hnmax;h++)
		{
			if(hcpn[h]<hnmax && temphr[hcpn[h]]<1.)
			{
				//cout << "horis substituted" << endl;
				neck[h]=temphr[hcpn[h]];
			}
			else
			{
				neck[h]=10.;
			}

			if(neck[h]>1.)
			{
				fout 
				<< " "  << setw(20) << "-";	//2
			}
			else
			{
				fout 
				<< " "  << setw(20) << neck[h];
			}
		}
		fout << endl;
		
		return;
	}

	////////////////////////////////////////
	//  flag setting func.
	////////////////////////////////////////
	
	void set_excflags()
	{
		set_bflag_zero();
		int j=0;
		int k=0;

		set_bflag(lli+exg-1,k,j)=1;
		set_bflag(lli+exg-2,k,j)=2;
		set_bflag(lli+exg-3,k,j)=3;
		set_bflag(lli+exg-4,k,j)=4;
		set_bflag(lli+exg-5,k,j)=5;
		set_bflag(lli+exg-6,k,j)=6;

		for(int l=lmin;l<lli+exg-6;l++)
		{
			set_bflag(l,k,j)=-1;
		}
	}
	
	////////////////////////////////////////
	//  fluid func.
	////////////////////////////////////////
	
	double pres(double rho);
	double dpres(double rho);
	void get_rhoGam(double Ene, double S,double& rho,double& Gam);
	double minmod(double a,double b);
	double sign(double A);
	void dyntoprim();
	
	////////////////////////////////////////
	//  Interpolation func.
	////////////////////////////////////////
	
	double ipol( double rr,double *xx,double *yy,int order );
	double bv_ipol(int lc,double zc,int order,int number);
	double bv0_ipol(int lc,double zc,int order,int number);
	double bv1_ipol(int lc,double zc,int order,int number);
	
	void put_ddy(int i);									// second derivatives for spline interpolation
	double bv_spline_ipol(int lc,double zc,int number);		// spline interpolation
	void set_luvec();										// LU decomposition for spline interpolation
	
	////////////////////////////////////////
	//  BSSN func.
	////////////////////////////////////////
	
	void BSSN_adv();
	void BSSN(int itype);
	void enforce_const();
	void KOdiss(int lup);
	void flux_fill();
	void excision();
	void excision_quadfunc();
	
	////////////////////////////////////////
	//  Initial func.
	////////////////////////////////////////

	void initial_continue(ifstream& fcontinue);
	double funcf(double X);
	double ifuncf(double X);
	double df(double X);
	double ddf(double X);
	double dddf(double X);
	void set_flat();
	void set_refz();
	void set_Gam();
	void set_enemomini();
	void initial_params(double cfli,double etaai,double etabi,double etabbi,double lambdai,double dt0i,double dtpi,double dtppi,double ti,double tinii,double Hbi,double KOepi,int exgi,double fluidwi,double scalarmi,double kap_MUSCLi,double b_minmodi);
	
	////////////////////////////////////////
	//  OUTPUT func.
	////////////////////////////////////////

	void print_z(ofstream& fn, int j,int k);
	void print_bz(ofstream& fn, int j,int k);
	void print_Kremax(ofstream& fout);
	void print_const(ofstream& fout);
	void print_all(ofstream& fout);
	void print_ipolz(ofstream& fout);
	

	////////////////////////////////////////
	//  Potential func.
	////////////////////////////////////////
	
	double funcV(double p)
	{
		double w=0.5*pow(scalarm*p,2);
		
		return(w);
	}
	
	double funcdV(double p)
	{
		double w=pow(scalarm,2)*p;
		
		return(w);
	}

};

class Fmv : public Fmv0{
// private: 
// 	int laymax;
// 	bool *mrflag;
// 	int *lbs;
// 	double *alp_fmr;

public:
	Fmv(int tabs,int tabsx,int lupper,
	double xupper,double zupper,double am,bool fld, bool scl, bool cuev) : Fmv0(tabs,tabsx,lupper,
	xupper,zupper,am,fld, scl, cuev){
		layn=0;
		// laymax=laynum;
		
		// mrflag = new bool[laymax+1];
		// lbs = new int[laymax+1];
		// alp_fmr = new double[laymax+1];
		
		// for(int n=0;n<=laymax;n++)
		// {
		// 	mrflag[n]=false;
		// 	lbs[n]=0;
		// }
	}

 public:
// 	bool get_mrflag(int n) const{
// 		return mrflag[n];
// 	}
// 	int get_lbs(int n) const{
// 		return lbs[n];
// 	}
// 	double get_alp_fmr(int n) const{
// 		return alp_fmr[n];
// 	}
	
// 	bool& set_mrflag(int n){
// 		return mrflag[n];
// 	}
// 	int& set_lbs(int n){
// 		return lbs[n];
// 	}
// 	double& set_alp_fmr(int n){
// 		return alp_fmr[n];
// 	}

	////////////////////////////////////////
	//  BOUNDARY func.
	////////////////////////////////////////
	
	void asymcond(int l,int i,double bgv0,double bgv,double dt,int itype);
	void asymcond(int l,int i,double bgv);
	void boundary_asym(int itype);

	////////////////////////////////////////
	//  func. for specific initial conditions
	////////////////////////////////////////
	void set_initial_final();
	void base_initial_continue(ifstream& fcontinue);

	void initial_iso_longwave(double mu,double k,double inr,double L);
	void initial_longwave(double mu,double k,double inr,double L);
	double Phi(double r,double mu,double kk, double inr, double L);
	double dzPhi(double r,double mu,double kk, double inr, double L);
	double ddzPhi(double r,double mu,double kk, double inr, double L);

	double Psi(double r,double mu,double kk, double inr, double L);
	double dzPsi(double r,double mu,double kk, double inr, double L);
	double ddzPsi(double r,double mu,double kk, double inr, double L);
	double dddzPsi(double r,double mu,double kk, double inr, double L);

};


class Fmv1 : public Fmv0{
private: 
	int lb;
	Fmv0* llay;
	Fmv1* ulay;

public:
	Fmv1(int tabs,int tabsx, int lupper,double xupper,double zupper,double am,bool fld, bool scl, bool cuev, Fmv0* lolay)
	 : Fmv0(tabs,tabsx, lupper,xupper,zupper,am,fld, scl, cuev)
	{
		layn=lolay->get_layn()+1;
		llay=lolay;
		
		cout << "Fmv1 for layer #" << layn << " constructer done" << endl;
	}

public:
	int get_lb() const{
		return lb;
	}
	
	void set_lb(int p){
		lb=p;
		return;
	}
	void set_ulay(Fmv1* uplay){
		ulay=uplay;
		return;
	}
	
	void boundary_fmr();
	void set_boundary(int btype);
	void set_fmr_initial();
	void evolve();
	void refine_llay();
	void onestep(int btype);
	void fmr_initial_continue(ifstream& fcontinue);

};

#endif
