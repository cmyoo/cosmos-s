/**
 * @file cosmos_s.cpp
 * @brief Main driver program for the COSMOS-S simulation code.
 *
 * This file contains the main entry point (`main`) for the simulation,
 * parameter reading functions (`initial_read`, `initial_fmr`),
 * continuation file handling (`check_continue_file`, `output_params`),
 * and orchestrates the overall simulation flow including initialization,
 * time stepping, FMR layer management, output, and finalization.
 */
//**********************************************************************************************//
//----------------------------------------------------------------------------------------------//
//                     @@@@@   @@@@    @@@@    @@@ @@    @@@@    @@@@                           //
//                   @@@      @@   @  @@   @  @@@ @@ @  @@   @  @@   @                          //
//                   @@       @@   @   @@@    @@  @  @  @@   @   @@@                            //
//                   @@       @@   @ @    @@  @@     @  @@   @ @    @@                          //
//                     @@@@@   @@@@   @@@@@   @@     @   @@@@   @@@@@                           //
//----------------------------------------------------------------------------------------------//
//**********************************************************************************************//
//----------------------------------------------------------------------------------------------//
//                 cosmos-s                                                                     //
//                       ver. 1.00 coded by Chulmoon Yoo                                        //
//----------------------------------------------------------------------------------------------//
// Main file   :: cosmos_s.cpp                                                                  //
// header      :: cosmos_s.h                                                                    //
// definition  :: cosmos_s_bssn.cpp cosmos_s_initial.cpp cosmos_s_output.cpp                    //
//                cosmos_s_boundary.cpp cosmos_s_ipol.cpp cosmos_s_fluid.cpp                    //
//                cosmos_s_fmr.cpp                                                              //
// how to make :: makefile                                                                      //
//----------------------------------------------------------------------------------------------//
//    Compile :: $ make                                                                         //
//               Choose compiler by changing makefile.                                          //
//    OpenMP  :: $ export OMP_NUM_THREADS=#                                                     //
//               Set the number of CPU cores you use.                                           //
//    Execute :: $ taskset -c 0-# ./cosmos_s                                                    //
//               Set CPU cores directly on the machine.                                         //
//    Clean   :: $ make clean                                                                   //
//               Delete .o files.                                                               //
//**********************************************************************************************//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cstring>
#include "cosmos_s.h"						/**< Header for the COSMOS-S BSSN evolution classes. */

using namespace std;


void check_continue_file(					//checking parameter consistency with the continue file
ifstream& fcontinue, 								//initial parameter file
bool& fld,									//for fluid evolution
bool& scl,									//for scalar evolution
int& exc,									//excision
int& exg,									//excision grid
int& tab, 									//tab num for buf grids
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *lbs									//grid number for fmr region on z-axis
);

void output_params(							//output all needed parameters for continue
ofstream& fcontinue, 						//initial parameter file
bool& fld,									//for fluid evolution
bool& scl,									//for scalar evolution
int& exc,									//excision
int& exg,									//excision grid
int& tab, 									//tab num for buf grids
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *lbs									//grid number for fmr region on z-axis
);


void initial_read(							//reading inital para
ifstream& fin, 								//initial parameter file
long int& mstep, 							//max step num for time evo
double& tmax, 								//max time to evolve
int& tab, 									//tab num for buf grids
double& amp, 								//amplitude for inhomogeneous grid
int& nzmax, 								//max grid num of z
double& xmax, 								//max coord val of x
double& zmax, 								//max coord val of z
double& cfl, 								//CFL para
double& cdt, 								//factor for cosmological time scale(dt=cfl*cdt*1/H)
double& etaa, 								//etaa for gauge
double& etab, 								//etab for gauge
double& etabb, 								//etabb for gauge
double& KOep, 								//factor for Kreiss-Oliger disspation term
int& exg,									//excision grid num
int& contn,									//1 for continue
char *file_continue,						//continue file
double& inr, 								//inner radius of matching region for the profile
double& outr, 								//outer radius of matching region for the profile
double& mu, 								//amplitude of the perturbation
double& kk, 								//scale of the perturbation
double& mus, 								//amplitude of the perturbation of scalar field
double& kks, 								//scale of the perturbation of scalar field
double& Hb, 								//initial Hubble
double& fluidw, 							//fluidw
double& Mkap, 								//kappa in MUSCL
double& bminmod, 							//b in minmod function
double& ptintval1,							//1st part print interval boundary time
double& ptintval2,							//2nd
double& changept,							//printing interval change time
int& horicheckintv,							//horizon formation check interval
int& constoutintv
);


void initial_fmr(							//reading inital para for fmr
ifstream& fin, 								//initial parameter file for fmr
int& laymax, 								//max fmr layer number
int *lbs, 									//grid number for fmr region
double *alp_fmr								//values of the lapse for starting fmr
);


/**
 * @brief Main entry point for the COSMOS-S simulation.
 *
 * This function orchestrates the entire simulation process. It reads initial parameters,
 * sets up the simulation grid(s) potentially with Fixed Mesh Refinement (FMR),
 * handles continuation from previous runs, executes the main time evolution loop,
 * manages output, and performs final cleanup.
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line argument strings (currently unused).
 * @return 0 on successful completion, non-zero on error.
 */
int main(int argc,char* argv[])
{
	bool	fld,							//for fluid evolution
		scl,								//for scalar evolution
		cuev;								//for curvature evaluation

	long int mstep; 						//max step num for time evo
	int 	tab, 							//tab num for buf grids
		nzmax; 								//max grid num of z

	int ln=0;

	double amp,								//for inhom grid parameter
		xmax, 								//max coord val of x
		zmax, 								//max coord val of z
		t, 									//time
		tmax, 								//max time to evolve
		cfl, 								//CFL para
		cdt, 								//factor for cosmological time scale(dt=cfl*cdt*1/H)
		etaa, 								//etaa for gauge
		etab, 								//etab for gauge
		etabb, 								//etabb for gauge
		KOep,								//factor for Kreiss-Oliger dissipation term
		inr, 								//inner radius of the matching region for the profile
		outr, 								//outer radius of the matching region for the profile
		mu, 								//amplitude of the perturbation
		kk, 								//scale of the perturbation
		mus, 								//amplitude of the perturbation
		kks, 								//scale of the perturbation
		Hb,									//initial Hubble
		fluidw, 							//w for fluid EOS
		Mkap,								//kappa for MUSCL
		bminmod,							//b for minmod
		hammax,								//maximum value of hamiltonian constraint violation
		mommax;								//maximum value of momentum constraint violation

		int exg;							//excision grid num

		int contn;							//1(true) for continue

		char file_continue[100]; 			//data file of geometry initial


	double ptintval1;						//1st part print interval boundary time
	double ptintval2;						//2nd
	double changept;						//printing interval change time
	double nexttimeprint;

	int horicheckintv,constoutintv;

	int exc=0;								//excision(0:false,1:true)
	int laymax;								//max number of layers
	int lbs[10];							//grid numbers for fmr region
	double alp_fmr[10];						//lapse values for starting fmr

	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);

	//reading initial file start
	ifstream fin("par_ini.d");				//open par file
	initial_read(fin,
		mstep,tmax, tab, amp,
		nzmax,
		xmax, zmax,
		cfl, cdt, etaa, etab, etabb,
		KOep, exg, contn, file_continue, inr,outr,mu,kk,mus,kks, Hb, fluidw, Mkap, bminmod,
		ptintval1,ptintval2,changept,horicheckintv,constoutintv);
	fin.close();


	ifstream finfmr("par_fmr.d");			//open par file
	initial_fmr(finfmr,
		laymax,lbs,alp_fmr);
	fin.close();

	cout << "mstep=" << mstep
		<< " cfl=" << cfl
		<< " KOep=" << KOep
		<< " exg=" << exg
		<< " mu=" << mu
		<< " kk=" << kk
		<< " mus=" << mus
		<< " kks=" << kks
		<< " changept=" << changept
		<<endl;


	//setting for bools
	fld=true;								// fluid evolution -> true/false
	scl=true;								// scalar evolution -> true/false
	cuev=false;								// curvature evaluation -> true/false

	//main class for bssn
	xmax=zmax/double(nzmax)/(1.+amp)*tab*xmax;
	Fmv *fmv=new Fmv(tab,tab,nzmax,xmax,zmax,amp,fld,scl,cuev);
	Fmv1 **fmv1;
	Fmv0 **fmv0;

	fmv0 = new Fmv0 *[laymax+1];
	fmv1 = new Fmv1 *[laymax];

	fmv0[0]=fmv;

	//initial parameter setting start
	double dr=fmv->get_dz();
	double dtini=cfl*dr;
	double lambda=0.;
	double scalarm=0.;
	double tini=2./(3.*(1+fluidw)*Hb);
	t=tini;

	fmv->initial_params(cfl,etaa,etab,etabb,lambda,dtini,dtini,dtini,t,t,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
	//initial parameter setting end

	//output files start
	ofstream filez("out_jkz.dat");							//as functions of z on (j,k)
	filez << "## nmax=" << nzmax << " zmax=" << zmax <<  endl;

	ofstream fileconst("out_const.dat");					//constrain evolution

	ofstream filehorizon[laymax+1];							//horizon evolution for each layer
	filehorizon[0].open("out_horizon_00.dat");				//horizon evolution for the first layer

	//ofstream fileneck[laymax+1];							//neck evolution for each layer
	//fileneck[0].open("out_neck_00.dat");					//neck evolution for the first layer

	ofstream fileall;										//for all variables to continue

	//ofstream fconv("out_conv.dat",ios::app);				//for convergence check
	//fconv.setf(ios_base::fixed, ios_base::floatfield);
	//fconv.precision(16);
	//output files end

	//reading continue or setting initial date start
	if(contn)
	{

		ifstream fcontinue(file_continue);
		if(fcontinue.fail()){
			cout << "Initial data file to continue is not found." << endl;
			exit(1);
		}

		//parameter consistency check
		check_continue_file(fcontinue,fld,scl,exc,exg,tab,nzmax,laymax,ln,lbs);

		fmv->base_initial_continue(fcontinue);
		// fmv->boundary_asym(0);
		// fmv->dyntoprim();

		cout << "continue" << endl;

		string buf;
		getline(fcontinue, buf);
		getline(fcontinue, buf);

		for(int i=0;i<ln;i++)
		{
			fmv0[i]->set_llmin(lbs[i]);
			xmax*=0.5;
			fmv1[i]=new Fmv1(tab,tab,2*lbs[i],xmax,fmv0[i]->get_z(lbs[i]),amp,fld, scl, cuev, fmv0[i]);
			fmv0[i+1]=fmv1[i];

			//initial setting for the upper layer start
			double deltat=0.5*fmv0[i]->get_dt0();
			fmv1[i]->initial_params(cfl,etaa,etab,etabb,lambda,deltat,deltat,deltat,t,tini,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
			fmv1[i]->fmr_initial_continue(fcontinue);
			// fmv1[i]->boundary_fmr();
			// fmv1[i]->dyntoprim();
			cout << "initial setting done" << endl;
			//initial setting for the upper layer end

			char buff[100];
			snprintf(buff,sizeof(buff),"out_horizon_%02d.dat",i+1);
			filehorizon[i+1].open(buff);
			//snprintf(buff,sizeof(buff),"out_neck_%02d.dat",i+1);
			//fileneck[i+1].open(buff);

			getline(fcontinue, buf);
			getline(fcontinue, buf);

			if(i>0)
			{
				fmv1[i-1]->set_mrf(true);
				fmv1[i-1]->set_ulay(fmv1[i]);
			}
		}

		//preparation for using excision
		if(exc)
		{
			fmv1[ln-1]->set_exc(true);
			fmv1[ln-1]->set_excflags();
		}

		fcontinue.close();

		t=fmv->get_t();
		fmv->print_z(filez,0,0);

		for(int i=0;i<ln;i++)
		{
			fmv1[i]->print_bz(filez,0,0);
		}
	}
	else
	{
		cout << "initial start" << endl;
		fmv->initial_longwave(mu,kk,inr,outr);
		//fmv->initial_iso_longwave(mus,kks,inr,outr);
		//fmv->print_z(filez,0,0);
		cout << "initial end" << endl;
		//exit(0);
	}
	//reading continue or setting initial date end

	cout << endl
	<< "          +---------------------------------------+" << endl
	<< "          |      Set up initial condition         |" << endl
	<< "          +---------------------------------------+" << endl << endl;

	cout << " time=" << fmv->get_t() << " : step=" << 0 << " " << endl;

	cout.setf(ios_base::scientific, ios_base::floatfield);
	cout.precision(10);

	//time step settings start
	double dtl=cfl*dr;
	double dtt=cdt*cfl/abs(fmv->get_bv(nzmax,0,0,20));

	cout << "dtt=" << dtt << endl;

	if(!contn)
	{
		fmv->set_dtpp(fmv->get_dt0());			//one before previous time step
		fmv->set_dtp(fmv->get_dt0());			//previous time step
		fmv->set_dt0(min(dtl,dtt));
	}
	//time step settings end

	//other settings for main loop
	if(t+ptintval1>changept)
	 nexttimeprint=int(t/ptintval2)*ptintval2;
	else
	 nexttimeprint=int(t/ptintval1)*ptintval1;

	//main loop start
	cout << "enter the main loop" << endl;

	bool changedt=false;						//if dt changed or not
	bool hform=false;							//horizon formation

	for(int step=1;step<mstep+1;step++)
	{
		cout << "step=" << step << endl;
		fmv->set01();				//xxx1=xxx0
		#pragma omp barrier
		fmv->setv0();				//xxx0=xxxv
		#pragma omp barrier
		fmv->set_zero_r();			//xxxr=0
		#pragma omp barrier

		//each step start
		fmv->BSSN(1);
		#pragma omp barrier

		fmv->boundary_asym(1);
		#pragma omp barrier

		//checking constraint start
		#pragma omp barrier
		fmv->check_const();
		#pragma omp barrier
		//fmv->print_const(fileconst);
		hammax=fmv->get_hammax();
		mommax=fmv->get_mommax();
		#pragma omp barrier
		cout << " layer number=0"
		<< " time=" << t
		<< " step=" << step+1 << endl
		<< " alp=" << fmv->get_bv(0,0,0,0)
		<< endl
		<< " ham=" << fmv->get_ham()
		<< " hammax=" << fmv->get_hammax()
		<< "  l="
		<< fmv->get_lhm()
		<< " r=" << fmv->get_z(fmv->get_lhm())
		<< endl
		<< " mom=" << fmv->get_mom()
		<< " mommax=" << fmv->get_mommax()
		<< "  l="
		<< fmv->get_lmm()
		<< " r=" << fmv->get_z(fmv->get_lmm())
		<< endl << endl;
		//checking constraint end

		fmv->BSSN(2);
		#pragma omp barrier
		fmv->boundary_asym(2);
		#pragma omp barrier

		fmv->BSSN(3);
		#pragma omp barrier
		fmv->boundary_asym(3);
		#pragma omp barrier

		fmv->BSSN(4);
		#pragma omp barrier
		//each step end
		fmv->boundary_asym(4);
		#pragma omp barrier

		fmv->set_dtpp(fmv->get_dtp());			//one before previous time step
		fmv->set_dtp(fmv->get_dt0());			//previous time step

		if(ln!=0)
		{
			fmv1[0]->evolve();
		}
		#pragma omp barrier

		//checking constraint in the highest layer start
		if(ln!=0)
		{
			//int lni=ln;
			for(int lni=1;lni<=ln;lni++)
			{
				hammax=max(hammax,fmv0[lni]->get_hammax());
				mommax=max(mommax,fmv0[lni]->get_mommax());

				cout << " layer number=" << lni << endl
				<< " ham=" << fmv0[lni]->get_ham()
				<< " hammax=" << fmv0[lni]->get_hammax()
				<< "  l="
				<< fmv0[lni]->get_lhm()
				<< " r=" << fmv0[lni]->get_z(fmv0[lni]->get_lhm())
				<< endl
				<< " mom=" << fmv0[lni]->get_mom()
				<< " mommax=" << fmv0[lni]->get_mommax()
				<< "  l="
				<< fmv0[lni]->get_lmm()
				<< " r=" << fmv0[lni]->get_z(fmv0[lni]->get_lmm())
				<< endl << endl;
			}
		}

		if(step%constoutintv==0)
		{
			fileconst << setw(20) << t					//1
			<< " "  << setw(20) << hammax				//2
			<< " "  << setw(20) << mommax				//3
			<< endl;
		}
		//checking constraint in the highest layer end

		//time forward start
		t=fmv->get_t()+fmv->get_dt0();
		fmv->set_t(t);
		#pragma omp barrier

		dtt=cdt*cfl/abs(fmv->get_bv(nzmax,0,0,20));
		fmv->set_dt0(min(dtl,dtt));

		if(changedt==false)
		{
			if(dtl<dtt)
			{
				cout << "dtl<dtt" << " t="<< t <<  endl;
				changedt=true;
			}
		}
		//time forward end

		//output judge start
		bool printflag=false;

		if(t>nexttimeprint-1.0e-10)
		{
			cout << "t=" << t << "nexttimeprint=" << nexttimeprint << endl;
			printflag=true;

			if(t+ptintval1>changept)
			nexttimeprint+=ptintval2;
			else
			nexttimeprint+=ptintval1;

			cout << "nexttimeprint=" << nexttimeprint << endl;
		}
		//output judge end

		//checking horizon formation (future outer trapped region) and excision start
		if(step%horicheckintv==0)
		{
			//NOTE: it might be found in not highest layer
			//      then the higher layer will be removed
			for(int i=ln;i>=0;i--)
			{
				fmv0[i]->check_horizon(filehorizon[i]);
				//fmv0[i]->check_neck(fileneck[i]);

				if(!hform)
				{
					if(fmv0[i]->get_hform())
					{
						hform=true;
						printflag=true;

						for(int ii=ln;ii>i;ii--)
						{
							delete fmv0[ii];
							cout << "fmv0-" << ii << "delete" << endl;
						}

						ln=i;							//reset the highest layer
						laymax=ln;

						if(i!=0)
						{
							fmv1[ln-1]->set_mrf(false);		//reset the mesh-refinement flag
							fmv1[ln-1]->set_llmin(fmv1[ln-1]->get_lli());
						}

						cout << "future trapped region found t="<< t << endl
						<< "layer number=" << i << endl;
					}

					if(fmv0[i]->get_exc())
					{
						if(!exc)
						{
							exc=true;

							cout << "excision start t="<< t << endl
							<< "excision layer number=" << i << endl;

						}
					}
				}
			}
		}
		//checking horizon formation (future outer trapped region) and excision end

		//escape judge
		if(t>tmax)
		 break;
		//if(hform)
		//break;

		//adding another layer if criterion is met start
		if(ln<laymax && exc==0)
		{
			//criterion (value of the lapse at the center)
			if(fmv->get_bv(0,0,0,0)<alp_fmr[ln])
			{
				// fmv->print_bz(filez,0,0);

				// for(int i=0;i<ln;i++)
				// {
				// 	fmv1[i]->print_bz(filez,0,0);
				// }

				fmv0[ln]->set_llmin(lbs[ln]);
				xmax*=0.5;
				fmv1[ln]=new Fmv1(tab,tab,2*lbs[ln],xmax,fmv0[ln]->get_z(lbs[ln]),amp,fld, scl, cuev, fmv0[ln]);
				fmv0[ln+1]=fmv1[ln];

				//initial setting for the upper layer start
				double deltat=0.5*fmv0[ln]->get_dt0();
				fmv1[ln]->initial_params(cfl,etaa,etab,etabb,lambda,deltat,deltat,deltat,t,tini,Hb,KOep,exg,fluidw,scalarm,Mkap,bminmod);
				fmv1[ln]->set_fmr_initial();
				cout << "initial setting done" << endl;
				//initial setting for the upper layer end

				char buff[100];
				snprintf(buff,sizeof(buff),"out_horizon_%02d.dat",ln+1);
				filehorizon[ln+1].open(buff);
				//snprintf(buff,sizeof(buff),"out_neck_%02d.dat",ln+1);
				//fileneck[ln+1].open(buff);
				cout << "initial setting done" << endl;

				if(ln>0)
				{
					fmv1[ln-1]->set_mrf(true);
					fmv1[ln-1]->set_ulay(fmv1[ln]);
				}
				ln++;
				cout << "fmr layer #" << ln << " start" << endl;
			}
		}
		//adding another layer if criterion is met end

		//printing all data files start
		if(printflag)
		{
			fileall.open("out_all.dat", ios::out );
			output_params(fileall,fld,scl,exc,exg,tab,nzmax,laymax,ln,lbs);

 			for(int i=0;i<=ln;i++)
			 fmv0[i]->print_all(fileall);
			fileall.close();

			for(int i=0;i<=ln;i++)
			{
				fmv0[i]->print_bz(filez,0,0);
			}
		}
		//printing all data files end
	}
	//main loop end

	//final print

	//final print
	if(mstep!=0)
	{
		for(int i=0;i<=ln;i++)
		{
			fmv0[i]->print_bz(filez,0,0);
		}
	}

	//convergence plot
	//fconv << setw(19) << fmv->get_dx() << " "  	/// 1
	//	<< setw(20) << fmv->get_bv(0,0,0,25) 		/// 2
	//	<< endl;

	//continue plot start
	fileall.open("out_all.dat", ios::out );
	output_params(fileall,fld,scl,exc,exg,tab,nzmax,laymax,ln,lbs);

	for(int i=0;i<=ln;i++)
	{
		fmv0[i]->setv0();
		fmv0[i]->print_all(fileall);
	}
	fileall.close();
	//continue plot end

	//finalize
	for(int i=0;i<=ln;i++)
	{
		delete fmv0[i];
		cout << "fmv0-"<< i << "delete" << endl;
	}

	delete[] fmv0;
	cout << "fmv0 delete" << endl;

	delete[] fmv1;
	cout << "fmv1 delete" << endl;
	//finalize end

	return 0;
}

/**
 * @brief Reads initial simulation parameters from a file ("par_ini.d").
 *
 * Parses the parameter file line by line, extracting values for simulation control,
 * grid setup, physics models, gauge conditions, numerical methods, initial data,
 * and output settings.
 *
 * @param fin Input file stream connected to "par_ini.d".
 * @param mstep Maximum number of time steps.
 * @param tmax Maximum simulation time.
 * @param tab Number of buffer grid points (ghost zones) in z-direction.
 * @param amp Amplitude parameter for inhomogeneous grid stretching.
 * @param nzmax Number of grid points in the z-direction for the base layer.
 * @param xmax Maximum coordinate value in x (related to zmax and grid setup).
 * @param zmax Maximum coordinate value in z for the base layer.
 * @param cfl CFL factor for time step calculation.
 * @param cdt Factor for cosmological time scale adjustment in dt calculation.
 * @param etaa Gauge parameter eta for lapse evolution.
 * @param etab Gauge parameter eta for shift evolution (related to Gamma driver).
 * @param etabb Gauge parameter eta for B^i evolution (Gamma driver).
 * @param KOep Coefficient for Kreiss-Oliger dissipation.
 * @param exg Grid point index for excision boundary (if enabled).
 * @param contn Flag indicating whether to continue from a previous run (1=yes, 0=no).
 * @param file_continue Name of the file to read data from for continuation.
 * @param inr Inner radius for matching initial perturbation profile.
 * @param outr Outer radius for matching initial perturbation profile.
 * @param mu Amplitude of the primary (adiabatic) perturbation.
 * @param kk Wavenumber/scale of the primary perturbation.
 * @param mus Amplitude of the scalar field (isocurvature) perturbation.
 * @param kks Wavenumber/scale of the scalar field perturbation.
 * @param Hb Initial Hubble parameter.
 * @param fluidw Equation of state parameter (w = p/rho) for the fluid.
 * @param Mkap Kappa parameter for MUSCL reconstruction (fluid dynamics).
 * @param bminmod Beta parameter for the minmod limiter (fluid dynamics).
 * @param ptintval1 Time interval for output during the first phase.
 * @param ptintval2 Time interval for output during the second phase.
 * @param changept Time at which the output interval changes from ptintval1 to ptintval2.
 * @param horicheckintv Step interval for checking apparent horizon formation.
 * @param constoutintv Step interval for outputting constraint violation information.
 */
void initial_read(ifstream& fin,
long int& mstep,
double& tmax,
int& tab,
double& amp,
int& nzmax,
double& xmax,
double& zmax,
double& cfl,
double& cdt,
double& etaa,
double& etab,
double& etabb,
double& KOep,
int& exg,
int& contn,
char *file_continue,
double& inr,
double& outr,
double& mu,
double& kk,
double& mus,
double& kks,
double& Hb,
double& fluidw,
double& Mkap,
double& bminmod,
double& ptintval1,
double& ptintval2,
double& changept,
int& horicheckintv,
int& constoutintv
){
	string buf;
	char cp[100];
	//cp = NULL;

	if(fin.fail()){
		cout << "Parameter File does not exist." << endl;
		abort();
	}

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%ld %s",&mstep,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&tmax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&tab,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&amp,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&nzmax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&xmax,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&zmax,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&cfl,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&cdt,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etaa,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etab,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&etabb,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&KOep,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&exg,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&contn,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%s %s",file_continue,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&inr,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&outr,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&mu,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&kk,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&mus,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&kks,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&Hb,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&fluidw,cp);
	}
	getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&Mkap,cp);
	}
	getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&bminmod,cp);
	}
	getline(fin, buf);

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ptintval1,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&ptintval2,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%lf %s",&changept,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&horicheckintv,cp);
	}
	getline(fin, buf);
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&constoutintv,cp);
	}
	return;
}


/**
 * @brief Reads Fixed Mesh Refinement (FMR) parameters from a file ("par_fmr.d").
 *
 * @param fin Input file stream connected to "par_fmr.d".
 * @param laymax Maximum number of FMR layers allowed (0 means no FMR).
 * @param lbs Array storing the grid index (in the parent layer's z-coordinate)
 *            defining the boundary of each refinement region.
 * @param alp_fmr Array storing the threshold values of the central lapse function (`alpha`)
 *                that trigger the creation of each new FMR layer.
 */
void initial_fmr(
ifstream& fin, 								//initial parameter file for fmr
int& laymax, 								//max fmr layer number
int *lbs, 									//grid number for fmr region
double *alp_fmr								//values of the lapse for starting fmr
){
	string buf;
	char cp[100];
	//cp = NULL;

	if(fin.fail()){
		cout << "Parameter File does not exist." << endl;
		abort();
	}

	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);
	getline(fin, buf);

	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&laymax,cp);
	}
	getline(fin, buf);

	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%d %s",&lbs[n],cp);
		}
	}
	getline(fin, buf);

	for(int n=0;n<=laymax;n++)
	{
		getline(fin, buf);
		if((buf[0]!='#')&&(!buf.empty())){
			sscanf(buf.c_str(),"%lf %s",&alp_fmr[n],cp);
		}
	}
	return;
}


/**
 * @brief Checks parameter consistency when continuing a simulation.
 *
 * Reads header information from the continuation file and compares key parameters
 * (like fluid/scalar flags, excision settings, grid setup) against the current
 * configuration read from parameter files. Exits if inconsistencies are found.
 *
 * @param fcontinue Input file stream connected to the continuation data file.
 * @param fld Current fluid evolution flag (to be checked).
 * @param scl Current scalar evolution flag (to be checked).
 * @param exc Excision flag read from the continuation file (output).
 * @param exg Current excision grid setting (to be checked against file).
 * @param tab Current buffer grid size (to be checked).
 * @param nzmax Current base layer z-grid size (to be checked).
 * @param laymax Current maximum allowed FMR layers (to be checked against file).
 * @param ln Number of FMR layers found in the continuation file (output).
 * @param lbs Array of FMR boundary indices (to be checked).
 */
void check_continue_file(					//checking parameter consistency with the continue file
ifstream& fcontinue, 						//initial parameter file
bool& fld,									//for fluid evolution
bool& scl,									//for scalar evolution
int& exc,									//excision
int& exg,									//excision grid
int& tab, 									//tab num for buf grids
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *lbs									//grid number for fmr region on z-axis
){
		string buf;
		//file preparateion end
		int cpar;
		//get the parameters of the continue file
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##fld=%d",&cpar);
		if(cpar!=fld)
		{
			cout << "fld is different from the initial data file" << endl
				<< "fld in data file=" << cpar << endl
				<< "fld in par_ini.d=" << tab << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##scl=%d",&cpar);
		if(cpar!=scl)
		{
			cout << "scl is different from the initial data file" << endl
				<< "scl in data file=" << cpar << endl
				<< "scl in par_ini.d=" << tab << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##exc=%d",&cpar);
		exc=cpar;
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##exg=%d",&cpar);
		if(cpar>exg)
		{
			cout << "exg is larger in initial data file" << endl
				<< "exg in data file=" << cpar << endl
				<< "exg in par_ini.d=" << exg << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##tab=%d",&cpar);
		if(cpar!=tab)
		{
			cout << "tab is different from the initial data file" << endl
				<< "tab in data file=" << cpar << endl
				<< "tab in par_ini.d=" << tab << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##nzmax=%d",&cpar);
		if(cpar!=nzmax)
		{
			cout << "nzmax is different from the initial data file" << endl
				<< "nzmax in data file=" << cpar << endl
				<< "nzmax in par_ini.d=" << tab << endl;
			exit(1);
		}
		getline(fcontinue, buf);
		cout << buf << endl;
		sscanf(buf.data(),"##ln=%d",&cpar);
		if(cpar>laymax)
		{
			cout << "ln is larger than laymax" << endl
				<< "ln in data file=" << cpar << endl
				<< "laymax in par_fmr.d=" << laymax << endl;
			exit(1);
		}
		ln=cpar;
		for(int i=0;i<ln;i++)
		{
			getline(fcontinue, buf);
			cout << buf << endl;
			sscanf(buf.data(),"##fmrzgnum=%d",&cpar);
			if(cpar!=lbs[i])
			{
				cout << i+1 << "-th fmrzgnum is different from the initial data file" << endl
					<< "fmrzgnum in data file=" << cpar << endl
					<< "fmrzgnum in par_ini.d=" << tab << endl;
			exit(1);
			}
		}
	return;
}


/**
 * @brief Writes essential parameters to the header of the continuation output file.
 *
 * This allows `check_continue_file` to verify consistency when restarting.
 *
 * @param fileall Output file stream for the continuation data ("out_all.dat").
 * @param fld Fluid evolution flag.
 * @param scl Scalar evolution flag.
 * @param exc Excision flag (0 or 1).
 * @param exg Excision grid index.
 * @param tab Buffer grid size.
 * @param nzmax Base layer z-grid size.
 * @param laymax Maximum allowed FMR layers.
 * @param ln Current number of active FMR layers.
 * @param lbs Array of FMR boundary indices for active layers.
 */
void output_params(							//output all needed parameters for continue
ofstream& fileall, 							//initial parameter file
bool& fld,									//for fluid evolution
bool& scl,									//for scalar evolution
int& exc,									//excision
int& exg,									//excision grid
int& tab, 									//tab num for buf grids
int& nzmax, 								//max grid num of z
int& laymax,								//max fmr layer number
int& ln,									//fmr layer number
int *lbs									//grid number for fmr region on z-axis
){
	fileall.setf(ios_base::scientific, ios_base::floatfield);
	fileall.precision(10);
	fileall
	<< "##fld="<< fld << endl
	<< "##scl="<< scl << endl
	<< "##exc="<< exc << endl
	<< "##exg="<< exg << endl
	<< "##tab="<< tab << endl
	<< "##nzmax="<< nzmax << endl
	<< "##ln="<< ln << endl;
	for(int i=0;i<ln;i++)
	fileall << "##fmrzgnum="<< lbs[i] << endl;

}
