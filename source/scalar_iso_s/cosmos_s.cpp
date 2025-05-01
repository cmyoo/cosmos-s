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
#include "../cosmos_s.h"						//header for Einstein solver

using namespace std;

void check_jkz(Fmv *fmv);

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
		//fmv->initial_longwave(mu,kk,inr,outr);
		fmv->initial_iso_longwave(mus,kks,inr,outr);
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
		if(step%horicheckintv==1)
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
		// if(hform)
		// break;


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

	ofstream filecheck("out_check_jkz.dat");						//as functions of x on (k,l) for check
	filecheck << "## nmax=" << nzmax << " zmax=" << zmax <<  endl;
	fmv->print_bz(filecheck,0,0);

	//comparing out_jkz.dat with exp_jkz.dat 
	check_jkz(fmv);

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

//function to check difference of out_xkl.dat from exp_xkl.dat
void check_jkz(Fmv *fmv)
{
	ifstream expfile("exp_jkz.dat");				//open exp file
	if(expfile.fail()){
		cout << "exp_jkz.dat is not found." << endl;
		abort();
	}

	ifstream outfile("out_check_jkz.dat");				//open out file
	if(outfile.fail()){
		cout << "out_check_jkz.dat is not found." << endl;
		abort();
	}

	string bufexp,bufout;

	cout << "header texts in exp_jkz.dat" << endl;
	getline(expfile, bufexp);
	cout << bufexp << endl;
	getline(expfile, bufexp);
	cout << bufexp << endl;
	
	cout << "header texts in out_check_jkz.dat" << endl;
	getline(outfile, bufout);
	cout << bufout << endl;
	getline(outfile, bufout);
	cout << bufout << endl;

	double vdiff[50];
	for(int n=0;n<50;n++)
	vdiff[n]=0.;

	for(int l=fmv->get_llmin();l<=fmv->get_lmax();l++)
	{
		double vexp[50];
		double vout[50];

		getline(outfile, bufout);
		getline(expfile, bufexp);
		
		if(fmv->get_fluidevo() && fmv->get_scalarevo())
		{
			sscanf(bufexp.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vexp[0],&vexp[1],&vexp[2],&vexp[3],&vexp[4],&vexp[5],&vexp[6],&vexp[7],&vexp[8],&vexp[9],
			&vexp[10],&vexp[11],&vexp[12],&vexp[13],&vexp[14],&vexp[15],&vexp[16],&vexp[17],&vexp[18],&vexp[19],
			&vexp[20],&vexp[21],&vexp[22],&vexp[23],&vexp[24],&vexp[25],&vexp[26],&vexp[27],&vexp[28],&vexp[29],
			&vexp[30],&vexp[31],&vexp[32],&vexp[33],&vexp[34],&vexp[35],&vexp[36],&vexp[37],&vexp[38],&vexp[39],
			&vexp[40],&vexp[41],&vexp[42],&vexp[43],&vexp[44],&vexp[45],&vexp[46],&vexp[47],&vexp[48],&vexp[49]);

			sscanf(bufout.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vout[0],&vout[1],&vout[2],&vout[3],&vout[4],&vout[5],&vout[6],&vout[7],&vout[8],&vout[9],
			&vout[10],&vout[11],&vout[12],&vout[13],&vout[14],&vout[15],&vout[16],&vout[17],&vout[18],&vout[19],
			&vout[20],&vout[21],&vout[22],&vout[23],&vout[24],&vout[25],&vout[26],&vout[27],&vout[28],&vout[29],
			&vout[30],&vout[31],&vout[32],&vout[33],&vout[34],&vout[35],&vout[36],&vout[37],&vout[38],&vout[39],
			&vout[40],&vout[41],&vout[42],&vout[43],&vout[44],&vout[45],&vout[46],&vout[47],&vout[48],&vout[49]);
			
			for(int n=0;n<50;n++)
			{
				vdiff[n]+=abs(vexp[n]-vout[n]);
			}
		}
		else if(fmv->get_fluidevo())
		{
			sscanf(bufexp.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vexp[0],&vexp[1],&vexp[2],&vexp[3],&vexp[4],&vexp[5],&vexp[6],&vexp[7],&vexp[8],&vexp[9],
			&vexp[10],&vexp[11],&vexp[12],&vexp[13],&vexp[14],&vexp[15],&vexp[16],&vexp[17],&vexp[18],&vexp[19],
			&vexp[20],&vexp[21],&vexp[22],&vexp[23],&vexp[24],&vexp[25],&vexp[26],&vexp[27],&vexp[28],&vexp[29],
			&vexp[30],&vexp[31],&vexp[32],&vexp[33],&vexp[34],&vexp[35],&vexp[36],&vexp[37],&vexp[38],&vexp[39],
			&vexp[40],&vexp[41],&vexp[42],&vexp[43],&vexp[44],&vexp[45]);

			sscanf(bufout.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vout[0],&vout[1],&vout[2],&vout[3],&vout[4],&vout[5],&vout[6],&vout[7],&vout[8],&vout[9],
			&vout[10],&vout[11],&vout[12],&vout[13],&vout[14],&vout[15],&vout[16],&vout[17],&vout[18],&vout[19],
			&vout[20],&vout[21],&vout[22],&vout[23],&vout[24],&vout[25],&vout[26],&vout[27],&vout[28],&vout[29],
			&vout[30],&vout[31],&vout[32],&vout[33],&vout[34],&vout[35],&vout[36],&vout[37],&vout[38],&vout[39],
			&vout[40],&vout[41],&vout[42],&vout[43],&vout[44],&vout[45]);
			
			for(int n=0;n<46;n++)
			{
				vdiff[n]+=abs(vexp[n]-vout[n]);
			}
		}
		else
		{
			sscanf(bufexp.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vexp[0],&vexp[1],&vexp[2],&vexp[3],&vexp[4],&vexp[5],&vexp[6],&vexp[7],&vexp[8],&vexp[9],
			&vexp[10],&vexp[11],&vexp[12],&vexp[13],&vexp[14],&vexp[15],&vexp[16],&vexp[17],&vexp[18],&vexp[19],
			&vexp[20],&vexp[21],&vexp[22],&vexp[23],&vexp[24],&vexp[25],&vexp[26],&vexp[27],&vexp[28],&vexp[29],
			&vexp[30],&vexp[31],&vexp[32],&vexp[33],&vexp[34],&vexp[35],&vexp[36],&vexp[37],&vexp[38],&vexp[39]);

			sscanf(bufout.data(),"%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf",
			&vout[0],&vout[1],&vout[2],&vout[3],&vout[4],&vout[5],&vout[6],&vout[7],&vout[8],&vout[9],
			&vout[10],&vout[11],&vout[12],&vout[13],&vout[14],&vout[15],&vout[16],&vout[17],&vout[18],&vout[19],
			&vout[20],&vout[21],&vout[22],&vout[23],&vout[24],&vout[25],&vout[26],&vout[27],&vout[28],&vout[29],
			&vout[30],&vout[31],&vout[32],&vout[33],&vout[34],&vout[35],&vout[36],&vout[37],&vout[38],&vout[39]);
			
			for(int n=0;n<40;n++)
			{
				vdiff[n]+=abs(vexp[n]-vout[n]);
			}
		}	
	}
	
	int gridnum=fmv->get_lmax()-fmv->get_llmin()+1;

	ofstream fdiff("out_diff.dat");
	fdiff.precision(16);
	fdiff.setf(ios_base::scientific, ios_base::floatfield);

	fdiff << "## difference between out_check_xkl.dat and exp_xkl.dat" << endl;
	fdiff 
	<< "diff x=" <<		vdiff[0]/gridnum << endl		//1
	<< "diff alp=" <<  	vdiff[1]/gridnum << endl		//2
	<< "diff by=" <<		vdiff[2]/gridnum << endl		//4
	<< "diff bz=" <<		vdiff[3]/gridnum << endl		//5
	<< "diff bx=" <<   	vdiff[4]/gridnum << endl		//3
	<< "diff bbx=" <<  	vdiff[5]/gridnum << endl		//6
	<< "diff bby=" <<	vdiff[6]/gridnum << endl		//7
	<< "diff bbz=" <<	vdiff[7]/gridnum << endl		//8
	<< "diff gxx=" <<	vdiff[8]/gridnum << endl		//9
	<< "diff gyy=" <<	vdiff[9]/gridnum << endl		//10
	<< "diff gzz=" <<	vdiff[10]/gridnum << endl		//11
	<< "diff gxy=" <<	vdiff[11]/gridnum << endl		//12
	<< "diff gxz=" <<	vdiff[12]/gridnum << endl		//13
	<< "diff gyz=" <<	vdiff[13]/gridnum << endl		//14
	<< "diff wa=" <<		vdiff[14]/gridnum << endl		//15
	<< "diff akxx=" <<	vdiff[15]/gridnum << endl		//16
	<< "diff akyy=" <<	vdiff[16]/gridnum << endl		//17
	<< "diff akzz=" <<	vdiff[17]/gridnum << endl		//18
	<< "diff akxy=" <<	vdiff[18]/gridnum << endl		//19
	<< "diff akxz=" <<	vdiff[19]/gridnum << endl		//20
	<< "diff akyz=" <<	vdiff[20]/gridnum << endl		//21
	<< "diff ek=" <<	vdiff[21]/gridnum << endl		//22
	<< "diff zgx=" <<	vdiff[22]/gridnum << endl		//23
	<< "diff zgy=" <<	vdiff[23]/gridnum << endl		//24
	<< "diff zgz=" <<	vdiff[24]/gridnum << endl		//25
	<< "diff hamn=" <<	vdiff[25]/gridnum << endl		//26
	<< "diff ham=" <<	vdiff[26]/gridnum << endl		//27
	<< "diff nM_z=" <<	vdiff[27]/gridnum << endl		//28
	<< "diff M_z=" <<	vdiff[28]/gridnum << endl		//29
	<< "diff dGamz=" <<	vdiff[29]/gridnum << endl		//30
	<< "diff Mass=" <<	vdiff[30]/gridnum << endl		//31
	<< "diff Kinv=" <<	vdiff[31]/gridnum << endl		//32
	<< "diff Arad=" <<	vdiff[32]/gridnum << endl		//33
	<< "diff Comp=" <<	vdiff[33]/gridnum << endl		//34
	<< "diff null_exp_p=" <<	vdiff[34]/gridnum << endl		//35
	<< "diff null_exp_m=" <<	vdiff[35]/gridnum << endl;		//36
	
	if(fmv->get_fluidevo())
	{
		fdiff 
		<<"diff Ene=" << vdiff[36]/gridnum << endl		//37
		<<"diff px=" << vdiff[37]/gridnum << endl		//38
		<<"diff py=" << vdiff[38]/gridnum << endl		//39
		<<"diff pz=" << vdiff[39]/gridnum << endl		//40
		<<"diff Den=" << vdiff[40]/gridnum << endl		//41
		<<"diff rho/rhob=" << vdiff[41]/gridnum << endl		//42
		<<"diff Vx=" << vdiff[42]/gridnum << endl		//43
		<<"diff Vy=" << vdiff[43]/gridnum << endl		//44
		<<"diff Vz=" << vdiff[44]/gridnum << endl		//45
		<<"diff eps=" << vdiff[45]/gridnum << endl;		//46

		if(fmv->get_scalarevo())
		{
			fdiff 
			<<"diff phii=" << vdiff[46]/gridnum << endl	//47
			<<"diff phii=" << vdiff[47]/gridnum << endl	//48
			<<"diff phii=" << vdiff[48]/gridnum << endl	//49
			<<"diff Pi=" << vdiff[49]/gridnum << endl;	//50
		}
	}
	else if(fmv->get_scalarevo())
	{
		fdiff 
		<<"diff phii=" << vdiff[36]/gridnum << endl	//37
		<<"diff phii=" << vdiff[37]/gridnum << endl	//38
		<<"diff phii=" << vdiff[38]/gridnum << endl	//39
		<<"diff Pi=" << vdiff[39]/gridnum << endl;	//40
	}

	fdiff << "## diff file end" << endl;

	return;
}
