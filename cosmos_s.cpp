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
#include "cosmos_s.h"						//header for Einstein solver

using namespace std;



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
int& exc,									//initial excision setting
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
double& changept							//printing interval change time
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
		bminmod;							//b for minmod 
		
		int exg;							//excision grid num
		
		int contn;							//1(true) for continue
		
		char file_continue[100]; 			//data file of geometry initial

		
	double ptintval1;						//1st part print interval boundary time
	double ptintval2;						//2nd
	double changept;						//printing interval change time
	double nexttimeprint;
	
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
		KOep, exg, exc, contn, file_continue, inr,outr,mu,kk,mus,kks, Hb, fluidw, Mkap, bminmod,
		ptintval1,ptintval2,changept);
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
	Fmv *fmv=new Fmv(tab,tab,nzmax,xmax,zmax,amp,fld,scl,cuev,1);
	Fmv1 **fmv1;
	Fmv0 **fmv0;
	
	fmv0 = new Fmv0 *[laymax+1];
	fmv1 = new Fmv1 *[laymax];
	
	fmv0[0]=fmv;
	
	//initial parameter setting start
	double dr=fmv->get_dz();
	double dtini=cfl*dr;
	fmv->set_cfl(cfl);
	fmv->set_etaa(etaa);
	fmv->set_etab(etab);
	fmv->set_etabb(etabb);
	fmv->set_lambda(0.);
	fmv->set_tmax(tmax);
	fmv->set_dt0(dtini);
	fmv->set_dtp(dtini);
	fmv->set_dtpp(dtini);
	fmv->set_fluidw(fluidw);
	t=2./(3.*(1+fluidw)*Hb);
	fmv->set_t(t);
	fmv->set_Hb(Hb);
	fmv->set_tini(t);
	fmv->set_KOep(KOep);
	fmv->set_exg(exg);
	fmv->set_Mkap(Mkap);
	fmv->set_b(bminmod);
	//initial parameter setting end
	
	
	//output files start
	ofstream filez("out_jkz.dat");							//as functions of z on (j,k)
	filez << "## nmax=" << nzmax << " zmax=" << zmax <<  endl;

	ofstream fileconst("out_const.dat");					//constrain evolution

	ofstream filehorizon("out_horizon.dat");				//horizon evolution
	
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
			abort();
		}
		
		string buf;
		getline(fcontinue, buf);
		cout << buf << endl;
		
		sscanf(buf.data(),"##ln=%d",
				&ln);

		cout << "continue" << endl;
		fmv->initial_continue(fcontinue);
		getline(fcontinue, buf);
		getline(fcontinue, buf);
		
		for(int i=0;i<ln;i++)
		{
			fmv0[i]->set_llmin(lbs[i]);
			xmax*=0.5;
			fmv1[i]=new Fmv1(tab,tab,2*lbs[i],xmax,fmv0[i]->get_z(lbs[i]),amp,fld, scl, cuev, fmv0[i]);
			fmv0[i+1]=fmv1[i];
			fmv1[i]->set_fmr_initial();
			cout << "initial setting done" << endl;
			
			fmv1[i]->initial_continue(fcontinue);
			getline(fcontinue, buf);
			getline(fcontinue, buf);

			if(i>0)
			{
				fmv1[i-1]->set_mrf(true);
				fmv1[i-1]->set_ulay(fmv1[i]);
			}
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

	//preparation for using excision
	if(exc)
	{
		fmv->set_exc(true);
		fmv->set_excflags();
	}

	//main loop start
	cout << "enter the main loop" << endl;
	
	bool changedt=false;						//if dt changed or not
	bool hform=false;							//horizon formation
	bool excf=false;							//excision is acting or not
	
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
		fmv->print_const(fileconst);
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
		
		//escape judge 
		if(t>tmax)
		{
			step=mstep+1;
			//fmv->print_z(filez,0,0);
			break;
		}	
		
		
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

		//printing 
		if(printflag)
		{
			fmv->print_z(filez,0,0);
			
			for(int i=0;i<ln;i++)
			{
				fmv1[i]->print_bz(filez,0,0);
			}
		}

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
		
		if(ln!=0)
		{
			int lni=ln;
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
		
		//checking horizon formation and excision
		if(!hform || step%100==0)
		{
			//NOTE: it might be found in not highest layer
			//      then the higher layer will be removed
			for(int i=ln;i>=0;i--)
			{
				fmv0[i]->check_horizon(filehorizon);
			
				if(fmv0[i]->get_hform())
				{
					if(!hform)
					{
						fileall.open(file_continue, ios::out );
						fileall.setf(ios_base::scientific, ios_base::floatfield);
						fileall.precision(10);
						
						for(int i=0;i<=ln;i++)
						{	
							fmv0[i]->setv0();
							fmv0[i]->print_all(fileall);
						}
									
						fileall.close();
						hform=true;
					}
					break;
				}
			}

			if(!excf)
			{
				for(int i=ln;i>=0;i--)
				{
					//if horizon has been found and excision has been done 
					if(fmv0[i]->get_exc())
					{
						excf=true;
						ln=i;							//reset the highest layer

						if(i!=0)
						fmv1[ln-1]->set_mrf(false);		//reset the mesh-refinement flag
					}
				}
			}
		}

		//adding another layer if criterion is met
		if(ln<laymax && !excf)
		{	
			//criterion (value of the laps at the center
			if(fmv->get_bv(0,0,0,0)<alp_fmr[ln])
			{
				fmv->print_bz(filez,0,0);
			
				for(int i=0;i<ln;i++)
				{
					fmv1[i]->print_bz(filez,0,0);
				}

				fmv0[ln]->set_llmin(lbs[ln]);
				xmax*=0.5;
				fmv1[ln]=new Fmv1(tab,tab,2*lbs[ln],xmax,fmv0[ln]->get_z(lbs[ln]),amp,fld, scl, cuev, fmv0[ln]);
				fmv0[ln+1]=fmv1[ln];
				fmv1[ln]->set_fmr_initial();
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
		
		if(printflag)
		{
			fileall.open(file_continue, ios::out );
 			fileall.setf(ios_base::scientific, ios_base::floatfield);
			fileall.precision(10);
			fileall << "##ln=" << ln << endl;
 			cout <<"##ln=" << ln << endl;
 			for(int i=0;i<=ln;i++)
			 fmv0[i]->print_all(fileall);
			fileall.close();
		}
	}
	//main loop end
	
	//final print 
	if(mstep!=0)
	{
		fmv->setv0();
		fmv->print_z(filez,0,0);
	}
	
	//convergence plot 
	//fconv << setw(19) << fmv->get_dx() << " "  	/// 1
	//	<< setw(20) << fmv->get_bv(0,0,0,25) 		/// 2
	//	<< endl;
	
	//continue plot start
	fileall.open(file_continue, ios::out );
	fileall.setf(ios_base::scientific, ios_base::floatfield);
	fileall.precision(10);
	for(int i=0;i<=ln;i++)
	{	
		fmv0[i]->setv0();
		fmv0[i]->print_all(fileall);
	}
	fileall.close();
	//continue plot end
		
	//finalize
	delete fmv;
	
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
int& exc, 
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
double& changept
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
	if((buf[0]!='#')&&(!buf.empty())){
		sscanf(buf.c_str(),"%d %s",&exc,cp);
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

