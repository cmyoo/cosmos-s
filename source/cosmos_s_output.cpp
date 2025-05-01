/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* OUTPUT FUNCTIONS :: BSSN evolution Class of COSMOS_S                                                  */
/*                                             ver. 1.00           coded by Chulmoon Yoo                 */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos_s.h"
#include <iomanip>
#include <fstream>

//output as functions of z, j:xgrid num, k:y grid num
void Fmv0::print_z(ofstream& fout, int j, int k)
{
	if(j<jmin)
	{
		fout << "error : j<jmin" << endl;
		return;

	}
	else if(j>jmax)
	{
		fout << "error : j>jmax" << endl;
		return;
	}
	else if(k<kmin)
	{
		fout << "error : k<kmin" << endl;
		return;
	}
	else if(k>kmax)
	{
		fout << "error : k>kmax" << endl;
		return;
	}

	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##layer number=" << layn 
		<< "time=" << tt 
		<< " trK="  << get_bv0(lui,kui,jui,20)
		<< " j="  << j
		<< " k="  << k
		<< endl;

	double psibg=1./(3.*(1.+fluidw))*log(tt/tini);
	double trkbg=-2./(1.+fluidw)/tt;
	double rhobg=1./(24.*M_PI)*trkbg*trkbg;
	
	for(int l=llmin;l<=lmax;l++)
	//for(int l=llmin;l<=lui;l++)
	//for(int l=lmin;l<=lui;l++)
	{
		if(get_bflag(l,k,j)!=0)
		continue;

		double z = get_z(l);
		
		/////////////// output for geometry start ////////////////////////////////////////////////
		double gxx,gyy,gzz,gxy,gxz,gyz,wa;														//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double M_z,nM_z,dGamz;
		double alp,hamn,ham;
		double Kinv;
		double Arad;
		// double dM,dR;
		double null_exp_p,null_exp_m;
		// double eneratio;
		double Mass,Comp;
		
		alp= get_bv0(l,k,j, 0);
		bx =  get_bv0(l,k,j, 1);
		by =  get_bv0(l,k,j, 2);
		bz =  get_bv0(l,k,j, 3);
		bbx=  get_bv0(l,k,j, 4);
		bby=  get_bv0(l,k,j, 5);
		bbz=  get_bv0(l,k,j, 6);
		gxx=  get_bv0(l,k,j, 7)+1.;
		gyy=  get_bv0(l,k,j, 8)+1.;
		gzz=  get_bv0(l,k,j, 9)+get_flat_df2z(l);
		gxy=  get_bv0(l,k,j,10);
		gxz=  get_bv0(l,k,j,11);
		gyz=  get_bv0(l,k,j,12);
		wa=   get_bv0(l,k,j,13);
		akxx= get_bv0(l,k,j,14);
		akyy= get_bv0(l,k,j,15);
		akzz= get_bv0(l,k,j,16);
		akxy= get_bv0(l,k,j,17);
		akxz= get_bv0(l,k,j,18);
		akyz= get_bv0(l,k,j,19);
		ek=   get_bv0(l,k,j,20);
		zgx=  get_bv0(l,k,j,21);
		zgy=  get_bv0(l,k,j,22);
		zgz=  get_bv0(l,k,j,23);

				
		hamn  = get_con(l,k,j,0);
		ham = get_con(l,k,j,1);
		nM_z=get_con(l,k,j,2);
		M_z=get_con(l,k,j,3);
		dGamz=get_con(l,k,j,4);
		
		Kinv=get_outv(l,k,j,0);
		
		//dM=get_outv(l,k,j,4);
		//dR=get_outv(l,k,j,6);
		null_exp_p=get_outv(l,k,j,4);
		null_exp_m=get_outv(l,k,j,6);
		Mass=get_outv(l,k,j,9);
		
		if(l<=lui)
		Arad=get_outv(l,k,j,5);
		else
		Arad=get_outv(lui,k,j,5);
		
		
		if(Arad<=1.0e-15)
		Comp=0.;
		else
		Comp=2.*Mass/Arad;
		
		// if(l>lui)
		// eneratio=0.;
		// else
		// eneratio=(dM/dR)/(pi4*pow(Arad,2)*rhobg)-1.;
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << z 						//1
		<< " "  << setw(20) << alp  				//2
		<< " "  << setw(20) << bx   				//3
		<< " "  << setw(20) << by					//4
		<< " "  << setw(20) << bz					//5
		<< " "  << setw(20) << bbx  				//6
		<< " "  << setw(20) << bby					//7
		<< " "  << setw(20) << bbz					//8
		<< " "  << setw(20) << gxx					//9
		<< " "  << setw(20) << gyy					//10
		<< " "  << setw(20) << gzz					//11
		<< " "  << setw(20) << gxy					//12
		<< " "  << setw(20) << gxz					//13
		<< " "  << setw(20) << gyz					//14
		<< " "  << setw(20) << exp(wa)/exp(psibg)	//15
		<< " "  << setw(20) << akxx					//16
		<< " "  << setw(20) << akyy					//17
		<< " "  << setw(20) << akzz					//18
		<< " "  << setw(20) << akxy					//19
		<< " "  << setw(20) << akxz					//20
		<< " "  << setw(20) << akyz					//21
		<< " "  << setw(20) << ek					//22
		<< " "  << setw(20) << zgx					//23
		<< " "  << setw(20) << zgy					//24
		<< " "  << setw(20) << zgz					//25
		<< " "  << setw(20) << hamn					//26
		<< " "  << setw(20) << ham					//27
		<< " "  << setw(20) << nM_z					//28
		<< " "  << setw(20) << M_z					//29
		<< " "  << setw(20) << dGamz				//30
		<< " "  << setw(20) << Mass					//31
		<< " "  << setw(20) << Kinv					//32
		<< " "  << setw(20) << Arad					//33
		<< " "  << setw(20) << Comp					//34
		<< " "  << setw(20) << null_exp_p			//35
		<< " "  << setw(20) << null_exp_m;			//36
																								//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)																			//
		{
			double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
			
			Ene=  get_bv0(l,k,j,24);
			px=   get_bv0(l,k,j,25);
			py=   get_bv0(l,k,j,26);
			pz=   get_bv0(l,k,j,27);
			Den=  get_bv0(l,k,j,28);
			
			rho=  get_primv(l,k,j,0);
			Vx=   get_primv(l,k,j,1);
			Vy=   get_primv(l,k,j,2);
			Vz=   get_primv(l,k,j,3);
			eps=  get_primv(l,k,j,4);
			
			fout 
			<< " "  << setw(20) << Ene			//37
			<< " "  << setw(20) << px			//38
			<< " "  << setw(20) << py			//39
			<< " "  << setw(20) << pz			//40
			<< " "  << setw(20) << Den			//41
			<< " "  << setw(20) << rho/rhobg	//42
			<< " "  << setw(20) << Vx			//43
			<< " "  << setw(20) << Vy			//44
			<< " "  << setw(20) << Vz			//45
			<< " "  << setw(20) << eps;			//46
		}																						//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)																			//
		{
			double phii=get_bv0(l,k,j,nsc);
			double Pi=get_bv0(l,k,j,nscp);
			//double sEne=get_outv(l,k,j,7);
			double rEne=get_outv(l,k,j,10);
			double sp_z=get_outv(l,k,j,8);
			
			fout
			<< " "  << setw(20) << phii				//47	//37
			<< " "  << setw(20) << Pi				//48	//38
			<< " "  << setw(20) << rEne				//49	//39
			<< " "  << setw(20) << sp_z;			//50	//40
		}																						//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of z, j:xgrid num, k:y grid num
void Fmv0::print_ipolz(ofstream& fout)
{
	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##layer number=" << layn 
		<< "time=" << tt 
		<< " trK="  << get_bv0(lui,kui,jui,20)
		<< endl;

	double psibg=1./(3.*(1.+fluidw))*log(tt/tini);
	
	//for(int l=lmin;l<=lmax;l++)
	//for(int l=lli;l<=lui;l++)
	double z=get_z(lmin);
	
	#pragma omp parallel for
	for(int i=0;i<nn;i++)
	put_ddy(i);
	
	for(int l=lmin;l<=2*lui;l++,z+=0.5*dz)
	{
		/////////////// output for geometry start ////////////////////////////////////////////////
		double gxx,gyy,gzz,gxy,gxz,gyz,wa;							//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double alp;
		
		int lc=int(l/2);
		//alp= get_bv0(l,k,j, 0);
		alp= bv_spline_ipol(lc,z,0);
		bx = bv_spline_ipol(lc,z, 1);
		by = bv_spline_ipol(lc,z, 2);
		bz = bv_spline_ipol(lc,z, 3);
		bbx= bv_spline_ipol(lc,z, 4);
		bby= bv_spline_ipol(lc,z, 5);
		bbz= bv_spline_ipol(lc,z, 6);
		gxx= bv_spline_ipol(lc,z, 7)+1.;
		gyy= bv_spline_ipol(lc,z, 8)+1.;
		gzz= bv_spline_ipol(lc,z, 9)+get_flat_df2z(l);
		gxy= bv_spline_ipol(lc,z,10);
		gxz= bv_spline_ipol(lc,z,11);
		gyz= bv_spline_ipol(lc,z,12);
		wa=  bv_spline_ipol(lc,z,13);
		akxx=bv_spline_ipol(lc,z,14);
		akyy=bv_spline_ipol(lc,z,15);
		akzz=bv_spline_ipol(lc,z,16);
		akxy=bv_spline_ipol(lc,z,17);
		akxz=bv_spline_ipol(lc,z,18);
		akyz=bv_spline_ipol(lc,z,19);
		ek=  bv_spline_ipol(lc,z,20);
		zgx= bv_spline_ipol(lc,z,21);
		zgy= bv_spline_ipol(lc,z,22);
		zgz= bv_spline_ipol(lc,z,23);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << z 						//1
		<< " "  << setw(20) << alp  				//2
		<< " "  << setw(20) << bx   				//3
		<< " "  << setw(20) << by					//4
		<< " "  << setw(20) << bz					//5
		<< " "  << setw(20) << bbx  				//6
		<< " "  << setw(20) << bby					//7
		<< " "  << setw(20) << bbz					//8
		<< " "  << setw(20) << gxx					//9
		<< " "  << setw(20) << gyy					//10
		<< " "  << setw(20) << gzz					//11
		<< " "  << setw(20) << gxy					//12
		<< " "  << setw(20) << gxz					//13
		<< " "  << setw(20) << gyz					//14
		<< " "  << setw(20) << exp(wa)/exp(psibg)	//15
		<< " "  << setw(20) << akxx					//16
		<< " "  << setw(20) << akyy					//17
		<< " "  << setw(20) << akzz					//18
		<< " "  << setw(20) << akxy					//19
		<< " "  << setw(20) << akxz					//20
		<< " "  << setw(20) << akyz					//21
		<< " "  << setw(20) << ek					//22
		<< " "  << setw(20) << zgx					//23
		<< " "  << setw(20) << zgy					//24
		<< " "  << setw(20) << zgz;					//25
																								//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)																			//
		{
			double Ene,px,py,pz,Den;
			
			Ene=  bv_spline_ipol(lc,z,24);
			px=   bv_spline_ipol(lc,z,25);
			py=   bv_spline_ipol(lc,z,26);
			pz=   bv_spline_ipol(lc,z,27);
			Den=  bv_spline_ipol(lc,z,28);
			
			fout 
			<< " "  << setw(20) << Ene		//36
			<< " "  << setw(20) << px		//37
			<< " "  << setw(20) << py		//38
			<< " "  << setw(20) << pz		//39
			<< " "  << setw(20) << Den;		//40
		}																						//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)																			//
		{
			double phii=bv_spline_ipol(lc,z,nsc);
			double Pi=bv_spline_ipol(lc,z,nscp);
			
			fout
			<< " "  << setw(20) << phii			//46	//36
			<< " "  << setw(20) << Pi;			//47	//37
		}																						//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

//output as functions of z, j:xgrid num, k:y grid num
void Fmv0::print_bz(ofstream& fout, int j, int k)
{

	if(j<jmin)
	{
		fout << "error : j<jmin" << endl;
		return;

	}
	else if(j>jmax)
	{
		fout << "error : j>jmax" << endl;
		return;
	}
	else if(k<kmin)
	{
		fout << "error : k<kmin" << endl;
		return;
	}
	else if(k>kmax)
	{
		fout << "error : k>kmax" << endl;
		return;
	}

	double tt = get_t();

	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(8);
	fout << "##layer number=" << layn 
		<< " time=" << tt 
		<< " trK="  << get_bv(lui,kui,jui,20)
		<< " j="  << j
		<< " k="  << k
		<< endl;

	double psibg=1./(3.*(1.+fluidw))*log(tt/tini);
	double trkbg=-2./(1.+fluidw)/tt;
	double rhobg=1./(24.*M_PI)*trkbg*trkbg;
	
	//for(int l=lmin;l<=lmax;l++)
	//for(int l=lli;l<=lui;l++)
	for(int l=llmin;l<=lui;l++)
	{
		if(get_bflag(l,k,j)!=0)
		continue;

		double z = get_z(l);
		
		/////////////// output for geometry start ////////////////////////////////////////////////
		double gxx,gyy,gzz,gxy,gxz,gyz,wa;														//
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double M_z,nM_z,dGamz;
		double alp,hamn,ham;
		double Kinv;
		double Arad;
		// double dM,dR;
		double null_exp_p,null_exp_m;
		// double eneratio;
		double Mass,Comp;
		
		alp= get_bv(l,k,j, 0);
		bx =  get_bv(l,k,j, 1);
		by =  get_bv(l,k,j, 2);
		bz =  get_bv(l,k,j, 3);
		bbx=  get_bv(l,k,j, 4);
		bby=  get_bv(l,k,j, 5);
		bbz=  get_bv(l,k,j, 6);
		gxx=  get_bv(l,k,j, 7)+1.;
		gyy=  get_bv(l,k,j, 8)+1.;
		gzz=  get_bv(l,k,j, 9)+get_flat_df2z(l);
		gxy=  get_bv(l,k,j,10);
		gxz=  get_bv(l,k,j,11);
		gyz=  get_bv(l,k,j,12);
		wa=   get_bv(l,k,j,13);
		akxx= get_bv(l,k,j,14);
		akyy= get_bv(l,k,j,15);
		akzz= get_bv(l,k,j,16);
		akxy= get_bv(l,k,j,17);
		akxz= get_bv(l,k,j,18);
		akyz= get_bv(l,k,j,19);
		ek=   get_bv(l,k,j,20);
		zgx=  get_bv(l,k,j,21);
		zgy=  get_bv(l,k,j,22);
		zgz=  get_bv(l,k,j,23);
				
		hamn  = get_con(l,k,j,0);
		ham = get_con(l,k,j,1);
		nM_z=get_con(l,k,j,2);
		M_z=get_con(l,k,j,3);
		dGamz=get_con(l,k,j,4);
		
		Kinv=get_outv(l,k,j,0);
		
		null_exp_p=get_outv(l,k,j,4);
		null_exp_m=get_outv(l,k,j,6);
		// dM=get_outv(l,k,j,4);
		// dR=get_outv(l,k,j,6);
		Mass=get_outv(l,k,j,9);
		
		if(l<=lui)
		Arad=get_outv(l,k,j,5);
		else
		Arad=get_outv(lui,k,j,5);
		
		if(l==lli)
		Comp=0.;
		else
		Comp=2.*Mass/Arad;
		
		// if(l>lui)
		// eneratio=0.;
		// else
		// eneratio=(dM/dR)/(pi4*pow(Arad,2)*rhobg)-1.;
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(19) << z 						//1
		<< " "  << setw(20) << alp  				//2
		<< " "  << setw(20) << bx   				//3
		<< " "  << setw(20) << by					//4
		<< " "  << setw(20) << bz					//5
		<< " "  << setw(20) << bbx  				//6
		<< " "  << setw(20) << bby					//7
		<< " "  << setw(20) << bbz					//8
		<< " "  << setw(20) << gxx					//9
		<< " "  << setw(20) << gyy					//10
		<< " "  << setw(20) << gzz					//11
		<< " "  << setw(20) << gxy					//12
		<< " "  << setw(20) << gxz					//13
		<< " "  << setw(20) << gyz					//14
		<< " "  << setw(20) << exp(wa)/exp(psibg)	//15
		<< " "  << setw(20) << akxx					//16
		<< " "  << setw(20) << akyy					//17
		<< " "  << setw(20) << akzz					//18
		<< " "  << setw(20) << akxy					//19
		<< " "  << setw(20) << akxz					//20
		<< " "  << setw(20) << akyz					//21
		<< " "  << setw(20) << ek					//22
		<< " "  << setw(20) << zgx					//23
		<< " "  << setw(20) << zgy					//24
		<< " "  << setw(20) << zgz					//25
		<< " "  << setw(20) << hamn					//26
		<< " "  << setw(20) << ham					//27
		<< " "  << setw(20) << nM_z					//28
		<< " "  << setw(20) << M_z					//29
		<< " "  << setw(20) << dGamz				//30
		<< " "  << setw(20) << Mass					//31
		<< " "  << setw(20) << Kinv					//32
		<< " "  << setw(20) << Arad					//33
		<< " "  << setw(20) << Comp					//34
		<< " "  << setw(20) << null_exp_p			//35
		<< " "  << setw(20) << null_exp_m;			//36
																								//
		/////////////// output for geometry end //////////////////////////////////////////////////

		/////////////// output for fluid start ///////////////////////////////////////////////////
		if(fluidevo)
		{
			double Ene,px,py,pz,Den,rho,Vx,Vy,Vz,eps;
			
			Ene=  get_bv(l,k,j,24);
			px=   get_bv(l,k,j,25);
			py=   get_bv(l,k,j,26);
			pz=   get_bv(l,k,j,27);
			Den=  get_bv(l,k,j,28);
			
			rho=  get_primv(l,k,j,0);
			Vx=   get_primv(l,k,j,1);
			Vy=   get_primv(l,k,j,2);
			Vz=   get_primv(l,k,j,3);
			eps=  get_primv(l,k,j,4);
			
			fout 
			<< " "  << setw(20) << Ene		//37
			<< " "  << setw(20) << px		//38
			<< " "  << setw(20) << py		//39
			<< " "  << setw(20) << pz		//40
			<< " "  << setw(20) << Den		//41
			<< " "  << setw(20) << rho/rhobg	//42
			<< " "  << setw(20) << Vx		//43
			<< " "  << setw(20) << Vy		//44
			<< " "  << setw(20) << Vz		//45
			<< " "  << setw(20) << eps;		//46
		}											//
		/////////////// output for fluid end   ///////////////////////////////////////////////////
		
		/////////////// output for scalar start///////////////////////////////////////////////////
		if(scalarevo)										//
		{
			double phii=get_bv(l,k,j,nsc);
			double Pi=get_bv(l,k,j,nscp);
			double sEne=get_outv(l,k,j,7);
			double sp_z=get_outv(l,k,j,8);
			
			fout
			<< " "  << setw(20) << phii			//47	//37
			<< " "  << setw(20) << Pi			//48	//38
			<< " "  << setw(20) << sEne/rhobg		//49	//39
			<< " "  << setw(20) << sp_z;			//50	//40
		}											//
		/////////////// output for scalar end  ///////////////////////////////////////////////////
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}

void Fmv0::print_Kremax(ofstream& fout)
//printing curvature invariants 
{
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << setw(20) << get_t()					//1
	<< " "  << setw(20) << get_Kremax()				//2
	<< " "  << setw(20) << get_z(lkm)				//3
	<< " "  << setw(20) << get_Weylmax()				//4
	<< " "  << setw(20) << get_z(lwm)				//5
	<< endl;
	return;
}

void Fmv0::print_const(ofstream& fout)
//printing constraints
{
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << setw(20) << get_t()					//1
	<< " "  << setw(20) << get_ham()				//2
	<< " "  << setw(20) << get_hammax()				//3
	<< " "  << setw(20) << get_mom()				//4
	<< " "  << setw(20) << get_mommax()				//5
	<< " "  << setw(20) << get_z(lhm)				//6
	<< endl;
	return;
}

//output all variables
void Fmv0::print_all(ofstream& fout)
{
	double tt = get_t();
	cout << "print_all:time=" << tt 
		<< endl;
	
	fout.setf(ios_base::fixed, ios_base::floatfield);
	fout.precision(16);
	fout << "##" 
		<< "time=" << tt 
		<< " dtp=" << dtp 
		<< " dtpp=" << dtpp 
		<< endl;
	
	int k=0;
	int j=0;
	double hxx,hyy,hzz,gxy,gxz,gyz,wa;
	double akxx,akyy,akzz,akxy,akxz,akyz,ek;
	double zgx,zgy,zgz;
	double bx,by,bz,bbx,bby,bbz;
	double alp;
		
	for(int l=lli;l<=lui+tab;l++)
	{
		alp= get_bv(l,k,j, 0);
		bx =  get_bv(l,k,j, 1);
		by =  get_bv(l,k,j, 2);
		bz =  get_bv(l,k,j, 3);
		bbx=  get_bv(l,k,j, 4);
		bby=  get_bv(l,k,j, 5);
		bbz=  get_bv(l,k,j, 6);
		hxx=  get_bv(l,k,j, 7);
		hyy=  get_bv(l,k,j, 8);
		hzz=  get_bv(l,k,j, 9);
		gxy=  get_bv(l,k,j,10);
		gxz=  get_bv(l,k,j,11);
		gyz=  get_bv(l,k,j,12);
		wa=   get_bv(l,k,j,13);
		akxx= get_bv(l,k,j,14);
		akyy= get_bv(l,k,j,15);
		akzz= get_bv(l,k,j,16);
		akxy= get_bv(l,k,j,17);
		akxz= get_bv(l,k,j,18);
		akyz= get_bv(l,k,j,19);
		ek=   get_bv(l,k,j,20);
		zgx=  get_bv(l,k,j,21);
		zgy=  get_bv(l,k,j,22);
		zgz=  get_bv(l,k,j,23);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(20) << alp 	 		//1
		<< " "  << setw(20) << bx   		//2
		<< " "  << setw(20) << by			//3
		<< " "  << setw(20) << bz			//4
		<< " "  << setw(20) << bbx  		//5
		<< " "  << setw(20) << bby			//6
		<< " "  << setw(20) << bbz			//7
		<< " "  << setw(20) << hxx			//8
		<< " "  << setw(20) << hyy			//9
		<< " "  << setw(20) << hzz			//10
		<< " "  << setw(20) << gxy			//11
		<< " "  << setw(20) << gxz			//12
		<< " "  << setw(20) << gyz			//13
		<< " "  << setw(20) << wa			//14
		<< " "  << setw(20) << akxx			//15
		<< " "  << setw(20) << akyy			//16
		<< " "  << setw(20) << akzz			//17
		<< " "  << setw(20) << akxy			//18
		<< " "  << setw(20) << akxz			//19
		<< " "  << setw(20) << akyz			//20
		<< " "  << setw(20) << ek			//21
		<< " "  << setw(20) << zgx			//22
		<< " "  << setw(20) << zgy			//23
		<< " "  << setw(20) << zgz;			//24
		
		if(fluidevo)
		{
			double Ene,px,py,pz,Den;

			Ene=  get_bv(l,k,j,24);
			px=   get_bv(l,k,j,25);
			py=   get_bv(l,k,j,26);
			pz=   get_bv(l,k,j,27);
			Den=  get_bv(l,k,j,28);
			
			fout 
			<< " "  << setw(20) << Ene		//25
			<< " "  << setw(20) << px		//26
			<< " "  << setw(20) << py		//27
			<< " "  << setw(20) << pz		//28
			<< " "  << setw(20) << Den;		//29
		}	
		
		if(scalarevo)
		{
			
			double phii,Pi;

			phii=get_bv(l,k,j,nsc);
			Pi=get_bv(l,k,j,nscp);
			
			fout
			<< " "  << setw(20) << phii		//30	//25
			<< " "  << setw(20) << Pi;		//31	//26
		}
		
		fout << endl;
	}
	fout << endl << endl;
	
	for(int l=lli;l<=lui+tab;l++)
	{
		int k=0;
		int j=0;
		double hxx,hyy,hzz,gxy,gxz,gyz,wa;
		double akxx,akyy,akzz,akxy,akxz,akyz,ek;
		double zgx,zgy,zgz;
		double bx,by,bz,bbx,bby,bbz;
		double alp;
		
		alp= get_bv0(l,k,j, 0);
		bx =  get_bv0(l,k,j, 1);
		by =  get_bv0(l,k,j, 2);
		bz =  get_bv0(l,k,j, 3);
		bbx=  get_bv0(l,k,j, 4);
		bby=  get_bv0(l,k,j, 5);
		bbz=  get_bv0(l,k,j, 6);
		hxx=  get_bv0(l,k,j, 7);
		hyy=  get_bv0(l,k,j, 8);
		hzz=  get_bv0(l,k,j, 9);
		gxy=  get_bv0(l,k,j,10);
		gxz=  get_bv0(l,k,j,11);
		gyz=  get_bv0(l,k,j,12);
		wa=   get_bv0(l,k,j,13);
		akxx= get_bv0(l,k,j,14);
		akyy= get_bv0(l,k,j,15);
		akzz= get_bv0(l,k,j,16);
		akxy= get_bv0(l,k,j,17);
		akxz= get_bv0(l,k,j,18);
		akyz= get_bv0(l,k,j,19);
		ek=   get_bv0(l,k,j,20);
		zgx=  get_bv0(l,k,j,21);
		zgy=  get_bv0(l,k,j,22);
		zgz=  get_bv0(l,k,j,23);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(20) << alp 	 		//1
		<< " "  << setw(20) << bx   			//2
		<< " "  << setw(20) << by			//3
		<< " "  << setw(20) << bz			//4
		<< " "  << setw(20) << bbx  			//5
		<< " "  << setw(20) << bby			//6
		<< " "  << setw(20) << bbz			//7
		<< " "  << setw(20) << hxx			//8
		<< " "  << setw(20) << hyy			//9
		<< " "  << setw(20) << hzz			//10
		<< " "  << setw(20) << gxy			//11
		<< " "  << setw(20) << gxz			//12
		<< " "  << setw(20) << gyz			//13
		<< " "  << setw(20) << wa			//14
		<< " "  << setw(20) << akxx			//15
		<< " "  << setw(20) << akyy			//16
		<< " "  << setw(20) << akzz			//17
		<< " "  << setw(20) << akxy			//18
		<< " "  << setw(20) << akxz			//19
		<< " "  << setw(20) << akyz			//20
		<< " "  << setw(20) << ek			//21
		<< " "  << setw(20) << zgx			//22
		<< " "  << setw(20) << zgy			//23
		<< " "  << setw(20) << zgz;			//24
		
		if(fluidevo)
		{
			double Ene,px,py,pz,Den;

			Ene=  get_bv0(l,k,j,24);
			px=   get_bv0(l,k,j,25);
			py=   get_bv0(l,k,j,26);
			pz=   get_bv0(l,k,j,27);
			Den=  get_bv0(l,k,j,28);
			
			fout 
			<< " "  << setw(20) << Ene		//25
			<< " "  << setw(20) << px		//26
			<< " "  << setw(20) << py		//27
			<< " "  << setw(20) << pz		//28
			<< " "  << setw(20) << Den;		//29
		}	
		
		if(scalarevo)
		{
			
			double phii,Pi;

			phii=get_bv0(l,k,j,nsc);
			Pi=get_bv0(l,k,j,nscp);
			
			fout
			<< " "  << setw(20) << phii		//30	//25
			<< " "  << setw(20) << Pi;		//31	//26
		}
		
		fout << endl;
	}
	fout << endl << endl;
	
	for(int l=lli;l<=lui+tab;l++)
	{
		alp= get_bv1(l,k,j, 0);
		bx =  get_bv1(l,k,j, 1);
		by =  get_bv1(l,k,j, 2);
		bz =  get_bv1(l,k,j, 3);
		bbx=  get_bv1(l,k,j, 4);
		bby=  get_bv1(l,k,j, 5);
		bbz=  get_bv1(l,k,j, 6);
		hxx=  get_bv1(l,k,j, 7);
		hyy=  get_bv1(l,k,j, 8);
		hzz=  get_bv1(l,k,j, 9);
		gxy=  get_bv1(l,k,j,10);
		gxz=  get_bv1(l,k,j,11);
		gyz=  get_bv1(l,k,j,12);
		wa=   get_bv1(l,k,j,13);
		akxx= get_bv1(l,k,j,14);
		akyy= get_bv1(l,k,j,15);
		akzz= get_bv1(l,k,j,16);
		akxy= get_bv1(l,k,j,17);
		akxz= get_bv1(l,k,j,18);
		akyz= get_bv1(l,k,j,19);
		ek=   get_bv1(l,k,j,20);
		zgx=  get_bv1(l,k,j,21);
		zgy=  get_bv1(l,k,j,22);
		zgz=  get_bv1(l,k,j,23);
		
		fout.setf(ios_base::fixed, ios_base::floatfield);
		fout.precision(16);
		fout << setw(20) << alp 	 		//1
		<< " "  << setw(20) << bx   			//2
		<< " "  << setw(20) << by			//3
		<< " "  << setw(20) << bz			//4
		<< " "  << setw(20) << bbx  			//5
		<< " "  << setw(20) << bby			//6
		<< " "  << setw(20) << bbz			//7
		<< " "  << setw(20) << hxx			//8
		<< " "  << setw(20) << hyy			//9
		<< " "  << setw(20) << hzz			//10
		<< " "  << setw(20) << gxy			//11
		<< " "  << setw(20) << gxz			//12
		<< " "  << setw(20) << gyz			//13
		<< " "  << setw(20) << wa			//14
		<< " "  << setw(20) << akxx			//15
		<< " "  << setw(20) << akyy			//16
		<< " "  << setw(20) << akzz			//17
		<< " "  << setw(20) << akxy			//18
		<< " "  << setw(20) << akxz			//19
		<< " "  << setw(20) << akyz			//20
		<< " "  << setw(20) << ek			//21
		<< " "  << setw(20) << zgx			//22
		<< " "  << setw(20) << zgy			//23
		<< " "  << setw(20) << zgz;			//24
		
		if(fluidevo)
		{
			double Ene,px,py,pz,Den;

			Ene=  get_bv1(l,k,j,24);
			px=   get_bv1(l,k,j,25);
			py=   get_bv1(l,k,j,26);
			pz=   get_bv1(l,k,j,27);
			Den=  get_bv1(l,k,j,28);
			
			fout 
			<< " "  << setw(20) << Ene		//25
			<< " "  << setw(20) << px		//26
			<< " "  << setw(20) << py		//27
			<< " "  << setw(20) << pz		//28
			<< " "  << setw(20) << Den;		//29
		}	
		
		if(scalarevo)
		{
			
			double phii,Pi;

			phii=get_bv1(l,k,j,nsc);
			Pi=get_bv1(l,k,j,nscp);
			
			fout
			<< " "  << setw(20) << phii		//30	//25
			<< " "  << setw(20) << Pi;		//31	//26
		}
		
		fout << endl;
	}
	fout << endl << endl;
	return;
}



