/*********************************************************************************************************/
/*-------------------------------------------------------------------------------------------------------*/
/* BSSN EVOLUTION :: BSSN evolution Class of COSMOS_S                                                    */
/*                                             ver. 1.00          coded by Chulmoon Yoo                  */
/*-------------------------------------------------------------------------------------------------------*/
/*********************************************************************************************************/
#include "cosmos_s.h"
#include <algorithm>

void Fmv0::BSSN(int itype)
{
	//time deviation and fraction for Rungekutta start
	double frac;
	double tt;
	
	switch(itype){
	case 1:
		dt=0.5*dt0;
		frac=1./6.;
		tt=get_t();
		break;
	case 2:
		dt=0.5*dt0;
		frac=1./3.;
		tt=get_t()+dt;
		break;
	case 3:
		dt=dt0;
		frac=1./3.;
		tt=get_t()+0.5*dt;
		break;
	case 4:
		dt=dt0;
		frac=1./6.;
		tt=get_t()+dt;
		break;
	default:
		frac=1./6.;
		tt=get_t();
	}
	//time deviation and fraction for Rungekutta end
		
	// @@@@@@@@@@@@   variables (BSSN +GAUGE +fluid +scalar)   @@@@@@@@@@@
	//  0:alpha , 1:bx  ,  2:by  ,  3:bz   ,  -> Gauge
	//  4:bbx ,  5:bby  ,  6:bbz , -> for hyper-Gamma driver
	//  7:gxx   , 8:gyy ,  9:gzz , 10:gxy  , 11:gxz , 12:gyz  ,  -> tilde gamma
	//  13:wa , -> psi 
	//  14:akxx  ,15:akyy, 16:akzz, 17:akxy , 18:akxz, 19:akyz , -> tilde A_ij
	//  20:ek  , -> trK
	// 21:zgx   ,22:zgy , 23:zgz -> tilde Gamma
	// 24:E, 25:S_x, 26:S_y, 27:S_z, 28:N_B(D), 
	// 29(25):phi, 30(25) :Pi -> scalar, if no fluid, the numbers are 24 and 25
	
	if(exc)
	 excision();
	
	//set_zero_d_exc();
	set_zero_d();
	
	#pragma omp barrier

	BSSN_adv();
	//advection terms with \beta^i
	#pragma omp barrier
	
	if(fluidevo)
	{
		dyntoprim();
		#pragma omp barrier
		flux_fill();
		#pragma omp barrier
	}
	
	
	#pragma omp parallel for
	for(int l=lli;l<=lui;l++)
	{
		int k=0;
		int j=0;

		if(get_bflag(l,k,j)!=0)
		 continue;
		
		//definitions of variables start
		double alpha_p,								//lapse
		gxx_p,gyy_p,gzz_p,							//conformal metric
		wa_p,										//conformal factor(e^(4 wa_p))
		akxx_p,akzz_p,ek_p;							//extrinsic curvature
		
		double bz_p;								//shift
		
		double zgz_p;								//\tilde Gamma^j(see https://arxiv.org/abs/gr-qc/0703035v1)
		
		//x-derivatives
		double bx_x_p,
		dxz_x_p,dzx_x_p,
		zgx_x_p;

		//y-derivatives
		double by_y_p;
		
		//z-derivatives
		double alpha_z_p,bz_z_p,
		dxx_z_p,dzz_z_p,wa_z_p,
		ek_z_p,
		zgz_z_p;
		
		//second derivatives of lapse
		double alpha_xx,alpha_zz;
		
		//second derivatives of shift
		double bx_xz;
		double by_yz;
		double bz_xx,bz_yy,bz_zz;
		
		//second derivatives of conformal factor
		double wa_xx,wa_zz;
		
		//second derivatives of conformal metric
		double gxx_xx,gxx_yy,gxx_zz;
		double gyy_xx;
		double gzz_xx,gzz_yy,gzz_zz;
		
		//conformal metric determinant and inverse
		double det_p,det_pi;
		
		//inverse conformal metric
		double gixx_p,giyy_p,gizz_p;
		
		//\tilde A_i^j
		double akx_ux_p,aky_uy_p,akz_uz_p;
		
		//\tilde A^{ij}
		double ak_ixx,ak_iyy,ak_izz;
		
		//  \tilde{\Gamma}_{ijk}=\tilde{\gamma}_{im}\tilde{\Gamma}^m_{jk} -> crd{i}_{jk}
		double crdx_xz;
		double crdz_xx,crdz_zz;
		
		//  \tilde{\Gamma}^i_{jk} -> cr{i}_{jk}
		double crx_xz,cry_yz,crz_xx,crz_yy,crz_zz;

		//tilde Gamma0^i =- cal D_j tilde gamma ^ji=tilde gamma^jk Delta^i_jk		//(eq(9.80) in https://arxiv.org/abs/gr-qc/0703035v1
		double gamma0_z;
		
		// \tilde D_i D_j \psi = \psi,ij -\tilde{\Gamma}^k_ij \psi_k
		// \tilde Laplacian \psi
		// (\tilde D_i \psi)(\tilde D^i \psi)
		double wa_cdxx,wa_cdyy,wa_cdzz,wa_lap,wawa;
		
		// exp(-4\psi)* R_{jk} -> rc_{ij}; 
		double rc_xx_p,rc_yy_p,rc_zz_p;
		
		//\tilde A_i^j \tilde A_j^i
		// R (not \tilde R)
		double aaaa,ricci;
		
		//D_i beta^i and 2/3 D_i beta^i
		double divbeta,divbeta23;
		
		//time forward values
		double fwa,fgxx,fgyy,fgzz;
		
		// \tilde{\gamma}^{ij} alpha_i \phi_j -> alphawa_ip
		// exp(-4*\phi) D_i D_j \alpha
		// Lap \alpha
		double alphawa_ip,alpha_cdxx4,alpha_cdyy4,alpha_cdzz4,
			alpha_cdtr;
		
		//temporary variables for fak
		double ftmp_xx,ftmp_yy,ftmp_zz,ftmptr,trs;
		
		//time forward values for extrinsic curvature
		double fakxx,fakyy,fakzz,fek;
		
		//Laplacian for shift and derivative of divergence of shift
		double bz_lap,divb_z;
		
		//time forward values for \tilde \Gamma^i
		double fzgz;
		
		//time forward values for gauge
		double fbz,fbbz,falpha;
		
		//stress-energy tensor components
		double Ene=0.;
		double pz;
		double p_z=0.;
		double sxx=0.;
		double syy=0.;
		double szz=0.;


		//reference metric variables
		double fzz,dfzz,dfz,Z;
		
		//bar Gamma_ui_jk for inhomogeneous grid 
		double Gam_uz_zz;
		
		//double delG_1_u2_34
		double delG_z_uz_zz;
		
		//cal D_i tilde gamma_jk
		double Dgam_x_xz,Dgam_x_zx,Dgam_y_yz,Dgam_z_xx,Dgam_z_zz;
		
		//Delta_ui_jk=tilde Gamma - bar Gamma
		double Del_ux_xz,Del_ux_zx,Del_uy_yz,Del_uy_zy,Del_uz_xx,Del_uz_yy,Del_uz_zz;
		
		//cal D_i tilde gamma^jk
		double Dgam_x_uxz,Dgam_z_uxx,Dgam_z_uyy,Dgam_z_uzz;
		
		//\tilde gamma^kl \del_k \del_l \tilde \gamma_ij
		double lapgam_xx,lapgam_zz;
		
		//\tilde gamma^kl \cal D_k \cal D_l \tilde \gamma_ij
		double gDDg_xx,gDDg_zz;
		
		//\cal D_i \tilde Gamma^j
		double DGam_x_ux,DGam_z_uz;
		
		//exp(-4*wa), exp(-8*wa), exp(-12*wa)
		double ewa4i,ewa8i,ewa12i;
		
		//fij for inhomogeneous grid 
		fzz=get_flat_df2z(l);
		dfz=sqrt(fzz);
		Z=get_coordZ(l);

		//substitution of geometrical variables start
		alpha_p=get_bv(l,k,j, 0);
		bz_p=   get_bv(l,k,j, 3);
		gxx_p=  get_bv(l,k,j, 7)+1.;
		gyy_p=  get_bv(l,k,j, 8)+1.;
		gzz_p=  get_bv(l,k,j, 9)+fzz;
		wa_p=   get_bv(l,k,j,13);
		
		ewa4i=exp(-4.*wa_p);
		ewa8i=exp(-8.*wa_p);
		ewa12i=exp(-12.*wa_p);

		akxx_p= get_bv(l,k,j,14);
		akzz_p= get_bv(l,k,j,16);
		ek_p=   get_bv(l,k,j,20);
		
		zgz_p=  get_bv(l,k,j,23);
		//substitution of geometrical variables end
		
		//bar Gamma_ui_jk for inhomogeneous grid 
		Gam_uz_zz=get_flat_Gamz(l);
		
		//calculation of derivatives of flat Christoffel
		//double delG_1_u2_34
		delG_z_uz_zz=get_flat_dGamz(l);
		
		// first derivatives 4-th order
		// \del_x Gauge (\alpha \beta^i)
		bx_x_p=get_f_x(l,k,j,1);
		
		dfzz=2.*fzz*Gam_uz_zz;
		
		// \del_x \tilde{gamma}_{ij} *0.5 
		dxz_x_p=get_f_x(l,k,j,11)*0.5;
		dzx_x_p=dxz_x_p;
		
		//\del_x \Gamma^i
		zgx_x_p=get_f_x(l,k,j,21);
		
		// \del_y Gauge (\alpha \beta^i)
		by_y_p=bx_x_p;
		
		// \del_z Gauge (\alpha \beta^i)
		alpha_z_p=get_f_z(l,k,j,0);
		bz_z_p=get_f_z(l,k,j,3);
		
		// \del_z \tilde \gamma_{ij} * 0.5
		dxx_z_p=get_f_z(l,k,j,7)*0.5;
		dzz_z_p=(get_f_z(l,k,j,9)+dfzz)*0.5;

		// \del_z \psi
		wa_z_p=get_f_z(l,k,j,13);
		
		// \del_z trK
		ek_z_p=get_f_z(l,k,j,20);
		
		// \del_z \Gamma^i
		zgz_z_p=get_f_z(l,k,j,23);
		
		alpha_xx=get_f_xx(l,k,j,0);
		alpha_zz=get_f_zz(l,k,j,0);

		// \del_i\del_j \beta^x
		bx_xz=get_f_xz(l,k,j,1);
		
		// \del_i\del_j \beta^y
		by_yz=bx_xz;
		
		// \del_i\del_j \beta^z
		bz_xx=get_f_xx(l,k,j,3);
		bz_yy=bz_xx;
		bz_zz=get_f_zz(l,k,j,3);

		//\del_x \del_x g_ij
		gxx_xx=get_f_xx(l,k,j,7);
		gyy_xx=get_f_xx(l,k,j,8);
		gzz_xx=get_f_xx(l,k,j,9);
		
		//\del_y \del_y g_ij
		gxx_yy=gyy_xx;
		gzz_yy=gzz_xx;
		
		//\del_z \del_z g_ij
		gxx_zz=get_f_zz(l,k,j,7);
		gzz_zz=get_f_zz(l,k,j,9)+2.*fzz*(delG_z_uz_zz+2.*pow(Gam_uz_zz,2));
		
		//\del_i\del_j \psi
		wa_xx=get_f_xx(l,k,j,13);
		wa_zz=get_f_zz(l,k,j,13);
		
		// tilde{gamma}^{ij} -> gi
		det_p=gxx_p*gyy_p*gzz_p;
		if(det_p<1.e-16) 
		 det_p=1.e-16;
		//	det_p=get_flat_df2x(j)*get_flat_df2y(k)*get_flat_df2z(l);
		det_pi=1./det_p;
		gixx_p= (gyy_p*gzz_p)*det_pi;
		giyy_p= gixx_p;
		gizz_p= (gxx_p*gyy_p)*det_pi;
		
		//  \tilde{A}_i^k=\tilde{A}_{ij}\tilde{\gamma}^{jk} -> ak{i}_u{k}
		akx_ux_p=akxx_p*gixx_p ;
		aky_uy_p=akx_ux_p;
		akz_uz_p=akzz_p*gizz_p;

		//  \tilde{A}^{ij} -> ak_i{ij}
		ak_ixx=gixx_p*akx_ux_p ;
		ak_iyy=ak_ixx;
		ak_izz=gizz_p*akz_uz_p;

		//  \tilde{\Gamma}_{ijk}=\tilde{\gamma}_{im}\tilde{\Gamma}^m_{jk} -> crd{i}_{jk}
		// NOTE:d.._._p is 0.5*derivatives
		crdx_xz=dxx_z_p+dxz_x_p -dxz_x_p;
		
		crdz_xx=dzx_x_p+dzx_x_p -dxx_z_p;
		crdz_zz=dzz_z_p+dzz_z_p -dzz_z_p;
		
		//  \tilde{\Gamma}^i_{jk} -> cr{i}_{jk}
		crx_xz=gixx_p*crdx_xz ;
		cry_yz=crx_xz;
		
		crz_xx=gizz_p*crdz_xx;
		crz_yy=crz_xx;
		crz_zz=gizz_p*crdz_zz;
		
		//cal D_i tilde gamma_jk

		//Dgam_1_12=2.*d12_1_p-Gam_u1_11*g12_p;
		Dgam_x_xz=2.*dxz_x_p;
		Dgam_x_zx=Dgam_x_xz;
		Dgam_y_yz=Dgam_x_xz;
		
		Dgam_z_xx=2.*dxx_z_p;
		Dgam_z_zz=2.*(dzz_z_p-Gam_uz_zz*gzz_p);
		
		//Delta_ui_jk=tilde Gamma - bar Gamma
		Del_ux_xz=crx_xz;
		Del_ux_zx=crx_xz;
		Del_uy_yz=cry_yz;
		Del_uy_zy=cry_yz;
		
		Del_uz_xx=crz_xx;
		Del_uz_yy=crz_yy;
		Del_uz_zz=crz_zz-Gam_uz_zz;
		
		//tilde Gamma0^i =- cal D_j tilde gamma ^ji=tilde gamma^jk Delta^i_jk
		gamma0_z=gixx_p*Del_uz_xx+giyy_p*Del_uz_yy+gizz_p*Del_uz_zz;
		
		// \tilde D_i D_j \psi = \psi,ij -\tilde{\Gamma}^k_ij \psi_k
		wa_cdxx=wa_xx-(crz_xx*wa_z_p);
		wa_cdyy=wa_cdxx;
		wa_cdzz=wa_zz-(crz_zz*wa_z_p);
		
		// \tilde Laplacian \psi
		wa_lap=gixx_p*wa_cdxx +giyy_p*wa_cdyy +gizz_p*wa_cdzz;
		
		// (\tilde D_i \psi)(\tilde D^i \psi)
		wawa=gizz_p*wa_z_p*wa_z_p;

		//cal D_i tilde gamma^jk
		
		//Dgam_1_u23=-Del_u2_1x*gix3_p-Del_u2_1y*giy3_p-Del_u2_1z*giz3_p-Del_u3_1x*gix2_p-Del_u3_1y*giy2_p-Del_u3_1z*giz2_p;
		Dgam_x_uxz=-Del_ux_xz*gizz_p-Del_uz_xx*gixx_p;
		
		Dgam_z_uxx=-Del_ux_zx*gixx_p-Del_ux_zx*gixx_p;
		Dgam_z_uyy=Dgam_z_uxx;
		Dgam_z_uzz=-Del_uz_zz*gizz_p-Del_uz_zz*gizz_p;
		
		//\tilde gamma^kl \del_k \del_l \tilde \gamma_ij
		
		lapgam_xx=gixx_p*gxx_xx +giyy_p*gxx_yy +gizz_p*gxx_zz;
		lapgam_zz=gixx_p*gzz_xx +giyy_p*gzz_yy +gizz_p*gzz_zz;
		
		
		//\tilde gamma^kl \cal D_k \cal D_l \tilde \gamma_ij
				
		//gDDg_11=lapgamma_11-2.*gi11_p*delG_1_u1_11*g11_p
		//			-8.*Gam_u1_11*(gix1_p*d11_x_p+giy1_p*d11_y_p+giz1_p*d11_z_p)+4.*gi11_p*pow(Gam_u1_11,2)*g11_p
		//			-(gizz_p*Gam_uz_zz*Dgam_z_11);
		gDDg_xx=lapgam_xx-(gizz_p*Gam_uz_zz*Dgam_z_xx);
		gDDg_zz=lapgam_zz-2.*gizz_p*delG_z_uz_zz*gzz_p
					-8.*Gam_uz_zz*(gizz_p*dzz_z_p)+4.*gizz_p*pow(Gam_uz_zz,2)*gzz_p
					-(gizz_p*Gam_uz_zz*Dgam_z_zz);
		
		//gDDg_12=lapgam_12-gi11_p*delG_1_u1_11*g12_p-gi22_p*delG_2_u2_22*g12_p
		//			-4.*Gam_u1_11*(gix1_p*d12_x_p+giy1_p*d12_y_p+giz1_p*d12_z_p)-4.*Gam_u2_22*(gix2_p*d12_x_p+giy2_p*d12_y_p+giz2_p*d12_z_p)
		//			+(gi11_p*pow(Gam_u1_11,2)+2.*gi12_p*Gam_u1_11*Gam_u2_22+gi22_p*pow(Gam_u2_22,2))*g12_p
		//			-(gizz_p*Gam_uz_zz*Dgam_z_12);
		
		//\cal D_i \tilde Gamma^j
		
		DGam_x_ux=zgx_x_p;
		DGam_z_uz=zgz_z_p+Gam_uz_zz*zgz_p;
		
		// exp(-4\psi)* R_{jk} -> rc_{ij};  !!! NOT R_{jk}
		// see other files for previous version
		//rc_12_p=0.5*(-gDDg_12
		//			+g1x_p*DGam_2_ux+g1y_p*DGam_2_uy+g1z_p*DGam_2_uz
		//			+g2x_p*DGam_1_ux+g2y_p*DGam_1_uy+g2z_p*DGam_1_uz)
		//		-0.5*(Dgam_x_x1*Dgam_2_uxx+Dgam_y_y1*Dgam_2_uyy+Dgam_z_z1*Dgam_2_uzz
		//			+Dgam_x_y1*Dgam_2_uxy+Dgam_x_z1*Dgam_2_uxz+Dgam_y_z1*Dgam_2_uyz
		//			+Dgam_y_x1*Dgam_2_uxy+Dgam_z_x1*Dgam_2_uxz+Dgam_z_y1*Dgam_2_uyz
		//			
		//			+Dgam_x_x2*Dgam_1_uxx+Dgam_y_y2*Dgam_1_uyy+Dgam_z_z2*Dgam_1_uzz
		//			+Dgam_x_y2*Dgam_1_uxy+Dgam_x_z2*Dgam_1_uxz+Dgam_y_z2*Dgam_1_uyz
		//			+Dgam_y_x2*Dgam_1_uxy+Dgam_z_x2*Dgam_1_uxz+Dgam_z_y2*Dgam_1_uyz
		//			
		//			-zgx_p*Dgam_x_12-zgy_p*Dgam_y_12-zgz_p*Dgam_z_12)
		//		
		//		-Del_ux_1x*Del_ux_x2-Del_uy_1y*Del_uy_y2-Del_uz_1z*Del_uz_z2
		//		-Del_ux_1y*Del_uy_x2-Del_ux_1z*Del_uz_x2-Del_uy_1z*Del_uz_y2
		//		-Del_uy_1x*Del_ux_y2-Del_uz_1x*Del_ux_z2-Del_uz_1y*Del_uy_z2;
		
		rc_xx_p=(0.5*(-gDDg_xx
					+gxx_p*DGam_x_ux
					+gxx_p*DGam_x_ux)
				-0.5*(Dgam_x_zx*Dgam_x_uxz
					+Dgam_z_xx*Dgam_x_uxz
					
					+Dgam_x_zx*Dgam_x_uxz
					+Dgam_z_xx*Dgam_x_uxz
					
					-zgz_p*Dgam_z_xx)
					//-gamma0_z*Dgam_z_xx)
				
				-Del_ux_xz*Del_uz_xx
				-Del_uz_xx*Del_ux_zx
				
				-2.*wa_cdxx -2.*wa_lap*gxx_p -4.*wawa*gxx_p )*ewa4i;
		rc_yy_p=rc_xx_p;

		rc_zz_p=(0.5*(-gDDg_zz+gzz_p*DGam_z_uz+gzz_p*DGam_z_uz)
				-0.5*(Dgam_x_xz*Dgam_z_uxx+Dgam_y_yz*Dgam_z_uyy+Dgam_z_zz*Dgam_z_uzz
					
					+Dgam_x_xz*Dgam_z_uxx+Dgam_y_yz*Dgam_z_uyy+Dgam_z_zz*Dgam_z_uzz
					
					-zgz_p*Dgam_z_zz)
					//-gamma0_z*Dgam_z_zz)
				
				-Del_ux_zx*Del_ux_xz-Del_uy_zy*Del_uy_yz-Del_uz_zz*Del_uz_zz
				
				-2.*wa_cdzz -2.*wa_lap*gzz_p -4.*wawa*gzz_p +4.*wa_z_p*wa_z_p)*ewa4i;
		
		// R 
		ricci=gixx_p*rc_xx_p +giyy_p*rc_yy_p +gizz_p*rc_zz_p;
		
		//\tilde A_i^j \tilde A_j^i
		aaaa=akx_ux_p*akx_ux_p 
			+aky_uy_p*aky_uy_p 
			+akz_uz_p*akz_uz_p;
		
						
		//D_i beta^i
		divbeta=bx_x_p +by_y_p +bz_z_p+Gam_uz_zz*bz_p;
		//divbeta=bx_x_p +by_y_p +bz_z_p;
		divbeta23=2./3.*divbeta;

		////////////// for scalar field ///////////////////
		if(scalarevo)
		{
			//time forward values for scalar field and conjugate momentum
			double fphi,fPi;
			
			//scalar field and conjugate momentum, those derivatives
			double phi,Pi,phi_z,phi_xx,phi_zz;
			
			// \tilde D_i D_j \phi = \phi,ij -\tilde{\Gamma}^k_ij \phi_k
			double phi_cdxx,phi_cdyy,phi_cdzz;
			
			// \tilde Laplacian \psi
			// (\tilde D_i \phi)(\tilde D^i \phi)
			// (\tilde D_i \psi)(\tilde D^i \phi)
			// (\tilde D_i \alpha)(\tilde D^i \phi)
			double phi_lap,dphidphi,dpsidphi,dalphadphi;
			
			//scalar field potential and derivative
			double Vpot,dVpot;
			
			//substitution start
			phi=get_bv(l,k,j,nsc);
			Pi=get_bv(l,k,j,nscp);

			phi_z=get_f_z(l,k,j,nsc);
			phi_xx=get_f_xx(l,k,j,nsc);
			phi_zz=get_f_zz(l,k,j,nsc);
	
			// \tilde D_i D_j \phi = \phi,ij -\tilde{\Gamma}^k_ij \phi_k
			// ok even for non-Cartesian
			phi_cdxx=phi_xx-(crz_xx*phi_z);
			phi_cdyy=phi_cdxx;
			phi_cdzz=phi_zz-(crz_zz*phi_z);
			
			// \tilde Laplacian \psi
			phi_lap=gixx_p*phi_cdxx +giyy_p*phi_cdyy +gizz_p*phi_cdzz;
			
			// (\tilde D_i \phi)(\tilde D^i \phi)
			dphidphi=gizz_p*phi_z*phi_z;
			
			// (\tilde D_i \psi)(\tilde D^i \phi)
			dpsidphi=gizz_p*wa_z_p*phi_z;
			
			// (\tilde D_i \palpha)(\tilde D^i \phi)
			dalphadphi=gizz_p*alpha_z_p*phi_z;
			
			Vpot=funcV(phi);
			dVpot=funcdV(phi);

			///////////////////////////////
			// Scalar field evolution 
			///////////////////////////////
			
			fphi=-alpha_p*Pi;
			fPi=-ewa4i*(alpha_p*(2.*dpsidphi+phi_lap)+dalphadphi)+alpha_p*(ek_p*Pi+dVpot);
			
			set_dbv(l,k,j,nsc)=get_dbv(l,k,j,nsc)+fphi;
			set_dbv(l,k,j,nscp)=get_dbv(l,k,j,nscp)+fPi;
			
			//////////////////////////////////////////////////////
			// scalar field contribution to stress energy tensor
			//////////////////////////////////////////////////////
			
			Ene=0.5*(pow(Pi,2)+ewa4i*dphidphi)+Vpot;

			p_z=Pi*phi_z;
			
			if(itype==1)
			{
				set_outv(l,k,j,7)=Ene;
				set_outv(l,k,j,8)=p_z;
			}
			
			sxx=((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gxx_p;
			syy=sxx;
			szz=phi_z*phi_z+((0.5*Pi*Pi-Vpot)/ewa4i-0.5*dphidphi)*gzz_p;
		}
		
		if(fluidevo)
		{
			//fluid stress-energy tensor components
			double fEne,fpz,fp_z,fsxx,fsyy,fszz;

			//E+P and P 
			double EpP,press;

			//\sqrt(\gamma)
			double sqgam;
			
			//time forward values for fluid variables
			double fS0,fSz,fDen;
			
			//extrinsic curvature componens
			double exKuu_xx,exKuu_yy,exKuu_zz;
			
			//S_ij K^ij
			double SK;
			
			//fluid dynamical variables
			double S0,Sz;
			
			//stress tensor componets
			double fsuxx,fsuyy,fsuzz;
			
			//\del_i \gamma^jk;
			double dz_gamxx,dz_gamyy,dz_gamzz;

			//substitution start
			sqgam=sqrt(det_p/ewa12i);
			
			//energy momentum tensor of fluid
			fEne=get_bv(l,k,j,24)/sqgam;
			fp_z=get_bv(l,k,j,27)/sqgam;
			
			if(itype==1)
			 set_outv(l,k,j,10)=Ene/fEne;
			
			//energy momentum tensor
			Ene+=fEne;
			p_z+=fp_z;
		
			fpz=(gizz_p*fp_z)*ewa4i;
		
			press=pres(get_primv(l,k,j,0));
			
			EpP=fEne+press;
			
			press/=ewa4i;
			fsxx=press*gxx_p;
			fsyy=fsxx;
			fszz=fp_z*fp_z/EpP+press*gzz_p;
			
			//stress tensor components
			sxx+=fsxx;
			syy+=fsyy;
			szz+=fszz;
			
			///////////////////////////////
			// fluid evolution
			///////////////////////////////
				
			fS0=-dxi*(get_flux_x(l,k,j+1,0)-get_flux_x(l,k,j,0));
			fSz=-dxi*(get_flux_x(l,k,j+1,1)-get_flux_x(l,k,j,1));
			fDen=-dxi*(get_flux_x(l,k,j+1,2)-get_flux_x(l,k,j,2));
			
			fS0=fS0-dyi*(get_flux_y(l,k+1,j,0)-get_flux_y(l,k,j,0));
			fSz=fSz-dyi*(get_flux_y(l,k+1,j,1)-get_flux_y(l,k,j,1));
			fDen=fDen-dyi*(get_flux_y(l,k+1,j,2)-get_flux_y(l,k,j,2));
			
			fS0=fS0-dzi*(get_flux_z(l+1,k,j,0)-get_flux_z(l,k,j,0));
			fSz=fSz-dzi*(get_flux_z(l+1,k,j,1)-get_flux_z(l,k,j,1));
			fDen=fDen-dzi*(get_flux_z(l+1,k,j,2)-get_flux_z(l,k,j,2));
			
			//extrinsic curvature
			exKuu_xx=(ak_ixx+gixx_p*ek_p/3.)*ewa4i;
			exKuu_yy=exKuu_xx;
			exKuu_zz=(ak_izz+gizz_p*ek_p/3.)*ewa4i;
			
			SK=fsxx*exKuu_xx+fsyy*exKuu_yy+fszz*exKuu_zz;
			
			S0=get_bv(l,k,j,24);
			Sz=get_bv(l,k,j,27);
			
			press*=ewa8i;
			fsuxx=press*gixx_p;
			fsuyy=fsuxx;
			fsuzz=fpz*fpz/EpP+press*gizz_p;
			
			dz_gamxx=(4.*gxx_p*wa_z_p+2.*dxx_z_p)/ewa4i;
			dz_gamyy=dz_gamxx;
			dz_gamzz=(4.*gzz_p*wa_z_p+2.*dzz_z_p)/ewa4i;
			
			fS0=fS0+sqgam*(-fpz*alpha_z_p+alpha_p*SK);
									
			fSz=fSz-S0*alpha_z_p+(Sz*bz_z_p)
					+0.5*alpha_p*sqgam*(fsuxx*dz_gamxx+fsuyy*dz_gamyy+fsuzz*dz_gamzz);
			set_dbv(l,k,j,24)=fS0;
			set_dbv(l,k,j,25)=0.;
			set_dbv(l,k,j,26)=0.;
			set_dbv(l,k,j,27)=fSz;
			set_dbv(l,k,j,28)=fDen;
		}
		
		pz=ewa4i*(gizz_p*p_z);

		///////////////////////////
		// Wa,\tilde{\gamma}_{ij}
		///////////////////////////
		
		fwa=-1./6.*(alpha_p*ek_p -divbeta);
		set_dbv(l,k,j,13)=get_dbv(l,k,j,13) +fwa;
		
		fgxx=-2.*alpha_p*akxx_p
			+gxx_p*bx_x_p  
			+gxx_p*bx_x_p  
			-gxx_p*divbeta23;
		fgyy=fgxx;
		
		fgzz=-2.*alpha_p*akzz_p 
			+gzz_p*bz_z_p 
			+gzz_p*bz_z_p 
			-gzz_p*divbeta23
			+bz_p*dfzz;
		
		set_dbv(l,k,j, 7) =get_dbv(l,k,j, 7) +fgxx;
		set_dbv(l,k,j, 8) =get_dbv(l,k,j, 8) +fgyy;
		set_dbv(l,k,j, 9) =get_dbv(l,k,j, 9) +fgzz;
		set_dbv(l,k,j,10) =0.;
		set_dbv(l,k,j,11) =0.;
		set_dbv(l,k,j,12) =0.;
		
		//////////////////////////////////////
		// \tilde{\gamma}^{ij} alpha_i \phi_j
		//////////////////////////////////////
		alphawa_ip=gizz_p*alpha_z_p*wa_z_p;

		// exp(-4*\phi) D_i D_j \alpha
		alpha_cdxx4=( alpha_xx-(crz_xx*alpha_z_p)-2.*( -gxx_p*alphawa_ip) )*ewa4i;
		alpha_cdyy4=alpha_cdxx4;
		alpha_cdzz4=(alpha_zz-( crz_zz*alpha_z_p)-2.*( alpha_z_p*wa_z_p+alpha_z_p*wa_z_p 
					-gzz_p*alphawa_ip))*ewa4i;
		
		ftmp_xx=-alpha_cdxx4 +alpha_p*(rc_xx_p-pi8*sxx*ewa4i);
		ftmp_yy=ftmp_xx;
		ftmp_zz=-alpha_cdzz4 +alpha_p*(rc_zz_p-pi8*szz*ewa4i);
		
		ftmptr=gixx_p*ftmp_xx +giyy_p*ftmp_yy +gizz_p*ftmp_zz;


		fakxx=ftmp_xx -1./3.*ftmptr*gxx_p
				+alpha_p*(-2.*( akx_ux_p*akxx_p ) +ek_p*akxx_p )
				+akxx_p*bx_x_p 
				+akxx_p*bx_x_p 
				-akxx_p*divbeta23;
		fakyy=fakxx;
		
		fakzz=ftmp_zz-1./3.*ftmptr*gzz_p 
				+alpha_p*(-2.*( akz_uz_p*akzz_p) +ek_p*akzz_p ) 
				+akzz_p*bz_z_p 
				+akzz_p*bz_z_p 
				-akzz_p*divbeta23;
		set_dbv(l,k,j,14)=get_dbv(l,k,j,14) +fakxx;
		set_dbv(l,k,j,15)=get_dbv(l,k,j,15) +fakyy;
		set_dbv(l,k,j,16)=get_dbv(l,k,j,16) +fakzz;
		set_dbv(l,k,j,17)=0.;
		set_dbv(l,k,j,18)=0.;
		set_dbv(l,k,j,19)=0.;

		// \Delta \alpha
		alpha_cdtr=gixx_p*alpha_cdxx4
					+giyy_p*alpha_cdyy4
					+gizz_p*alpha_cdzz4;

		trs=(gixx_p*sxx +giyy_p*syy +gizz_p*szz)*ewa4i;
		
		fek=-alpha_cdtr +alpha_p*(aaaa + pow(ek_p,2)/3.+pi4*(Ene+trs));
		set_dbv(l,k,j,20)=get_dbv(l,k,j,20) +fek;
		
		//added gi11_p*delG_1_u1_11*b1_p
		//				+2.*(Gam_u1_11*gi1x_p*b1_x_p+Gam_u1_11*gi1y_p*b1_y_p+Gam_u1_11*gi1z_p*b1_z_p)
		//				-Gam_uz_zz*gizz_p*b1_z_p;
		bz_lap=gixx_p*bz_xx +giyy_p*bz_yy +gizz_p*bz_zz
				
				+gizz_p*delG_z_uz_zz*bz_p
				+2.*(Gam_uz_zz*gizz_p*bz_z_p)
				-Gam_uz_zz*gizz_p*bz_z_p;
		
		//modified by chulmoon for non-Cartesian inhomogeneous grid
		//added delG_1_u1_11*b1_p+Gam_uz_zz*bz_1_p
		divb_z=bx_xz +by_yz +bz_zz
				+delG_z_uz_zz*bz_p+Gam_uz_zz*bz_z_p;
		
		
		//cr._.. -> Del_u._..
		
		fzgz=-2.*(ak_izz*alpha_z_p ) 
			+2.*alpha_p*( Del_uz_xx*ak_ixx  
			+Del_uz_yy*ak_iyy 
			+Del_uz_zz*ak_izz 
			-2./3.*(gizz_p*ek_z_p)
			+6.*( ak_izz*wa_z_p)    
			-pi8*pz/ewa4i   )
			//-(bz_z_p*gamma0_z) 
			//+divbeta23*gamma0_z 
			-(bz_z_p*zgz_p) 
			+divbeta23*zgz_p 
			+1./3.*(gizz_p*divb_z )   + bz_lap;
		
		set_dbv(l,k,j,21)=0.;
		set_dbv(l,k,j,22)=0.;
		set_dbv(l,k,j,23)=get_dbv(l,k,j,23) +fzgz;
		

		///////////////////////////////
		// Gauge
		///////////////////////////////
		
		//modified 1+log
		//falpha=-etaa*(ek_p-get_bv(lui,0,0,20))*alpha_p;
		falpha=-etaa*(ek_p+2./(1.+fluidw)/tt)*alpha_p;
		
		//synchronous
		//falpha=0.;
			
		//modified harmonic
		//falpha=-(ek_p-get_bv(lui,kui,jui,20))*pow(alpha_p,2);
		
		set_dbv(l,k,j,0)=get_dbv(l,k,j,0) +falpha;

		//Gamma driver
		fbz=etabb*get_bv(l,k,j,6);

		fbbz=fzgz -etab*get_bv(l,k,j,6);
		
		//zero shift gauge
		//fbz=0.;
		//fbbz=0.;
		
		set_dbv(l,k,j,1)=0.;
		set_dbv(l,k,j,2)=0.;
		set_dbv(l,k,j,3)=get_dbv(l,k,j,3) +fbz;
		//set_dbv(l,k,j,3)=fbz;

		set_dbv(l,k,j,4)=0.;
		set_dbv(l,k,j,5)=0.;
		set_dbv(l,k,j,6)=get_dbv(l,k,j,6) +fbbz;
		//set_dbv(l,k,j,6)=fbbz;
				
		if(itype==1)
		{
			////////////////////////////////////// constraints start /////////////////////////////////////
			double hamc,nham;			//hamiltonian constraint, normalized hamiltonian constraint		//

			//\tilde A_ij,k=\del_k \tilde A_ij for constraint and curvature invariants
			double daxz_x,dayz_y,daxx_z,dazz_z;
			
			//\cal D_i \tilde A_jk for constraints
			double Daxz_x;
			double Dayz_y;
			double Dazz_z;
			
			//\tilde D^1 \tilde A_x1 for constraints
			double Da_z;
			
			//momentum constraints and normalized ones
			double M_z,normM_z;
			
			//\tilde A_ij,k
			daxz_x=get_f_x(l,k,j,18);
			dayz_y=daxz_x;
			
			daxx_z=get_f_z(l,k,j,14);
			dazz_z=get_f_z(l,k,j,16);
			
			//D_3 \tilde A_12=Da12_3
			Daxz_x=daxz_x;
			Dayz_y=dayz_y;
			Dazz_z=dazz_z-2.*Gam_uz_zz*akzz_p;
			
			//\tilde D^1 \tilde A_x1=Da_x
			Da_z=gixx_p*Daxz_x+giyy_p*Dayz_y
				+gizz_p*Dazz_z
			
				-Del_ux_xz*akx_ux_p
				-Del_uy_yz*aky_uy_p
				-Del_uz_zz*akz_uz_p
				
				//-gamma0_z*akzz_p;
				-zgz_p*akzz_p;

			M_z=Da_z+6.*(akz_uz_p*wa_z_p)-2./3.*ek_z_p-pi8*p_z;
			
			normM_z=abs(Da_z)+6.*(abs(akz_uz_p*wa_z_p))+2./3.*abs(ek_z_p)+pi8*abs(p_z);

			//hamc = ricci + 2.*pow(ek_p,2)/3. - aaaa - 2.*lambda - pi16*get_enemom(l,k,j,0);
			hamc = ricci + 2.*pow(ek_p,2)/3. - aaaa - 2.*lambda-pi16*Ene;
			
			double ricci_psi=-8.*(wawa+wa_lap)*ewa4i;

			nham=2.*pow(ek_p,2)/3.+abs(ricci-ricci_psi)+(abs(8.*wawa)+abs(8.*wa_lap))*ewa4i+abs(aaaa)+pi16*abs(Ene);
			
			set_con(l,k,j,0)=hamc/nham;
			
			set_con(l,k,j,1)=hamc;
			set_con(l,k,j,2)=M_z/(normM_z+1.);
			set_con(l,k,j,3)=M_z;
			set_con(l,k,j,4)=gamma0_z-zgz_p;
			//set_con(l,k,j,4)=Del_uz_zz;
																											//
			////////////////////////////////////// constraints end ///////////////////////////////////////////
			
			/////////////////////////////// Kodama mass and Area rad start ///////////////////////////////////
			double gxx_z=get_f_z(l,k,j,7);																	//
			// double comfac=dfz/alpha_p*sqrt(gxx_p/gzz_p);
			// double Kt=-comfac*((2.*wa_z_p+0.5*gxx_z/gxx_p)*Z/dfz+1.);
			// double Kr=comfac*Z*(2.*get_dbv(l,k,j,13)+0.5*get_dbv(l,k,j,7)/gxx_p);
			// double dM=pi4*pow(Z,2)*exp(6.*wa_p)*gxx_p*sqrt(gzz_p)/dfz*(dfz*(-Ene*alpha_p+bz_p*p_z)*Kt+p_z*Kr);
			double gamZZ=gzz_p/(ewa4i*dfz*dfz);
			double gtt=-pow(alpha_p,-2);
			double gtr=dfz*bz_p*pow(alpha_p,-2);
			double grr=-pow(dfz*bz_p/alpha_p,2)+1./gamZZ;
			
			double dtR=exp(2.*wa_p)*Z*sqrt(gxx_p)*(2.*get_dbv(l,k,j,13)+0.5*get_dbv(l,k,j,7)/gxx_p);
			double drR=exp(2.*wa_p)*sqrt(gxx_p)*((2.*wa_z_p+0.5*gxx_z/gxx_p)*Z/dfz+1.);
			double Rad=Z*exp(2.*wa_p)*sqrt(gxx_p);
			double M=0.5*Rad*(1.-gtt*pow(dtR,2)-grr*pow(drR,2)-2.*gtr*dtR*drR);
			// double dR=Rad*(dfz/Z+2.*wa_z_p+0.5*gxx_z/gxx_p);
			
			double kt=exp(2.*wa_p)*sqrt(gzz_p);
			double kzp=(-kt*bz_p+alpha_p)*dfz;
			double kzm=(-kt*bz_p-alpha_p)*dfz;

			// set_outv(l,k,j,4)=dM;
			set_outv(l,k,j,4)=kt*dtR+kzp*drR;
			set_outv(l,k,j,5)=Rad;
			// set_outv(l,k,j,6)=dR;
			set_outv(l,k,j,6)=kt*dtR+kzm*drR;
			set_outv(l,k,j,9)=M;																			//
			/////////////////////////////// Kodama mass and Area rad end    //////////////////////////////////
			
			if(curveval)
			{
				
				//////////////////////////////////////////////
				//calculation of Kretschmann invariant start//
				//////////////////////////////////////////////
				
				double Kinv;
				//gamma gamma gamma gamma R4
				double tanRiem_xyxy,tanRiem_xzxz,tanRiem_yzyz;
				//upper indeces
				double utanRiem_xyxy,utanRiem_xzxz,utanRiem_yzyz;
				//gamma gamma gamma n R4
				double nRiem_xz_x,nRiem_yz_y;
				//upper indeces
				double unRiem_xz_x,unRiem_yz_y;
				//gamma n gamma n gamma R4
				double nnRiem_xx,nnRiem_yy,nnRiem_zz;
				//upper indeces
				double unnRiem_xx,unnRiem_yy,unnRiem_zz;
				//R3
				double Riem3_xyxy,Riem3_xzxz;
				
				//upper indeces
				double uRiem3_xyxy,uRiem3_xzxz;
				
				//K_ij 
				double exK_xx,exK_yy,exK_zz;
				double exK_xz_x,exK_zx_x,
						exK_xx_z;
				
				//Gamma^i_jk
				double Gamma_x_xz,Gamma_z_xx;
				
				//Gamma^l_ij K_lk
				double GammaK_xz_x,GammaK_zx_x,GammaK_xx_z;
				
				//D^i wa
				double wa_uz;
				//K^i_j
				double exKu_xx,exKu_zz;
				
				//Ki^k K_kj
				double KK_xx,KK_zz;
				
				//R3^ij
				double uR_xx,uR_yy,uR_zz;
				//K^ij
				double uK_xx,uK_yy,uK_zz;
				//upper RC_ij
				double uRC_xx,uRC_yy,uRC_zz;
				//RC_ij
				double RC_xx,RC_yy,RC_zz;
				//RC
				double RC;
				//RC_i
				double RC_z;
				double RijRij;
				double RR;
				
				//extrinsic curvature
				exK_xx=(akxx_p+gxx_p*ek_p/3.)/ewa4i;
				exK_yy=exK_xx;
				exK_zz=(akzz_p+gzz_p*ek_p/3.)/ewa4i;
				
				//3-dim Riemann tensor Riem3_1234=exp(8.*wa_p)*(g13_p*rc_42_p-g14_p*rc_32_p+g24_p*rc_31_p-g23_p*rc_41_p
				//								-0.5*ricci*(g13_p*g42_p-g14_p*g32_p))
				Riem3_xyxy=(gxx_p*rc_yy_p+gyy_p*rc_xx_p-0.5*ricci*(gxx_p*gyy_p))/ewa8i;
				Riem3_xzxz=(gxx_p*rc_zz_p+gzz_p*rc_xx_p-0.5*ricci*(gxx_p*gzz_p))/ewa8i;
				
				//K_ij,k=exK_12_3=4.*wa_3_p*exK_12+exp(4.*wa_p)*(da12_3+d12_3_p*ek_p/3.+g12_p*ek_3_p/3.)
				exK_xz_x=(daxz_x+2.*dxz_x_p*ek_p/3.)/ewa4i;
				exK_zx_x=exK_xz_x;
				
				exK_xx_z=4.*wa_z_p*exK_xx+(daxx_z+2.*dxx_z_p*ek_p/3.+gxx_p*ek_z_p/3.)/ewa4i;
				
				//wa_u1=gi1z_p*wa_z_p
				wa_uz=gizz_p*wa_z_p;
				
				//Gamma_1_23=cr1_23+2.*(delta13*wa_1_p+delta12*wa_2_p-g23_p*wa_u1)
				Gamma_x_xz=crx_xz+2.*(wa_z_p);
				
				Gamma_z_xx=crz_xx+2.*(-gxx_p*wa_uz);
				
				//GammaK_12_3=Gamma_x_12*exK_x3+Gamma_y_12*exK_y3+Gamma_z_12*exK_z3
				GammaK_xz_x=Gamma_x_xz*exK_xx;
				GammaK_zx_x=GammaK_xz_x;
		
				//GammaK_12_3=Gamma_x_12*exK_x3+Gamma_y_12*exK_y3+Gamma_z_12*exK_z3;
				GammaK_xx_z=Gamma_z_xx*exK_zz;
		
				//tanRiem_1234=Riem3_1234+exK_13*exK_24-exK_23*exK_14;
				tanRiem_xyxy=Riem3_xyxy+exK_xx*exK_yy;
				tanRiem_xzxz=Riem3_xzxz+exK_xx*exK_zz;
				tanRiem_yzyz=tanRiem_xzxz;
			
				//nRiem_12_3=-exK_23_1+exK_13_2+GammaK_13_2-GammaK_23_1;
				nRiem_xz_x=-exK_zx_x+exK_xx_z+GammaK_xx_z-GammaK_zx_x;
				nRiem_yz_y=nRiem_xz_x;
				
				//exKu_12=gi1x_p*exK_x2+gi1y_p*exKy2+gi1z_p*exKz2;
				exKu_xx=(gixx_p*exK_xx)*ewa4i;
				exKu_zz=(gizz_p*exK_zz)*ewa4i;
				
				//KK_12=exK_1x*exKu_x2+exK_1y*exKu_y2+exK_1z*exKu_z2;
				KK_xx=exK_xx*exKu_xx;
				KK_zz=exK_zz*exKu_zz;
				
				//nnRiem_12=-KK_12+exp(4.*wa_p)*rc_12_p+ek_p*exK_12+pi4*((trs-Ene)*g12_p-2.*s12);
				nnRiem_xx=-KK_xx+rc_xx_p/ewa4i+ek_p*exK_xx+pi4*((trs-Ene)*gxx_p-2.*sxx);
				nnRiem_yy=nnRiem_xx;
				nnRiem_zz=-KK_zz+rc_zz_p/ewa4i+ek_p*exK_zz+pi4*((trs-Ene)*gzz_p-2.*szz);
				
				//unnRiem_12=exp(-8.*wa_p)*(
				//			gi1x_p*gi2x_p*nnRiem_xx+gi1y_p*gi2y_p*nnRiem_yy+gi1z_p*gi2z_p*nnRiem_zz);
				unnRiem_xx=ewa8i*(gixx_p*gixx_p*nnRiem_xx);
				unnRiem_yy=unnRiem_xx;
				unnRiem_zz=ewa8i*(gizz_p*gizz_p*nnRiem_zz);
				
				//unRiem_12_3=exp(-12.*wa_p)*(
				//			+gi1x_p*gi2z_p*gi3x_p*nRiem_xz_x
				//			-gi1z_p*gi2x_p*gi3x_p*nRiem_xz_x
				//			+gi1y_p*gi2z_p*gi3y_p*nRiem_yz_y
				//			-gi1z_p*gi2y_p*gi3y_p*nRiem_yz_y);
				//unRiem_xy_x=0.;
				unRiem_xz_x=ewa12i*(gixx_p*gizz_p*gixx_p*nRiem_xz_x);
				unRiem_yz_y=unRiem_xz_x;
				
				uR_xx=ewa4i*(gixx_p*gixx_p*rc_xx_p);
				uR_yy=uR_xx;
				uR_zz=ewa4i*(gizz_p*gizz_p*rc_zz_p);

				uK_xx=ewa8i*(gixx_p*gixx_p*exK_xx);
				uK_yy=uK_xx;
				uK_zz=ewa8i*(gizz_p*gizz_p*exK_zz);
				
				//uRiem3_xyxy=exp(-4.*wa_p)*(gi13_p*uR_42-gi14_p*uR_32+gi24_p*uR_31-gi23_p*uR_41)
				//								-0.5*ricci*exp(-8.*wa_p)*(gi13_p*gi42_p-gi14_p*gi32_p);
				uRiem3_xyxy=ewa4i*(gixx_p*uR_yy+giyy_p*uR_xx)-0.5*ricci*ewa8i*(gixx_p*giyy_p);
				uRiem3_xzxz=ewa4i*(gixx_p*uR_zz+gizz_p*uR_xx)-0.5*ricci*ewa8i*(gixx_p*gizz_p);
				
				//utanRiem_1234=uRiem3_1234+uK_13*uK_24-uK_23*uK_14;
				utanRiem_xyxy=uRiem3_xyxy+uK_xx*uK_yy;
				utanRiem_xzxz=uRiem3_xzxz+uK_xx*uK_zz;
				utanRiem_yzyz=utanRiem_xzxz;
				
				Kinv=4.*(unnRiem_xx*nnRiem_xx+unnRiem_yy*nnRiem_yy+unnRiem_zz*nnRiem_zz)
					-4.*(2.*(unRiem_xz_x*nRiem_xz_x
							+unRiem_yz_y*nRiem_yz_y))
					+4.*(utanRiem_xyxy*tanRiem_xyxy+utanRiem_xzxz*tanRiem_xzxz+utanRiem_yzyz*tanRiem_yzyz);
				
				set_outv(l,k,j,0)=Kinv;
				
				//upper RC_ij
				//uRC_12=exp(4.*wa_p)*(gxx_p*utanRiem_x1x2+gyy_p*utanRiem_y1y2+gzz_p*utanRiem_z1z2
				//						+gxy_p*utanRiem_x1y2+gyz_p*utanRiem_y1z2+gzx_p*utanRiem_z1x2
				//						+gxy_p*utanRiem_y1x2+gyz_p*utanRiem_z1y2+gzx_p*utanRiem_x1z2);
				
				uRC_xx=(gyy_p*utanRiem_xyxy+gzz_p*utanRiem_xzxz)/ewa4i;
				uRC_yy=uRC_xx;
				uRC_zz=(gxx_p*utanRiem_xzxz+gyy_p*utanRiem_yzyz)/ewa4i;
				
				//RC_ij
				
				RC_xx=(gxx_p*gxx_p*uRC_xx)/ewa8i;
				RC_yy=RC_xx;
				RC_zz=(gzz_p*gzz_p*uRC_zz)/ewa8i;

				//RC
				
				RC=(gxx_p*unnRiem_xx+gyy_p*unnRiem_yy+gzz_p*unnRiem_zz)/ewa4i;
				
				//RC_i
				
				//RC_1=exp(-4.*wa_p)*(gixx_p*nRiem_x1_x+giyy_p*nRiem_y1_y+gizz_p*nRiem_z1_z);
				RC_z=ewa4i*(gixx_p*nRiem_xz_x+giyy_p*nRiem_yz_y);
				
				RijRij=(unnRiem_xx-uRC_xx)*(nnRiem_xx-RC_xx)
								+(unnRiem_yy-uRC_yy)*(nnRiem_yy-RC_yy)
								+(unnRiem_zz-uRC_zz)*(nnRiem_zz-RC_zz)
								-2.*ewa4i*(gizz_p*RC_z*RC_z)
								+RC*RC;
								
				RR=-2.*RC+(gxx_p*uRC_xx+gyy_p*uRC_yy+gzz_p*uRC_zz)/ewa4i;
				
				set_outv(l,k,j,1)=RijRij;
				set_outv(l,k,j,2)=RR*RR;
				set_outv(l,k,j,3)=Kinv-2.*RijRij+RR/3.;
			}
		}
	}
	
	
	#pragma omp barrier
	if(abs(KOep)>1.0e-10)
	KOdiss(lui);
	
	#pragma omp barrier
	//  4th order Runge-Kutta method
	//  add flux
	runge_kutta(dt0*frac);		//bvr+=dbv*dt
	
	#pragma omp barrier
	
	if(itype!=4)
	{
		new_bv(dt);				//bv=bv0+dbv*dt
	}
	else
	{
		new_bv4();				//bv=bv0+bvr
		
		if(get_exc())
		 excision_quadfunc();	//fill the excised region by a quadratic function
	}
	#pragma omp barrier

	enforce_const();
	//determinant of gamma=1 and tr A_ij=0
	#pragma omp barrier
	
	return;
}

void Fmv0::excision()
{
	int s=0,e=nn;
	
	#pragma omp parallel for 
	for(int l=lli+exg-6;l<lli+exg;l++)
	{
		int k=0;
		int j=0;
		
		for(int i=s;i<e;i++)
		 set_dbv(l,k,j,i)=get_dbv(l+1,k,j,i);
	}
	
	return;
}

void Fmv0::excision_quadfunc()
{
	#pragma omp parallel for
	for(int i=0;i<nn;i++)
	{
		if(!nonzero[i])
		continue;
		
		int j=0;
		int k=0;
		int l=lli+exg-6;

		double F=get_bv(l,k,j,i);
		double dF=-(3.*get_bv(l,k,j,i)-4.*get_bv(l+1,k,j,i)+get_bv(l+2,k,j,i))*dzi2;
		double x=get_z(l);
		
		if(!evod[i])
		{
			
			double a=(-F + dF*x)*pow(x,-2);
			double b=-dF + 2*F*pow(x,-1);
			
			for(int ll=lli;ll<lli+exg-6;ll++)
			{
				double z=get_z(ll);
				set_bv(ll,k,j,i)=a*pow(z,2)+b*z;
			}
		}
		else
		{
			double a=(dF*pow(x,-1))/2.;
			double c=F - (dF*x)/2.;
			
			for(int ll=lli;ll<lli+exg-6;ll++)
			{
				double z=get_z(ll);
				set_bv(ll,k,j,i)=a*pow(z,2)+c;
			}
		}
	}
	
	return;
}

//guarantee the determinant of gamma=1 and tracelessness of A_ij
//  -> psi(wa) and trK(ek) also change
void Fmv0::enforce_const()
{
	#pragma omp parallel for 
	for(int l=lli;l<=lui;l++)
	{
		int k=0;
		int j=0;
		
		if(get_bflag(l,k,j)!=0)
		 continue;

		double gxx_p,gyy_p,gzz_p;
		double akxx_p,akyy_p,akzz_p;
		double det_p,det_pi, det_p13;
		double gixx_p,giyy_p,gizz_p;
		double tmp_ek, tmp_ek3;
		 
		//  calculate determinant of tilde{\gamma_{ij}}
		gxx_p=get_bv(l,k,j, 7)+1.;
		gyy_p=get_bv(l,k,j, 8)+1.;
		gzz_p=get_bv(l,k,j, 9)+get_flat_df2z(l);
		det_p= gxx_p*gyy_p*gzz_p;

		if(det_p<1.e-16) 
		 det_p=1.e-16;
		det_pi=1./det_p;
		
		//det_p13=pow(det_p,(-1./3.));

		gixx_p= (gyy_p*gzz_p)*det_pi;
		giyy_p= gixx_p;
		gizz_p= (gxx_p*gyy_p)*det_pi;
		
		double fac=det_p/(get_flat_df2z(l));
		det_p13=pow(fac,(-1./3.));
		
		//set_bv(l,k,j,13)=get_bv(l,k,j,13)+log(det_p)/12.;
		set_bv(l,k,j,13)=get_bv(l,k,j,13)+log(fac)/12.;
		
		akxx_p=get_bv(l,k,j,14);
		akyy_p=get_bv(l,k,j,15);
		akzz_p=get_bv(l,k,j,16);

		tmp_ek=gixx_p*akxx_p+giyy_p*akyy_p+gizz_p*akzz_p;
		set_bv(l,k,j,20)=get_bv(l,k,j,20) +tmp_ek;
		
		set_bv(l,k,j, 7) = gxx_p*det_p13-1.;
		set_bv(l,k,j, 8) = gyy_p*det_p13-1.;
		set_bv(l,k,j, 9) = gzz_p*det_p13-get_flat_df2z(l);

		tmp_ek3=tmp_ek/3.;
		set_bv(l,k,j,14)=akxx_p*det_p13 -tmp_ek3*gxx_p*det_p13;
		set_bv(l,k,j,15)=akyy_p*det_p13 -tmp_ek3*gyy_p*det_p13;
		set_bv(l,k,j,16)=akzz_p*det_p13 -tmp_ek3*gzz_p*det_p13;
	}
	return;
}

void Fmv0::BSSN_adv()
{
	int s=0,e=24;

	#pragma omp parallel for 
	for(int l=lli+1;l<=lui;l++)
	{
		int k=0;
		int j=0;
		if(get_bflag(l,k,j)!=0)
		 continue;

		double bz_p;
		
		// beta^i
		bz_p=get_bv(l,k,j,3);
		
		for(int i=s;i<e;i++)
		{
			if(nonzero[i]!=true)
			continue;
			
			double adv_z;

			if(bz_p<0){
				adv_z=( -get_bv(l-3,k,j,i)
						+ 6.*get_bv(l-2,k,j,i)
						-18.*get_bv(l-1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l+1,k,j,i)
						)*dzi12*bz_p;
			}else{
				adv_z=-( -get_bv(l+3,k,j,i)
						+ 6.*get_bv(l+2,k,j,i)
						-18.*get_bv(l+1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l-1,k,j,i)
						)*dzi12*bz_p;
			}

			set_dbv(l,k,j,i)= adv_z;
		}
		
		if(scalarevo)
		{
			int i=nsc;
			double adv_z;
			
			if(bz_p<0){
				adv_z=( -get_bv(l-3,k,j,i)
						+ 6.*get_bv(l-2,k,j,i)
						-18.*get_bv(l-1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l+1,k,j,i)
						)*dzi12*bz_p;
			}else{
				adv_z=-( -get_bv(l+3,k,j,i)
						+ 6.*get_bv(l+2,k,j,i)
						-18.*get_bv(l+1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l-1,k,j,i)
						)*dzi12*bz_p;
			}

			set_dbv(l,k,j,i)= adv_z;
			
			i=nscp;
			
			if(bz_p<0){
				adv_z=( -get_bv(l-3,k,j,i)
						+ 6.*get_bv(l-2,k,j,i)
						-18.*get_bv(l-1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l+1,k,j,i)
						)*dzi12*bz_p;
			}else{
				adv_z=-( -get_bv(l+3,k,j,i)
						+ 6.*get_bv(l+2,k,j,i)
						-18.*get_bv(l+1,k,j,i)
						+10.*get_bv(l  ,k,j,i)
						+ 3.*get_bv(l-1,k,j,i)
						)*dzi12*bz_p;
			}

			set_dbv(l,k,j,i)= adv_z;
		}
	}// for-loop end

	return;
}

//Kreiss-Oliger dissipation term
void Fmv0::KOdiss(int lup)
{
	//only for geometry and scalar
	int s=0,e=24;
	
	#pragma omp parallel for 
	for(int l=lli;l<=lup;l++)
	{
		int k=0;
		int j=0;
		//double coef=KOep;
		
		if(get_bflag(l,k,j)!=0)
		 continue;
		
		double diss_z,diss_x,diss_y;
		
		for(int i=s;i<e;i++)
		{
			diss_x=( get_bv(l,k,j-3,i)
					 	-6.*get_bv(l,k,j-2,i)
							+ 15.*get_bv(l,k,j-1,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k,j+1,i)
						-6.*get_bv(l,k,j+2,i)
					+get_bv(l,k,j+3,i)
					)*dxi;
			diss_y=( get_bv(l,k-3,j,i)
					 	-6.*get_bv(l,k-2,j,i)
							+ 15.*get_bv(l,k-1,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k+1,j,i)
						-6.*get_bv(l,k+2,j,i)
					+get_bv(l,k+3,j,i)
					)*dyi;
			diss_z=( get_bv(l-3,k,j,i)
					 	-6.*get_bv(l-2,k,j,i)
							+ 15.*get_bv(l-1,k,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l+1,k,j,i)
						-6.*get_bv(l+2,k,j,i)
					+get_bv(l+3,k,j,i)
					)*dzi;
		        
			set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
		        
			//set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + coef*(diss_z)/64.;
		}
		
		//for scalar field
		if(scalarevo)
		{
			int i=nsc;
			diss_x=( get_bv(l,k,j-3,i)
					 	-6.*get_bv(l,k,j-2,i)
							+ 15.*get_bv(l,k,j-1,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k,j+1,i)
						-6.*get_bv(l,k,j+2,i)
					+get_bv(l,k,j+3,i)
					)*dxi;
			diss_y=( get_bv(l,k-3,j,i)
					 	-6.*get_bv(l,k-2,j,i)
							+ 15.*get_bv(l,k-1,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k+1,j,i)
						-6.*get_bv(l,k+2,j,i)
					+get_bv(l,k+3,j,i)
					)*dyi;
			diss_z=( get_bv(l-3,k,j,i)
					 	-6.*get_bv(l-2,k,j,i)
							+ 15.*get_bv(l-1,k,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l+1,k,j,i)
						-6.*get_bv(l+2,k,j,i)
					+get_bv(l+3,k,j,i)
					)*dzi;
		        
			set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
			
			
			i=nscp;
			diss_x=( get_bv(l,k,j-3,i)
					 	-6.*get_bv(l,k,j-2,i)
							+ 15.*get_bv(l,k,j-1,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k,j+1,i)
						-6.*get_bv(l,k,j+2,i)
					+get_bv(l,k,j+3,i)
					)*dxi;
			diss_y=( get_bv(l,k-3,j,i)
					 	-6.*get_bv(l,k-2,j,i)
							+ 15.*get_bv(l,k-1,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l,k+1,j,i)
						-6.*get_bv(l,k+2,j,i)
					+get_bv(l,k+3,j,i)
					)*dyi;
			diss_z=( get_bv(l-3,k,j,i)
					 	-6.*get_bv(l-2,k,j,i)
							+ 15.*get_bv(l-1,k,j,i)
								-20.*get_bv(l,k,j,i)
							+ 15.*get_bv(l+1,k,j,i)
						-6.*get_bv(l+2,k,j,i)
					+get_bv(l+3,k,j,i)
					)*dzi;
		        
			set_dbv(l,k,j,i)=get_dbv(l,k,j,i) + KOep*(diss_x+diss_y+diss_z)/64.;
		
		}
	}// for-loop end

	return;
}

void Fmv0::flux_fill()
{
	#pragma omp parallel for 
	for(int i=lli;i<=lui+1;i++)
	{
		int j,k,l=i;
		
		k=0;
		j=0;

		if(get_bflag(l,k,j)!=0)
		 continue;
		
		double zz;
		
		//metric variables
		double alp_ph,bx_ph,by_ph,bz_ph,wa_ph,gxx_ph,gyy_ph,gzz_ph,gxy_ph,gxz_ph,gyz_ph;
		double det,deti,sqdet;
		double gizz_ph,gixx_ph,giyy_ph;
		double fzz;
		
		//deviations of variable
		double Del2,Del1,Del0,Delpbar,Delmbar;
		
		//right and left dynamical variables 
		double S0L,S0R,SzL,SzR,DenL,DenR,rhoL,rhoR;
		double VxL,VxR,VyL,VyR,VzL,VzR,UxL,UxR,UyL,UyR,UzL,UzR,U2L,U2R;
		double cs2L,cs2R,ovfacL,ovfacR,lam1tL,lam2tL,lam1tR,lam2tR,lampL,lampR,lammL,lammR,as,prsL,prsR;

		zz=get_z(l)-0.5*dz;
	
		//substitution for geometrical variables start
		alp_ph=get_ipol_z_lower_mid(l,k,j,0);
		
		bz_ph=get_ipol_z_lower_mid(l,k,j,3);
		
		wa_ph=get_ipol_z_lower_mid(l,k,j,13);
		
		fzz=pow(df(zz),2);

		gxx_ph=get_ipol_z_lower_mid(l,k,j,7)+1.;
		gyy_ph=get_ipol_z_lower_mid(l,k,j,8)+1.;
		gzz_ph=get_ipol_z_lower_mid(l,k,j,9)+fzz;
		
		det=fzz;
		
		if(det<1.e-16) 
		 det=1.e-16;
		deti=1./det;
		gizz_ph= (gxx_ph*gyy_ph)*deti*exp(-4.*wa_ph);
		
		sqdet= sqrt(det)*exp(6.*wa_ph);
		//substitution for geometrical variables end

		//calculation of left and right S0 start
		Del2=get_bv(l+1,k,j,24)-get_bv(l,k,j,24);
		Del1=get_bv(l,k,j,24)-get_bv(l-1,k,j,24);
		Del0=get_bv(l-1,k,j,24)-get_bv(l-2,k,j,24);
		
		Delpbar=minmod(Del1,b_minmod*Del0);
		Delmbar=minmod(Del0,b_minmod*Del1);
		
		S0L=get_bv(l-1,k,j,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
		
		Delpbar=minmod(Del2,b_minmod*Del1);
		Delmbar=minmod(Del1,b_minmod*Del2);
		
		S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
		//calculation of left and right S0 end
		
		
		//calculation of left and right Sz start
		Del2=get_bv(l+1,k,j,27)-get_bv(l,k,j,27);
		Del1=get_bv(l,k,j,27)-get_bv(l-1,k,j,27);
		Del0=get_bv(l-1,k,j,27)-get_bv(l-2,k,j,27);
		
		Delpbar=minmod(Del1,b_minmod*Del0);
		Delmbar=minmod(Del0,b_minmod*Del1);
		
		SzL=get_bv(l-1,k,j,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
		
		Delpbar=minmod(Del2,b_minmod*Del1);
		Delmbar=minmod(Del1,b_minmod*Del2);
		
		SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
		//calculation of left and right Sz end
		
		//calculation of left and right Density start (Den:Density=\rho_* in my note)
		Del2=get_bv(l+1,k,j,28)-get_bv(l,k,j,28);
		Del1=get_bv(l,k,j,28)-get_bv(l-1,k,j,28);
		Del0=get_bv(l-1,k,j,28)-get_bv(l-2,k,j,28);
		
		Delpbar=minmod(Del1,b_minmod*Del0);
		Delmbar=minmod(Del0,b_minmod*Del1);
		
		DenL=get_bv(l-1,k,j,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
		
		Delpbar=minmod(Del2,b_minmod*Del1);
		Delmbar=minmod(Del1,b_minmod*Del2);
		
		DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
		//calculation of left and right Density end
		
		//calculation of left and right rho(primitive) start
		Del2=get_primv(l+1,k,j,0)-get_primv(l,k,j,0);
		Del1=get_primv(l,k,j,0)-get_primv(l-1,k,j,0);
		Del0=get_primv(l-1,k,j,0)-get_primv(l-2,k,j,0);
		
		Delpbar=minmod(Del1,b_minmod*Del0);
		Delmbar=minmod(Del0,b_minmod*Del1);
		
		rhoL=get_primv(l-1,k,j,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
		
		Delpbar=minmod(Del2,b_minmod*Del1);
		Delmbar=minmod(Del1,b_minmod*Del2);
		
		rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
		//calculation of left and right rho(primitive) end
	
		//calculation of left and right Vz(primitive) start
		Del2=get_primv(l+1,k,j,3)-get_primv(l,k,j,3);
		Del1=get_primv(l,k,j,3)-get_primv(l-1,k,j,3);
		Del0=get_primv(l-1,k,j,3)-get_primv(l-2,k,j,3);
		
		Delpbar=minmod(Del1,b_minmod*Del0);
		Delmbar=minmod(Del0,b_minmod*Del1);
		
		VzL=get_primv(l-1,k,j,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
		
		Delpbar=minmod(Del2,b_minmod*Del1);
		Delmbar=minmod(Del1,b_minmod*Del2);
		
		VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
		//calculation of left and right Vz(primitive) end

		UzL=(VzL+bz_ph)/alp_ph;
		UzR=(VzR+bz_ph)/alp_ph;
		
		U2L=(UzL*UzL*gzz_ph)*exp(4.*wa_ph);
		U2R=(UzR*UzR*gzz_ph)*exp(4.*wa_ph);

		//the speed should not exceed 1
		U2L=min(U2L,1.-1.0e-16);
		U2R=min(U2R,1.-1.0e-16);
		//calculation of left and right U_i end
		
		//sound velocity
		cs2L=dpres(rhoL);
		cs2R=dpres(rhoR);
		
		//calculation of characteristic speeds start
		ovfacL=alp_ph/(1.-U2L*cs2L);
		ovfacR=alp_ph/(1.-U2R*cs2R);
		
		lam1tL=UzL*(1-cs2L);
		lam2tL=sqrt(max(cs2L*(1-U2L)*(gizz_ph*(1-U2L*cs2L)-(1-cs2L)*UzL*UzL),1.0e-16));
		
		lam1tR=UzR*(1-cs2R);
		lam2tR=sqrt(max(cs2R*(1-U2R)*(gizz_ph*(1-U2R*cs2R)-(1-cs2R)*UzR*UzR),1.0e-16));
		
		lampL=ovfacL*(lam1tL+lam2tL)-bz_ph;
		lampR=ovfacR*(lam1tR+lam2tR)-bz_ph;
		
		lammL=ovfacL*(lam1tL-lam2tL)-bz_ph;
		lammR=ovfacR*(lam1tR-lam2tR)-bz_ph;
		
		as=max({abs(VzL),abs(VzR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
		//calculation of characteristic speeds end
		
		//pressure
		prsL=pres(rhoL);
		prsR=pres(rhoR);
		
		//fluxes
		set_flux_z(l,k,j,0)=0.5*(S0L*VzL+S0R*VzR+sqdet*(prsL*(VzL+bz_ph)+prsR*(VzR+bz_ph))-as*(S0R-S0L));
		set_flux_z(l,k,j,1)=0.5*(SzL*VzL+SzR*VzR+sqdet*alp_ph*(prsL+prsR)-as*(SzR-SzL));
		set_flux_z(l,k,j,2)=0.5*(DenL*VzL+DenR*VzR-as*(DenR-DenL));
		
		if(i==lui+1)
		continue;
		
		k=0;
		
		for(j=0;j<=1;j++)
		{
			//substitution for geometrical variables start
			alp_ph=get_ipol_x_lower_mid(l,k,j,0);
			
			bx_ph=get_ipol_x_lower_mid(l,k,j,1);
			by_ph=get_ipol_x_lower_mid(l,k,j,2);
			bz_ph=get_ipol_x_lower_mid(l,k,j,3);
			
			wa_ph=get_ipol_x_lower_mid(l,k,j,13);

			fzz=get_flat_df2z(l);
			
			gxx_ph=get_ipol_x_lower_mid(l,k,j,7)+1.;
			gyy_ph=get_ipol_x_lower_mid(l,k,j,8)+1.;
			gzz_ph=get_ipol_x_lower_mid(l,k,j,9)+fzz;
			
			gxy_ph=get_ipol_x_lower_mid(l,k,j,10);
			gxz_ph=get_ipol_x_lower_mid(l,k,j,11);
			gyz_ph=get_ipol_x_lower_mid(l,k,j,12);
		
			det=fzz;

			deti=1./det;
			gixx_ph= (gyy_ph*gzz_ph-gyz_ph*gyz_ph)*deti*exp(-4.*wa_ph);
			
			sqdet= sqrt(det)*exp(6.*wa_ph);
			//substitution for geometrical variables end
			
			//calculation of left and right S0 start
			Del2=get_bv(l,k,j+1,24)-get_bv(l,k,j,24);
			Del1=get_bv(l,k,j,24)-get_bv(l,k,j-1,24);
			Del0=get_bv(l,k,j-1,24)-get_bv(l,k,j-2,24);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			S0L=get_bv(l,k,j-1,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right S0 end
			
			//calculation of left and right Sz start
			Del2=get_bv(l,k,j+1,27)-get_bv(l,k,j,27);
			Del1=get_bv(l,k,j,27)-get_bv(l,k,j-1,27);
			Del0=get_bv(l,k,j-1,27)-get_bv(l,k,j-2,27);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			SzL=get_bv(l,k,j-1,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Sz end

			//calculation of left and right Density start (Den:Density=\rho_* in my note)
			Del2=get_bv(l,k,j+1,28)-get_bv(l,k,j,28);
			Del1=get_bv(l,k,j,28)-get_bv(l,k,j-1,28);
			Del0=get_bv(l,k,j-1,28)-get_bv(l,k,j-2,28);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			DenL=get_bv(l,k,j-1,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Density end
			
			//calculation of left and right rho(primitive) start
			Del2=get_primv(l,k,j+1,0)-get_primv(l,k,j,0);
			Del1=get_primv(l,k,j,0)-get_primv(l,k,j-1,0);
			Del0=get_primv(l,k,j-1,0)-get_primv(l,k,j-2,0);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			rhoL=get_primv(l,k,j-1,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right rho(primitive) end
			
			
			//calculation of left and right Vx(primitive) start
			Del2=get_primv(l,k,j+1,1)-get_primv(l,k,j,1);
			Del1=get_primv(l,k,j,1)-get_primv(l,k,j-1,1);
			Del0=get_primv(l,k,j-1,1)-get_primv(l,k,j-2,1);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VxL=get_primv(l,k,j-1,1)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VxR=get_primv(l,k,j,1)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vx(primitive) end
			
			//calculation of left and right Vy(primitive) start
			Del2=get_primv(l,k,j+1,2)-get_primv(l,k,j,2);
			Del1=get_primv(l,k,j,2)-get_primv(l,k,j-1,2);
			Del0=get_primv(l,k,j-1,2)-get_primv(l,k,j-2,2);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VyL=get_primv(l,k,j-1,2)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VyR=get_primv(l,k,j,2)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vx(primitive) end
			
			
			//calculation of left and right Vz(primitive) start
			Del2=get_primv(l,k,j+1,3)-get_primv(l,k,j,3);
			Del1=get_primv(l,k,j,3)-get_primv(l,k,j-1,3);
			Del0=get_primv(l,k,j-1,3)-get_primv(l,k,j-2,3);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VzL=get_primv(l,k,j-1,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vz(primitive) end

			//calculation of left and right epsilon(primitive) start
			//Del2=get_primv(l,k,j+2,4)-get_primv(l,k,j+1,4);
			//Del1=get_primv(l,k,j+1,4)-get_primv(l,k,j,4);
			//Del0=get_primv(l,k,j,4)-get_primv(l,k,j-1,4);
			
			//Delpbar=minmod(Del1,b_minmod*Del0);
			//Delmbar=minmod(Del0,b_minmod*Del1);
			
			//epsLx=get_primv(l,k,j,4)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			//Delpbar=minmod(Del2,b_minmod*Del1);
			//Delmbar=minmod(Del1,b_minmod*Del2);
			
			//epsRx=get_primv(l,k,j+1,4)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right epsilon(primitive) end
			
			//calculation of left and right U_i start
			UxL=(VxL+bx_ph)/alp_ph;
			UxR=(VxR+bx_ph)/alp_ph;
			UyL=(VyL+by_ph)/alp_ph;
			UyR=(VyR+by_ph)/alp_ph;
			UzL=(VzL+bz_ph)/alp_ph;
			UzR=(VzR+bz_ph)/alp_ph;
			
			U2L=(UxL*UxL*gxx_ph+UyL*UyL*gyy_ph+UzL*UzL*gzz_ph
				+2.*(UxL*UyL*gxy_ph+UxL*UzL*gxz_ph+UyL*UzL*gyz_ph))*exp(4.*wa_ph);
			U2R=(UxR*UxR*gxx_ph+UyR*UyR*gyy_ph+UzR*UzR*gzz_ph
				+2.*(UxR*UyR*gxy_ph+UxR*UzR*gxz_ph+UyR*UzR*gyz_ph))*exp(4.*wa_ph);
			
			//the speed should not exceed 1
			U2L=min(U2L,1.-1.0e-16);
			U2R=min(U2R,1.-1.0e-16);
			//calculation of left and right U_i end
			
			//sound velocity
			cs2L=dpres(rhoL);
			cs2R=dpres(rhoR);
			
			//calculation of characteristic speeds start
			ovfacL=alp_ph/(1.-U2L*cs2L);
			ovfacR=alp_ph/(1.-U2R*cs2R);
			
			lam1tL=UxL*(1-cs2L);
			lam2tL=sqrt(max(cs2L*(1-U2L)*(gixx_ph*(1-U2L*cs2L)-(1-cs2L)*UxL*UxL),1.0e-16));

			lam1tR=UxR*(1-cs2R);
			lam2tR=sqrt(max(cs2R*(1-U2R)*(gixx_ph*(1-U2R*cs2R)-(1-cs2R)*UxR*UxR),1.0e-16));
			
			lampL=ovfacL*(lam1tL+lam2tL)-bx_ph;
			lampR=ovfacR*(lam1tR+lam2tR)-bx_ph;
			
			lammL=ovfacL*(lam1tL-lam2tL)-bx_ph;
			lammR=ovfacR*(lam1tR-lam2tR)-bx_ph;
			
			as=max({abs(VxL),abs(VxR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
			//calculation of characteristic speeds end
			
			//pressure
			prsL=pres(rhoL);
			prsR=pres(rhoR);
			
			//fluxes
			set_flux_x(l,k,j,0)=0.5*(S0L*VxL+S0R*VxR+sqdet*(prsL*(VxL+bx_ph)+prsR*(VxR+bx_ph))-as*(S0R-S0L));
			set_flux_x(l,k,j,1)=0.5*(SzL*VxL+SzR*VxR-as*(SzR-SzL));
			set_flux_x(l,k,j,2)=0.5*(DenL*VxL+DenR*VxR-as*(DenR-DenL));
		}

		j=0;
		
		for(k=0;k<=1;k++)
		{
			//substitution for geometrical variables start
			alp_ph=get_ipol_y_lower_mid(l,k,j,0);
			
			bx_ph=get_ipol_y_lower_mid(l,k,j,1);
			by_ph=get_ipol_y_lower_mid(l,k,j,2);
			bz_ph=get_ipol_y_lower_mid(l,k,j,3);
			
			wa_ph=get_ipol_y_lower_mid(l,k,j,13);
			
			fzz=get_flat_df2z(l);
			
			gxx_ph=get_ipol_y_lower_mid(l,k,j,7)+1.;
			gyy_ph=get_ipol_y_lower_mid(l,k,j,8)+1.;
			gzz_ph=get_ipol_y_lower_mid(l,k,j,9)+fzz;
			gxy_ph=get_ipol_y_lower_mid(l,k,j,10);
			gxz_ph=get_ipol_y_lower_mid(l,k,j,11);
			gyz_ph=get_ipol_y_lower_mid(l,k,j,12);
		
			det=fzz;
			
			if(det<1.e-16) 
			 det=1.e-16;
			deti=1./det;
			giyy_ph= (gxx_ph*gzz_ph-gxz_ph*gxz_ph)*deti*exp(-4.*wa_ph);
			
			sqdet= sqrt(det)*exp(6.*wa_ph);
			//substitution for geometrical variables end
			
			//calculation of left and right S0 start
			Del2=get_bv(l,k+1,j,24)-get_bv(l,k,j,24);
			Del1=get_bv(l,k,j,24)-get_bv(l,k-1,j,24);
			Del0=get_bv(l,k-1,j,24)-get_bv(l,k-2,j,24);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			S0L=get_bv(l,k-1,j,24)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			S0R=get_bv(l,k,j,24)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right S0 end
			
			//calculation of left and right Sz start
			Del2=get_bv(l,k+1,j,27)-get_bv(l,k,j,27);
			Del1=get_bv(l,k,j,27)-get_bv(l,k-1,j,27);
			Del0=get_bv(l,k-1,j,27)-get_bv(l,k-2,j,27);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			SzL=get_bv(l,k-1,j,27)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			SzR=get_bv(l,k,j,27)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Sz end
			
			//calculation of left and right Density start (Den:Density=\rho_* in my note)
			Del2=get_bv(l,k+1,j,28)-get_bv(l,k,j,28);
			Del1=get_bv(l,k,j,28)-get_bv(l,k-1,j,28);
			Del0=get_bv(l,k-1,j,28)-get_bv(l,k-2,j,28);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			DenL=get_bv(l,k-1,j,28)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			DenR=get_bv(l,k,j,28)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Density end
			
			//calculation of left and right rho(primitive) start
			Del2=get_primv(l,k+1,j,0)-get_primv(l,k,j,0);
			Del1=get_primv(l,k,j,0)-get_primv(l,k-1,j,0);
			Del0=get_primv(l,k-1,j,0)-get_primv(l,k-2,j,0);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			rhoL=get_primv(l,k-1,j,0)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			rhoR=get_primv(l,k,j,0)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right rho(primitive) end
			
			
			//calculation of left and right Vx(primitive) start
			Del2=get_primv(l,k+1,j,1)-get_primv(l,k,j,1);
			Del1=get_primv(l,k,j,1)-get_primv(l,k-1,j,1);
			Del0=get_primv(l,k-1,j,1)-get_primv(l,k-2,j,1);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VxL=get_primv(l,k-1,j,1)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VxR=get_primv(l,k,j,1)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vx(primitive) end
			
			
			//calculation of left and right Vy(primitive) start
			Del2=get_primv(l,k+1,j,2)-get_primv(l,k,j,2);
			Del1=get_primv(l,k,j,2)-get_primv(l,k-1,j,2);
			Del0=get_primv(l,k-1,j,2)-get_primv(l,k-2,j,2);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VyL=get_primv(l,k-1,j,2)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VyR=get_primv(l,k,j,2)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vy(primitive) end
			
			
			//calculation of left and right Vz(primitive) start
			Del2=get_primv(l,k+1,j,3)-get_primv(l,k,j,3);
			Del1=get_primv(l,k,j,3)-get_primv(l,k-1,j,3);
			Del0=get_primv(l,k-1,j,3)-get_primv(l,k-2,j,3);
			
			Delpbar=minmod(Del1,b_minmod*Del0);
			Delmbar=minmod(Del0,b_minmod*Del1);
			
			VzL=get_primv(l,k-1,j,3)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			Delpbar=minmod(Del2,b_minmod*Del1);
			Delmbar=minmod(Del1,b_minmod*Del2);
			
			VzR=get_primv(l,k,j,3)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right Vz(primitive) end

			//calculation of left and right epsilon(primitive) start
			//Del2=get_primv(l,k,j+2,4)-get_primv(l,k,j+1,4);
			//Del1=get_primv(l,k,j+1,4)-get_primv(l,k,j,4);
			//Del0=get_primv(l,k,j,4)-get_primv(l,k,j-1,4);
			
			//Delpbar=minmod(Del1,b_minmod*Del0);
			//Delmbar=minmod(Del0,b_minmod*Del1);
			
			//epsLx=get_primv(l,k,j,4)+0.25*((1.-kap_MUSCL)*Delmbar+(1.+kap_MUSCL)*Delpbar);
			
			//Delpbar=minmod(Del2,b_minmod*Del1);
			//Delmbar=minmod(Del1,b_minmod*Del2);
			
			//epsRx=get_primv(l,k,j+1,4)-0.25*((1.-kap_MUSCL)*Delpbar+(1.+kap_MUSCL)*Delmbar);
			//calculation of left and right epsilon(primitive) end
			
			//calculation of left and right U_i start
			UxL=(VxL+bx_ph)/alp_ph;
			UxR=(VxR+bx_ph)/alp_ph;
			UyL=(VyL+by_ph)/alp_ph;
			UyR=(VyR+by_ph)/alp_ph;
			UzL=(VzL+bz_ph)/alp_ph;
			UzR=(VzR+bz_ph)/alp_ph;
			
			U2L=(UxL*UxL*gxx_ph+UyL*UyL*gyy_ph+UzL*UzL*gzz_ph
				+2.*(UxL*UyL*gxy_ph+UxL*UzL*gxz_ph+UyL*UzL*gyz_ph))*exp(4.*wa_ph);
			U2R=(UxR*UxR*gxx_ph+UyR*UyR*gyy_ph+UzR*UzR*gzz_ph
				+2.*(UxR*UyR*gxy_ph+UxR*UzR*gxz_ph+UyR*UzR*gyz_ph))*exp(4.*wa_ph);
			
			//the speed should not exceed 1
			U2L=min(U2L,1.-1.0e-16);
			U2R=min(U2R,1.-1.0e-16);
			//calculation of left and right U_i end

			//the speed should not exceed 1
			cs2L=dpres(rhoL);
			cs2R=dpres(rhoR);
			
			//calculation of characteristic speeds start
			ovfacL=alp_ph/(1.-U2L*cs2L);
			ovfacR=alp_ph/(1.-U2R*cs2R);
			
			lam1tL=UyL*(1-cs2L);
			lam2tL=sqrt(max(cs2L*(1-U2L)*(giyy_ph*(1-U2L*cs2L)-(1-cs2L)*UyL*UyL),1.0e-16));

			lam1tR=UyR*(1-cs2R);
			lam2tR=sqrt(max(cs2R*(1-U2R)*(giyy_ph*(1-U2R*cs2R)-(1-cs2R)*UyR*UyR),1.0e-16));
			
			lampL=ovfacL*(lam1tL+lam2tL)-by_ph;
			lampR=ovfacR*(lam1tR+lam2tR)-by_ph;
			
			lammL=ovfacL*(lam1tL-lam2tL)-by_ph;
			lammR=ovfacR*(lam1tR-lam2tR)-by_ph;
			
			as=max({abs(VyL),abs(VyR),abs(lampL),abs(lampR),abs(lammL),abs(lammR)});
			//calculation of characteristic speeds end
			
			//pressure
			prsL=pres(rhoL);
			prsR=pres(rhoR);
			
			//fluxes
			set_flux_y(l,k,j,0)=0.5*(S0L*VyL+S0R*VyR+sqdet*(prsL*(VyL+by_ph)+prsR*(VyR+by_ph))-as*(S0R-S0L));
			set_flux_y(l,k,j,1)=0.5*(SzL*VyL+SzR*VyR-as*(SzR-SzL));
			set_flux_y(l,k,j,2)=0.5*(DenL*VyL+DenR*VyR-as*(DenR-DenL));
		}
	}
}
