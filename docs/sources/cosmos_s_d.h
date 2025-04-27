/**
 * @file cosmos_s.h
 * @brief Header file defining the core classes for the COSMOS_S BSSN evolution code.
 * @version 1.00
 * @author Chulmoon Yoo
 *
 * @details
 * This file defines three main classes:
 * - Fmv0: The base class holding common grid data, simulation parameters, state variables,
 *         and methods for evolution, diagnostics, I/O, and interpolation. It manages a single grid level.
 * - Fmv: Derived from Fmv0, specialized for the base grid (Layer 0). Includes methods for
 *        asymptotic boundary conditions and specific initial data setups.
 * - Fmv1: Derived from Fmv0, specialized for refinement levels (Layer > 0). Includes methods
 *         for handling boundaries with parent/child layers in the FMR hierarchy.
 *
 * The code uses the BSSN formulation for evolving Einstein's equations, potentially coupled
 * with perfect fluid hydrodynamics and a scalar field. It supports Fixed Mesh Refinement (FMR).
 */

 #ifndef _COSMOS_S_H_
 #define _COSMOS_S_H_
 
 #include <cmath>
 #include <iostream>
 #include <iomanip>
 #include <fstream>
 #include <cstdio>
 #include <cstdlib>
 
 // Define standard constants if not already defined (common practice)
 #ifndef M_PI
 /** @brief Mathematical constant Pi. */
 #define M_PI 3.14159265358979323846
 #endif
 #ifndef M_E
 /** @brief Mathematical constant e (Euler's number). */
 #define M_E 2.71828182845904523536
 #endif
 
 using namespace std;
 
 // --- Forward declarations if needed ---
 class Fmv0;
 class Fmv;
 class Fmv1;
 
 /**
  * @class Fmv0
  * @brief Base class for managing a single grid level in the COSMOS_S simulation.
  *
  * @details
  * This class encapsulates the data and functionality for a single grid level.
  * It stores simulation parameters, grid geometry, state variables (metric, gauge,
  * fluid, scalar field), temporary storage for Runge-Kutta steps and interpolation,
  * and provides methods for initialization, evolution (BSSN, hydro, scalar),
  * boundary conditions (within the level), diagnostics (constraints, horizons),
  * interpolation, and basic I/O.
  *
  * Memory for grid functions (like bv, dbv, etc.) is dynamically allocated
  * in the constructor and deallocated in the destructor.
  * Grid functions are typically accessed using multi-dimensional arrays, often
  * with index offsets (e.g., `l-lmin`).
  */
 class Fmv0{
 
 protected:
     // --- General parameters and variables ---
     /** @brief Flag: Enable fluid evolution (true/false). */
     bool fluidevo;
     /** @brief Flag: Enable scalar field evolution (true/false). */
     bool scalarevo;
     /** @brief Flag: Enable curvature invariant evaluation (true/false). */
     bool curveval;
     /** @brief Flag: Horizon formation detected (true/false). */
     bool hform;
     /** @brief Flag: Excision is active (true/false). */
     bool exc;
     /** @brief Flag: Mesh refinement is active on this layer's parent (true if this layer has a child). */
     bool mrf; // mesh refinement flag (for the layer below this one)
 
     /** @brief Number of buffer cells in z-direction (radial/vertical). */
     int tab;
     /** @brief Number of buffer cells in x,y-directions (azimuthal/horizontal). */
     int tabx;
     /** @brief Max/min grid indices (including buffers) [j=x, k=y, l=z]. */
     int jmax,jmin,kmax,kmin,lmax,lmin;
     /** @brief Max/min grid indices for evolution domain (excluding buffers) [j=x, k=y, l=z]. */
     int jui,jli,kui,kli,lui,lli;
     /** @brief Total grid points (including buffers) in x, y, z. */
     int nx,ny,nz;
     /** @brief Total number of evolution variables. */
     int nn;
     /** @brief Number of constraint variables stored in `con`. */
     int nc;
     /** @brief Number of output variables stored in `outv`. */
     int nv;
     /** @brief Number of coefficients (?). Unused? */
     int ncoef;
     /** @brief Number of energy-momentum tensor components (?). Unused? */
     int nt;
     /** @brief Number of primitive fluid variables stored in `primv`. */
     int npr;
     /** @brief Index in `bv` for the scalar field variable phi. */
     int nsc;
     /** @brief Index in `bv` for the scalar field variable Pi. */
     int nscp;
     /** @brief Number of variables stored in `tmpipol` for interpolation. */
     int ntmpipol;
     /** @brief Grid index l of maximum Hamiltonian constraint violation. */
     int lhm;
     /** @brief Grid index l of maximum Momentum constraint violation. */
     int lmm;
     /** @brief Grid index l of maximum Kretschmann invariant value. */
     int lkm;
     /** @brief Grid index l of maximum Weyl invariant value. */
     int lwm;
     /** @brief Order of Lagrange interpolation used. */
     int ipolo;
     /** @brief Upper buffer size needed for interpolation stencil. */
     int ipolbufu;
     /** @brief Lower buffer size needed for interpolation stencil. */
     int ipolbufl;
     /** @brief Grid index offset from lli defining the inner boundary for excision. */
     int exg;
     /** @brief Grid index l just outside the apparent horizon. */
     int hl;
     /** @brief Number of apparent horizons found so far. */
     int hn;
     /** @brief Number of necks found so far (?). */
     int nneck;
     /** @brief Maximum number of horizons to track. */
     int hnmax;
 
     /** @brief Layer number in the FMR hierarchy (0 for base grid). */
     int layn;
     /** @brief Number of grid points to neglect near outer boundary for constraint checks. */
     int negb;
     /** @brief Minimum grid index l for constraint checking (can be adjusted, e.g., for excision). */
     int llmin;
 
     /** @brief Max/min physical coordinate values [x, y, z]. */
     double xu,xl, yu,yl, zu,zl;
     /** @brief Grid spacing (physical) in x, y, z and their inverses. */
     double dx,dy,dz,dxi,dyi,dzi;
     /** @brief Half inverses of grid spacing. */
     double dxi2,dyi2,dzi2;
     /** @brief Quarter inverses of grid spacing (?). */
     double dxi4,dyi4,dzi4;
     /** @brief Inverses divided by 12 (for 4th order finite difference). */
     double dxi12,dyi12,dzi12;
     /** @brief Inverses divided by 24 (for 4th order finite difference). */
     double dxi24,dyi24,dzi24;
     /** @brief Cell volume (dx*dy*dz). */
     double dvol;
     /** @brief Current simulation time. */
     double t;
     /** @brief Current time step size (potentially layer-dependent). */
     double dt;
     /** @brief Base time step size (usually of the coarsest level). */
     double dt0;
     /** @brief Previous time step size. */
     double dtp;
     /** @brief Time step size before previous. */
     double dtpp;
     /** @brief Maximum simulation time allowed. */
     double tmax;
     /** @brief CFL factor. */
     double cfl;
     /** @brief Gauge parameter eta for lapse (1+log slicing). */
     double etaa;
     /** @brief Gauge parameter eta for shift (Gamma driver B evolution). */
     double etab;
     /** @brief Gauge parameter eta_beta for shift (Gamma driver beta evolution). */
     double etabb;
     /** @brief Kreiss-Oliger dissipation coefficient. */
     double KOep;
     /** @brief Amplitude for inhomogeneous grid mapping (z-axis). */
     double amp;
     /** @brief Cosmological constant lambda. */
     double lambda;
     /** @brief 2*Pi, 4*Pi, 8*Pi, 16*Pi, 32*Pi constants. */
     double pi2,pi4,pi8,pi16,pi32;
 
     // --- Constraints and curvature invariants ---
     /** @brief Average Hamiltonian constraint violation (normalized). */
     double ham;
     /** @brief Maximum Hamiltonian constraint violation (normalized). */
     double hammax;
     /** @brief Average Momentum constraint violation (normalized, z-component). */
     double mom;
     /** @brief Maximum Momentum constraint violation (normalized, z-component). */
     double mommax;
     /** @brief Average deviation of Gamma constraint (Gamma^z - D_gamma^z). */
     double dGam;
     /** @brief Maximum deviation of Gamma constraint (Gamma^z - D_gamma^z). */
     double dGammax;
     /** @brief Maximum Kretschmann invariant value. */
     double Kremax;
     /** @brief Maximum Weyl invariant value (?). */
     double Weylmax;
 
     // --- Fluid parameters ---
     /** @brief Kappa parameter for MUSCL reconstruction. */
     double kap_MUSCL;
     /** @brief 'b' parameter for minmod slope limiter. */
     double b_minmod;
     /** @brief Fluid equation of state parameter (w = P/rho). */
     double fluidw;
 
     // --- Scalar field parameters ---
     /** @brief Scalar field mass parameter (m). Potential V = 0.5*(m*phi)^2. */
     double scalarm;
 
     // --- Parameters for initial data and settings ---
     /** @brief Initial Hubble parameter H_initial. */
     double Hb;
     /** @brief Initial simulation time t_initial. */
     double tini;
     /** @brief Initial value of (tr K)^2 (temporary/unused?). */
     double tk2;
     /** @brief Box size parameter for inhomogeneous grid mapping. */
     double boxL;
 
     // --- Coordinates and buffers ---
     /** @brief Upper grid index for buffer region `n` [j=x, k=y, l=z]. */
     int *ju,*jl,*ku,*kl,*lu,*ll;
     /** @brief Array of physical x-coordinates. Size nx. */
     double *x;
     /** @brief Array of physical y-coordinates. Size ny. */
     double *y;
     /** @brief Array of physical z-coordinates (computational grid). Size nz. */
     double *z;
 
     // --- Boundary flags ---
     /**
      * @brief Boundary flag array `bflag[l-lmin][k-kmin][j-jmin]`.
      * @details Stores information about grid point type (e.g., interior, physical boundary, symmetry boundary, excision boundary).
      *          Values might indicate distance to excision boundary (1 to 6) or -1 for inside excision. 0 for normal evolution points.
      */
     int ***bflag;
     /**
      * @brief Horizon flag array `hflag[l-lmin][k-kmin][j-jmin]`.
      * @details Indicates if a grid point is inside an apparent horizon (e.g., 1 if inside, 0 otherwise).
      */
     int ***hflag;
 
     // --- Dynamical variables ---
     /**
      * @brief Main state variables `bv[var_idx][z_idx][y_idx][x_idx]`.
      * @details Stores the current values of all evolved fields.
      * Indices `z_idx = l-lmin`, `y_idx = k-kmin`, `x_idx = j-jmin`.
      * Variable indices (`var_idx`) are:
      * - 0: alpha (lapse)
      * - 1: betax (shift x)
      * - 2: betay (shift y)
      * - 3: betaz (shift z)
      * - 4: Bx (Gamma-driver auxiliary var x)
      * - 5: By (Gamma-driver auxiliary var y)
      * - 6: Bz (Gamma-driver auxiliary var z)
      * - 7: gxx (conformal metric xx component - 1)
      * - 8: gyy (conformal metric yy component - 1)
      * - 9: gzz (conformal metric zz component - fzz)
      * - 10: gxy (conformal metric xy component)
      * - 11: gxz (conformal metric xz component)
      * - 12: gyz (conformal metric yz component)
      * - 13: wa (Log of conformal factor W = exp(wa) = chi^(-1/4))
      * - 14: Akxx (Trace-free extrinsic curvature xx)
      * - 15: Akyy (Trace-free extrinsic curvature yy)
      * - 16: Akzz (Trace-free extrinsic curvature zz)
      * - 17: Akxy (Trace-free extrinsic curvature xy)
      * - 18: Akxz (Trace-free extrinsic curvature xz)
      * - 19: Akyz (Trace-free extrinsic curvature yz)
      * - 20: trK (Trace of extrinsic curvature K)
      * - 21: Gammabarx (BSSN Gamma variable x)
      * - 22: Gammabary (BSSN Gamma variable y)
      * - 23: Gammabarz (BSSN Gamma variable z)
      * --- If `fluidevo` is true: ---
      * - 24: E (Conserved energy density T_nn * sqrt(gamma))
      * - 25: Sx (Conserved momentum density x, S_x * sqrt(gamma))
      * - 26: Sy (Conserved momentum density y, S_y * sqrt(gamma))
      * - 27: Sz (Conserved momentum density z, S_z * sqrt(gamma))
      * - 28: D (Conserved rest mass density rho_* * sqrt(gamma))
      * --- If `scalarevo` is true (indices depend on `fluidevo`): ---
      * - nsc (29 or 24): phi (Scalar field value)
      * - nscp (30 or 25): Pi (Scalar field momentum)
      */
     double ****bv;
     /** @brief RHS storage `dbv[var_idx][z_idx][y_idx][x_idx]`. Stores d(bv)/dt. */
     double ****dbv;
     /** @brief State variables at previous time step `bv0[var_idx][z_idx][y_idx][x_idx]`. */
     double ****bv0;
     /** @brief State variables two time steps before `bv1[var_idx][z_idx][y_idx][x_idx]`. */
     double ****bv1;
     /** @brief Storage for Runge-Kutta intermediate steps `bvr[var_idx][z_idx][y_idx][x_idx]`. */
     double ****bvr;
 
     // --- Temporary storage for interpolation ---
     /** @brief Storage for interpolated values `tmpipol[interp_var_idx][z_idx][y_idx][x_idx]`. */
     double ****tmpipol;
 
     // --- Fluid fluxes and primitive variables ---
     /** @brief Fluid fluxes in x-direction `flux_x[flux_idx][z_idx][y_idx][x_idx]`. Stored at i-1/2 interface. */
     double ****flux_x;
     /** @brief Fluid fluxes in y-direction `flux_y[flux_idx][z_idx][y_idx][x_idx]`. Stored at j-1/2 interface. */
     double ****flux_y;
     /** @brief Fluid fluxes in z-direction `flux_z[flux_idx][z_idx][y_idx][x_idx]`. Stored at k-1/2 interface. */
     double ****flux_z;
     /**
      * @brief Primitive fluid variables `primv[prim_idx][z_idx][y_idx][x_idx]`.
      * @details Indices `prim_idx`:
      * - 0: rho (rest mass density)
      * - 1: Vx (velocity x)
      * - 2: Vy (velocity y)
      * - 3: Vz (velocity z)
      * - 4: epsilon (specific internal energy)
      */
     double ****primv;
 
     // --- Reference coordinate values ---
     /**
      * @brief Reference z-coordinate `refz[z_idx][y_idx][x_idx]`.
      * @details Stores the corresponding z-coordinate on the axis (j=0, k=0) for interpolation
      *          in spherical symmetry or axisymmetry off-axis.
      */
     double ***refz;
 
     // --- Geometrical variables for inhomogeneous coordinate ---
     /** @brief Array of computational Z coordinates `coordZ[z_idx]`. Z = funcf(z). */
     double *coordZ;
     /** @brief Array of f_zz = (dZ/dz)^2 values `flat_df2z[z_idx]`. */
     double *flat_df2z;
     /** @brief Array of bar{Gamma}^z_{zz} values `flat_Gamz[z_idx]`. */
     double *flat_Gamz;
     /** @brief Array of d/dz (bar{Gamma}^z_{zz}) values `flat_dGamz[z_idx]`. */
     double *flat_dGamz;
 
     // --- Constraints and output storages ---
     /**
      * @brief Storage for constraint violations `con[con_idx][z_idx][y_idx][x_idx]`.
      * @details Indices `con_idx`:
      * - 0: Normalized Hamiltonian constraint
      * - 1: Hamiltonian constraint
      * - 2: Normalized Momentum constraint (z-component)
      * - 3: Momentum constraint (z-component)
      * - 4: Gamma constraint (Gamma^z - D_gamma^z)
      */
     double ****con;
     /**
      * @brief Storage for output variables `outv[out_idx][z_idx][y_idx][x_idx]`.
      * @details Used to store derived quantities for output or diagnostics (e.g., invariants, horizon info).
      * Indices `out_idx` (based on usage in code):
      * - 0: Kretschmann invariant (?)
      * - 1: R_ij R^ij (?)
      * - 2: R^2 (?)
      * - 3: Weyl invariant (?) C_ij C^ij (?)
      * - 4: Outgoing null expansion theta_+ (?)
      * - 5: Area radius R_areal
      * - 6: Ingoing null expansion theta_- (?)
      * - 7: Scalar field energy density (?)
      * - 8: Scalar field momentum density z (?)
      * - 9: Kodama mass M_Kodama
      * - 10: Ratio of scalar energy to fluid energy (?)
      * - 11: d(R_areal)/dr (?)
      */
     double ****outv;
 
     // --- Vectors for spline interpolation ---
     /** @brief Lower diagonal vector for spline matrix (even parity). */
     double *alp_e;
     /** @brief Upper diagonal vector for spline matrix (even parity). */
     double *bet_e;
     /** @brief Diagonal vector for spline matrix (even parity). */
     double *gam_e;
     /** @brief Lower diagonal vector for spline matrix (odd parity). */
     double *alp_o;
     /** @brief Upper diagonal vector for spline matrix (odd parity). */
     double *bet_o;
     /** @brief Diagonal vector for spline matrix (odd parity). */
     double *gam_o;
     /** @brief Storage for second derivatives `ddy[var_idx][z_idx]` used in spline interpolation. */
     double **ddy;
 
     // --- Supplements ---
     /** @brief Flag `nonzero[var_idx]` indicating if variable `i` is non-zero on the axis (used for advection/dissipation?). */
     bool *nonzero;
     /** @brief Flag `evod[var_idx]` indicating if variable `i` has even (true) or odd (false) parity at the center. */
     bool *evod;
     /** @brief Storage for apparent horizon radii `horis[horizon_idx][0=coord, 1=areal]`. */
     double **horis;
     /** @brief Storage for neck radii `neck[neck_idx]`. */
     double *neck;
 
 public:
     /**
      * @brief Constructor for the Fmv0 base class.
      * @param tabsz Number of buffer cells in z-direction.
      * @param tabsx Number of buffer cells in x,y-directions.
      * @param lupper Upper grid index (evolution domain) in z-direction.
      * @param xupper Max physical coordinate in x (used with zmax, nzmax etc to set dx).
      * @param zupper Max physical coordinate in z (evolution domain).
      * @param am Amplitude for inhomogeneous grid mapping.
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature invariant evaluation.
      */
     Fmv0(int tabsz,int tabsx,int lupper,
     double xupper,double zupper,double am,bool fld, bool scl, bool cuev); // Constructor definition follows...
 
     /**
      * @brief Virtual destructor for the Fmv0 base class.
      * @details Deallocates all dynamically allocated memory for grid functions and auxiliary arrays.
      */
     virtual ~Fmv0(); // Destructor definition follows...
 
 public:
 
     // --- GET functions ---
     // Provide read-only access to protected members.
 
     /** @brief Get fluid evolution flag. @return True if fluid evolution is enabled. */
     bool get_fluidevo() const { return fluidevo; }
     /** @brief Get horizon formation flag. @return True if a horizon has been detected. */
     bool get_hform() const { return hform; }
     /** @brief Get excision flag. @return True if excision is active. */
     bool get_exc() const { return exc; }
     /** @brief Get mesh refinement flag. @return True if this layer has a child layer. */
     bool get_mrf() const { return mrf; }
 
     /** @brief Get max grid index j (x-dir, incl. buffers). @return jmax. */
     int get_jmax() const { return jmax; }
     /** @brief Get min grid index j (x-dir, incl. buffers). @return jmin. */
     int get_jmin() const { return jmin; }
     /** @brief Get max grid index k (y-dir, incl. buffers). @return kmax. */
     int get_kmax() const { return kmax; }
     /** @brief Get min grid index k (y-dir, incl. buffers). @return kmin. */
     int get_kmin() const { return kmin; }
     /** @brief Get max grid index l (z-dir, incl. buffers). @return lmax. */
     int get_lmax() const { return lmax; }
     /** @brief Get min grid index l (z-dir, incl. buffers). @return lmin. */
     int get_lmin() const { return lmin; }
 
     /** @brief Get upper evolution index j (x-dir). @return jui. */
     int get_jui() const { return jui; }
     /** @brief Get lower evolution index j (x-dir). @return jli. */
     int get_jli() const { return jli; }
     /** @brief Get upper evolution index k (y-dir). @return kui. */
     int get_kui() const { return kui; }
     /** @brief Get lower evolution index k (y-dir). @return kli. */
     int get_kli() const { return kli; }
     /** @brief Get upper evolution index l (z-dir). @return lui. */
     int get_lui() const { return lui; }
     /** @brief Get lower evolution index l (z-dir). @return lli. */
     int get_lli() const { return lli; }
 
     /** @brief Get total number of grid points in x (incl. buffers). @return nx. */
     int get_nx() const { return nx; }
     /** @brief Get total number of grid points in y (incl. buffers). @return ny. */
     int get_ny() const { return ny; }
     /** @brief Get total number of grid points in z (incl. buffers). @return nz. */
     int get_nz() const { return nz; }
 
     /** @brief Get grid index l of max Hamiltonian constraint violation. @return lhm. */
     int get_lhm() const { return lhm; }
     /** @brief Get grid index l of max Momentum constraint violation. @return lmm. */
     int get_lmm() const { return lmm; }
     /** @brief Get grid index l of max Kretschmann invariant. @return lkm. */
     int get_lkm() const { return lkm; }
     /** @brief Get grid index l of max Weyl invariant. @return lwm. */
     int get_lwm() const { return lwm; }
     /** @brief Get excision grid index offset. @return exg. */
     int get_exg() const { return exg; }
     /** @brief Get grid index l just outside the horizon. @return hl. */
     int get_hl() const { return hl; }
 
     /** @brief Get the layer number (0 for base). @return layn. */
     int get_layn() const { return layn; }
     /** @brief Get the minimum l index for constraint checks. @return llmin. */
     int get_llmin() const { return llmin; }
 
     /** @brief Get current simulation time. @return t. */
     double get_t() const { return t; }
     /** @brief Get current time step dt for this layer. @return dt. */
     double get_dt() const { return dt; }
     /** @brief Get base time step dt0 (usually coarsest layer). @return dt0. */
     double get_dt0() const { return dt0; }
     /** @brief Get previous time step dtp. @return dtp. */
     double get_dtp() const { return dtp; }
     /** @brief Get time step before previous dtpp. @return dtpp. */
     double get_dtpp() const { return dtpp; }
     /** @brief Get maximum simulation time. @return tmax. */
     double get_tmax() const { return tmax; }
     /** @brief Get CFL factor. @return cfl. */
     double get_cfl() const { return cfl; }
     /** @brief Get fluid EOS parameter w. @return fluidw. */
     double get_fluidw() const { return fluidw; }
     /** @brief Get scalar field mass parameter m. @return scalarm. */
     double get_scalarm() const { return scalarm; }
 
     /** @brief Get grid spacing dx. @return dx. */
     double get_dx() const { return dx; }
     /** @brief Get grid spacing dy. @return dy. */
     double get_dy() const { return dy; }
     /** @brief Get grid spacing dz. @return dz. */
     double get_dz() const { return dz; }
     /** @brief Get cell volume dvol. @return dvol. */
     double get_dvol() const { return dvol; }
     /** @brief Get inverse grid spacing 1/dx. @return dxi. */
     double get_dxi() const { return dxi; }
     /** @brief Get inverse grid spacing 1/dy. @return dyi. */
     double get_dyi() const { return dyi; }
     /** @brief Get inverse grid spacing 1/dz. @return dzi. */
     double get_dzi() const { return dzi; }
     /** @brief Get 0.5/dx. @return dxi2. */
     double get_dxi2() const { return dxi2; }
     /** @brief Get 0.5/dy. @return dyi2. */
     double get_dyi2() const { return dyi2; }
     /** @brief Get 0.5/dz. @return dzi2. */
     double get_dzi2() const { return dzi2; }
     /** @brief Get 0.25/dx (?). @return dxi4. */
     double get_dxi4() const { return dxi4; }
     /** @brief Get 0.25/dy (?). @return dyi4. */
     double get_dyi4() const { return dyi4; }
     /** @brief Get 0.25/dz (?). @return dzi4. */
     double get_dzi4() const { return dzi4; }
 
     /** @brief Get max physical coordinate x. @return xu. */
     double get_xu() const { return xu; }
     /** @brief Get min physical coordinate x. @return xl. */
     double get_xl() const { return xl; }
     /** @brief Get max physical coordinate y. @return yu. */
     double get_yu() const { return yu; }
     /** @brief Get min physical coordinate y. @return yl. */
     double get_yl() const { return yl; }
     /** @brief Get max physical coordinate z. @return zu. */
     double get_zu() const { return zu; }
     /** @brief Get min physical coordinate z. @return zl. */
     double get_zl() const { return zl; }
 
     /** @brief Get average Hamiltonian constraint violation. @return ham. */
     double get_ham() const { return ham; }
     /** @brief Get max Hamiltonian constraint violation. @return hammax. */
     double get_hammax() const { return hammax; }
     /** @brief Get max Kretschmann invariant value. @return Kremax. */
     double get_Kremax() const { return Kremax; }
     /** @brief Get max Weyl invariant value. @return Weylmax. */
     double get_Weylmax() const { return Weylmax; }
     /** @brief Get average Momentum constraint violation (z). @return mom. */
     double get_mom() const { return mom; }
     /** @brief Get max Momentum constraint violation (z). @return mommax. */
     double get_mommax() const { return mommax; }
 
     /** @brief Get gauge parameter etaa. @return etaa. */
     double get_etaa() const { return etaa; }
     /** @brief Get gauge parameter etab. @return etab. */
     double get_etab() const { return etab; }
     /** @brief Get gauge parameter etabb. @return etabb. */
     double get_etabb() const { return etabb; }
     /** @brief Get cosmological constant lambda. @return lambda. */
     double get_lambda() const { return lambda; }
     /** @brief Get initial time tini. @return tini. */
     double get_tini() const { return tini; }
     /** @brief Get Kreiss-Oliger dissipation coefficient. @return KOep. */
     double get_KOep() const { return KOep; }
     /** @brief Get MUSCL kappa parameter. @return kap_MUSCL. */
     double get_Mkap() const { return kap_MUSCL; }
     /** @brief Get minmod b parameter. @return b_minmod. */
     double get_b() const { return b_minmod; }
     /** @brief Get initial Hubble parameter Hb. @return Hb. */
     double get_Hb() const { return Hb; }
 
     /**
      * @brief Get boundary flag at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Value of the boundary flag.
      */
     int get_bflag(int l,int k,int j) const { return bflag[l-lmin][k-kmin][j-jmin]; }
     /**
      * @brief Get horizon flag at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Value of the horizon flag.
      */
     int get_hflag(int l,int k,int j) const { return hflag[l-lmin][k-kmin][j-jmin]; }
 
     /** @brief Get physical x-coordinate at index j. @param j x-index. @return x[j-jmin]. */
     double get_x(int j) const { return x[j-jmin]; }
     /** @brief Get physical y-coordinate at index k. @param k y-index. @return y[k-kmin]. */
     double get_y(int k) const { return y[k-kmin]; }
     /** @brief Get physical z-coordinate at index l. @param l z-index. @return z[l-lmin]. */
     double get_z(int l) const { return z[l-lmin]; }
     /** @brief Calculate physical x-coordinate for any integer index j. @param j x-index. @return xl + dx*(j-jl[0]). */
     double get_ext_x(int j) const { return(xl +dx*(double(j-jl[0]))); }
     /** @brief Calculate physical y-coordinate for any integer index k. @param k y-index. @return yl + dy*(k-kl[0]). */
     double get_ext_y(int k) const { return(yl +dy*(double(k-kl[0]))); }
     /** @brief Calculate physical z-coordinate for any integer index l. @param l z-index. @return zl + dz*(l-ll[0]). */
     double get_ext_z(int l) const { return(zl +dz*(double(l-ll[0]))); }
 
     /**
      * @brief Get value of state variable `i` at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index (see bv documentation).
      * @return Value of bv[i][l-lmin][k-kmin][j-jmin].
      */
     double get_bv(int l,int k,int j,int i) const { return bv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of RHS variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Value of dbv[i][l-lmin][k-kmin][j-jmin]. */
     double get_dbv(int l,int k,int j,int i) const { return dbv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of previous state variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Value of bv0[i][l-lmin][k-kmin][j-jmin]. */
     double get_bv0(int l,int k,int j,int i) const { return bv0[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of state variable `i` two steps ago at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Value of bv1[i][l-lmin][k-kmin][j-jmin]. */
     double get_bv1(int l,int k,int j,int i) const { return bv1[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of Runge-Kutta sum variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Value of bvr[i][l-lmin][k-kmin][j-jmin]. */
     double get_bvr(int l,int k,int j,int i) const { return bvr[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of temporary interpolation variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Interpolation variable index. @return Value of tmpipol[i][l-lmin][k-kmin][j-jmin]. */
     double get_tmpipol(int l,int k,int j,int i) const { return tmpipol[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of constraint variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Constraint index (see con documentation). @return Value of con[i][l-lmin][k-kmin][j-jmin]. */
     double get_con(int l,int k,int j,int i) const { return con[i][l-lmin][k-kmin][j-jmin]; }
 
     /** @brief Get (dZ/dz)^2 at index l. @param l z-index. @return flat_df2z[l-lmin]. */
     double get_flat_df2z(int l) const { return flat_df2z[l-lmin]; }
     /** @brief Get computational coordinate Z = funcf(z) at index l. @param l z-index. @return coordZ[l-lmin]. */
     double get_coordZ(int l) const { return coordZ[l-lmin]; }
     /** @brief Get bar{Gamma}^z_{zz} at index l. @param l z-index. @return flat_Gamz[l-lmin]. */
     double get_flat_Gamz(int l) const { return flat_Gamz[l-lmin]; }
     /** @brief Get d/dz (bar{Gamma}^z_{zz}) at index l. @param l z-index. @return flat_dGamz[l-lmin]. */
     double get_flat_dGamz(int l) const { return flat_dGamz[l-lmin]; }
 
     /** @brief Get spline matrix lower diagonal element (even) at index l. @param l z-index. @return alp_e[l-lmin]. */
     double get_alp_e(int l) const { return alp_e[l-lmin]; }
     /** @brief Get spline matrix upper diagonal element (even) at index l. @param l z-index. @return bet_e[l-lmin]. */
     double get_bet_e(int l) const { return bet_e[l-lmin]; }
     /** @brief Get spline matrix diagonal element (even) at index l. @param l z-index. @return gam_e[l-lmin]. */
     double get_gam_e(int l) const { return gam_e[l-lmin]; }
     /** @brief Get spline matrix lower diagonal element (odd) at index l. @param l z-index. @return alp_o[l-lmin]. */
     double get_alp_o(int l) const { return alp_o[l-lmin]; }
     /** @brief Get spline matrix upper diagonal element (odd) at index l. @param l z-index. @return bet_o[l-lmin]. */
     double get_bet_o(int l) const { return bet_o[l-lmin]; }
     /** @brief Get spline matrix diagonal element (odd) at index l. @param l z-index. @return gam_o[l-lmin]. */
     double get_gam_o(int l) const { return gam_o[l-lmin]; }
     /** @brief Get spline second derivative for variable `i` at index l. @param l z-index. @param i Variable index. @return ddy[i][l-lmin]. */
     double get_ddy(int l,int i) const { return ddy[i][l-lmin]; }
 
     /** @brief Get value of output variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Output variable index (see outv documentation). @return Value of outv[i][l-lmin][k-kmin][j-jmin]. */
     double get_outv(int l,int k,int j,int i) const { return outv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference z-coordinate for grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @return Value of refz[l-lmin][k-kmin][j-jmin]. */
     double get_refz(int l,int k,int j) const { return refz[l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of primitive fluid variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Primitive variable index (see primv documentation). @return Value of primv[i][l-lmin][k-kmin][j-jmin]. */
     double get_primv(int l,int k,int j,int i) const { return primv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of fluid x-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Value of flux_x[i][l-lmin][k-kmin][j-jmin]. */
     double get_flux_x(int l,int k,int j,int i) const { return flux_x[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of fluid y-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Value of flux_y[i][l-lmin][k-kmin][j-jmin]. */
     double get_flux_y(int l,int k,int j,int i) const { return flux_y[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get value of fluid z-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Value of flux_z[i][l-lmin][k-kmin][j-jmin]. */
     double get_flux_z(int l,int k,int j,int i) const { return flux_z[i][l-lmin][k-kmin][j-jmin]; }
 
     /** @brief Get pointer to the 3D grid data for state variable `i`. @param i Variable index. @return Pointer bv[i]. */
     double*** get_bv(int i) const { return bv[i]; }
     /** @brief Get pointer to the 3D grid data for RHS variable `i`. @param i Variable index. @return Pointer dbv[i]. */
     double*** get_dbv(int i) const { return dbv[i]; }
 
     // --- GET derivative functions ---
     // Calculate finite differences using 4th-order centered stencils.
 
     /** @brief Calculate 4th-order centered 1st derivative d(bv[i])/dx at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_x(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered 1st derivative d(bv[i])/dy at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_y(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered 1st derivative d(bv[i])/dz at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_z(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered 2nd derivative d^2(bv[i])/dx^2 at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_xx(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered 2nd derivative d^2(bv[i])/dy^2 at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_yy(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered 2nd derivative d^2(bv[i])/dz^2 at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_zz(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered mixed derivative d^2(bv[i])/dxdy at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_xy(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered mixed derivative d^2(bv[i])/dxdz at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_xz(int l,int k,int j,int i) const;
     /** @brief Calculate 4th-order centered mixed derivative d^2(bv[i])/dydz at (l,k,j). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Derivative value. */
     double get_f_yz(int l,int k,int j,int i) const;
 
     /** @brief Interpolate bv[i] to the lower mid-point in x (j-1/2). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Interpolated value. */
     double get_ipol_x_lower_mid(int l,int k,int j,int i) const;
     /** @brief Interpolate bv[i] to the lower mid-point in y (k-1/2). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Interpolated value. */
     double get_ipol_y_lower_mid(int l,int k,int j,int i) const;
     /** @brief Interpolate bv[i] to the lower mid-point in z (l-1/2). @param l z-idx. @param k y-idx. @param j x-idx. @param i var-idx. @return Interpolated value. */
     double get_ipol_z_lower_mid(int l,int k,int j,int i) const;
 
     // --- SET functions ---
     // Provide write access to protected members.
 
     /** @brief Set fluid evolution flag. @param f New flag value. */
     void set_fluidevo(bool f) { fluidevo=f; }
     /** @brief Set scalar evolution flag. @param s New flag value. */
     void set_scalarevo(bool s) { scalarevo=s; }
     /** @brief Set mesh refinement flag. @param s New flag value. */
     void set_mrf(bool s) { mrf=s; }
     /** @brief Set current simulation time. @param time New time value. */
     void set_t(double time) { t=time; }
     /** @brief Set current time step dt for this layer. @param time New dt value. */
     void set_dt(double time) { dt=time; }
     /** @brief Set base time step dt0. @param time New dt0 value. */
     void set_dt0(double time) { dt0=time; }
     /** @brief Set previous time step dtp. @param time New dtp value. */
     void set_dtp(double time) { dtp=time; }
     /** @brief Set time step before previous dtpp. @param time New dtpp value. */
     void set_dtpp(double time) { dtpp=time; }
     /** @brief Set maximum simulation time tmax. @param time New tmax value. */
     void set_tmax(double time) { tmax=time; }
     /** @brief Set CFL factor. @param c New CFL value. */
     void set_cfl(double c) { cfl=c; }
     /** @brief Set gauge parameter etaa. @param e New etaa value. */
     void set_etaa(double e) { etaa=e; }
     /** @brief Set gauge parameter etab. @param e New etab value. */
     void set_etab(double e) { etab=e; }
     /** @brief Set gauge parameter etabb. @param e New etabb value. */
     void set_etabb(double e) { etabb=e; }
     /** @brief Set cosmological constant lambda. @param l New lambda value. */
     void set_lambda(double l) { lambda=l; }
     /** @brief Set initial Hubble parameter Hb. @param hb New Hb value. */
     void set_Hb(double hb) { Hb=hb; }
     /** @brief Set initial time tini. @param t New tini value. */
     void set_tini(double t) { tini=t; }
     /** @brief Set Kreiss-Oliger dissipation coefficient. @param l New KOep value. */
     void set_KOep(double l) { KOep=l; }
     /** @brief Set excision grid index offset. @param eg New exg value. */
     void set_exg(int eg) { exg=eg; }
     /** @brief Set grid index l just outside horizon. @param h New hl value. */
     void set_hl(int h) { hl=h; }
     /** @brief Set inhomogeneous grid amplitude. @param a New amp value. */
     void set_amp(double a) { amp=a; }
     /** @brief Set fluid EOS parameter w. @param fw New fluidw value. */
     void set_fluidw(double fw) { fluidw=fw; }
     /** @brief Set scalar field mass parameter m. @param sm New scalarm value. */
     void set_scalarm(double sm) { scalarm=sm; }
     /** @brief Set MUSCL kappa parameter. @param k New kap_MUSCL value. */
     void set_Mkap(double k) { kap_MUSCL=k; }
     /** @brief Set minmod b parameter. @param b New b_minmod value. */
     void set_b(double b) { b_minmod=b; }
     /** @brief Set excision flag. @param e New flag value. */
     void set_exc(bool e) { exc=e; }
     /** @brief Set minimum l index for constraint checks. @param lll New llmin value. */
     void set_llmin(int lll) { llmin=lll; }
 
     /**
      * @brief Get reference to boundary flag at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to bflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_bflag(int l,int k,int j) { return bflag[l-lmin][k-kmin][j-jmin]; }
     /**
      * @brief Get reference to horizon flag at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @return Reference to hflag[l-lmin][k-kmin][j-jmin].
      */
     int& set_hflag(int l,int k,int j) { return hflag[l-lmin][k-kmin][j-jmin]; }
 
     /**
      * @brief Get reference to state variable `i` at grid point (l,k,j).
      * @param l z-index.
      * @param k y-index.
      * @param j x-index.
      * @param i Variable index.
      * @return Reference to bv[i][l-lmin][k-kmin][j-jmin].
      */
     double& set_bv(int l,int k,int j,int i) { return bv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to RHS variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Reference to dbv[i][l-lmin][k-kmin][j-jmin]. */
     double& set_dbv(int l,int k,int j,int i) { return dbv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to previous state variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Reference to bv0[i][l-lmin][k-kmin][j-jmin]. */
     double& set_bv0(int l,int k,int j,int i) { return bv0[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to state variable `i` two steps ago at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Reference to bv1[i][l-lmin][k-kmin][j-jmin]. */
     double& set_bv1(int l,int k,int j,int i) { return bv1[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to Runge-Kutta sum variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Variable index. @return Reference to bvr[i][l-lmin][k-kmin][j-jmin]. */
     double& set_bvr(int l,int k,int j,int i) { return bvr[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to temporary interpolation variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Interpolation variable index. @return Reference to tmpipol[i][l-lmin][k-kmin][j-jmin]. */
     double& set_tmpipol(int l,int k,int j,int i) { return tmpipol[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to reference z-coordinate for grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @return Reference to refz[l-lmin][k-kmin][j-jmin]. */
     double& set_refz(int l,int k,int j) { return refz[l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to constraint variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Constraint index. @return Reference to con[i][l-lmin][k-kmin][j-jmin]. */
     double& set_con(int l,int k,int j,int i) { return con[i][l-lmin][k-kmin][j-jmin]; }
 
     /** @brief Get reference to computational coordinate Z at index l. @param l z-index. @return Reference to coordZ[l-lmin]. */
     double& set_coordZ(int l) { return coordZ[l-lmin]; }
     /** @brief Get reference to (dZ/dz)^2 at index l. @param l z-index. @return Reference to flat_df2z[l-lmin]. */
     double& set_flat_df2z(int l) { return flat_df2z[l-lmin]; }
     /** @brief Get reference to bar{Gamma}^z_{zz} at index l. @param l z-index. @return Reference to flat_Gamz[l-lmin]. */
     double& set_flat_Gamz(int l) { return flat_Gamz[l-lmin]; }
     /** @brief Get reference to d/dz (bar{Gamma}^z_{zz}) at index l. @param l z-index. @return Reference to flat_dGamz[l-lmin]. */
     double& set_flat_dGamz(int l) { return flat_dGamz[l-lmin]; }
 
     /** @brief Get reference to spline matrix lower diagonal element (even) at index l. @param l z-index. @return Reference to alp_e[l-lmin]. */
     double& set_alp_e(int l) { return alp_e[l-lmin]; }
     /** @brief Get reference to spline matrix upper diagonal element (even) at index l. @param l z-index. @return Reference to bet_e[l-lmin]. */
     double& set_bet_e(int l) { return bet_e[l-lmin]; }
     /** @brief Get reference to spline matrix diagonal element (even) at index l. @param l z-index. @return Reference to gam_e[l-lmin]. */
     double& set_gam_e(int l) { return gam_e[l-lmin]; }
     /** @brief Get reference to spline matrix lower diagonal element (odd) at index l. @param l z-index. @return Reference to alp_o[l-lmin]. */
     double& set_alp_o(int l) { return alp_o[l-lmin]; }
     /** @brief Get reference to spline matrix upper diagonal element (odd) at index l. @param l z-index. @return Reference to bet_o[l-lmin]. */
     double& set_bet_o(int l) { return bet_o[l-lmin]; }
     /** @brief Get reference to spline matrix diagonal element (odd) at index l. @param l z-index. @return Reference to gam_o[l-lmin]. */
     double& set_gam_o(int l) { return gam_o[l-lmin]; }
     /** @brief Get reference to spline second derivative for variable `i` at index l. @param l z-index. @param i Variable index. @return Reference to ddy[i][l-lmin]. */
     double& set_ddy(int l,int i) { return ddy[i][l-lmin]; }
 
     /** @brief Get reference to output variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Output variable index. @return Reference to outv[i][l-lmin][k-kmin][j-jmin]. */
     double& set_outv(int l,int k,int j,int i) { return outv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to primitive fluid variable `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Primitive variable index. @return Reference to primv[i][l-lmin][k-kmin][j-jmin]. */
     double& set_primv(int l,int k,int j,int i) { return primv[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to fluid x-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Reference to flux_x[i][l-lmin][k-kmin][j-jmin]. */
     double& set_flux_x(int l,int k,int j,int i) { return flux_x[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to fluid y-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Reference to flux_y[i][l-lmin][k-kmin][j-jmin]. */
     double& set_flux_y(int l,int k,int j,int i) { return flux_y[i][l-lmin][k-kmin][j-jmin]; }
     /** @brief Get reference to fluid z-flux `i` at grid point (l,k,j). @param l z-index. @param k y-index. @param j x-index. @param i Flux index. @return Reference to flux_z[i][l-lmin][k-kmin][j-jmin]. */
     double& set_flux_z(int l,int k,int j,int i) { return flux_z[i][l-lmin][k-kmin][j-jmin]; }
 
     /** @brief Get pointer to the 3D grid data for state variable `i`. @param i Variable index. @return Pointer bv[i]. */
     double*** set_bv(int i) const { return bv[i]; }
     /** @brief Get pointer to the 3D grid data for RHS variable `i`. @param i Variable index. @return Pointer dbv[i]. */
     double*** set_dbv(int i) const { return dbv[i]; }
     /** @brief Get pointer to the 3D grid data for previous state variable `i`. @param i Variable index. @return Pointer bv0[i]. */
     double*** set_bv0(int i) const { return bv0[i]; }
     /** @brief Get pointer to the 3D grid data for state variable `i` two steps ago. @param i Variable index. @return Pointer bv1[i]. */
     double*** set_bv1(int i) const { return bv1[i]; }
     /** @brief Get pointer to the 3D grid data for Runge-Kutta sum variable `i`. @param i Variable index. @return Pointer bvr[i]. */
     double*** set_bvr(int i) const { return bvr[i]; }
     /** @brief Get pointer to the 3D grid data for temporary interpolation variable `i`. @param i Interpolation variable index. @return Pointer tmpipol[i]. */
     double*** set_tmpipol(int i) const { return tmpipol[i]; }
 
     // --- SET zero or unity functions ---
 
     /** @brief Set all major grid arrays (bv, dbv, bv0, bvr, tmpipol, con, outv, primv, flux_*) and related 1D arrays to zero. */
     void set_zero_all();
     /** @brief Set primitive fluid variable array `primv` to zero. */
     void set_zero_primv();
     /** @brief Set current state variable array `bv` to zero. */
     void set_zero();
     /** @brief Set previous state variable array `bv0` to zero. */
     void set_zero_0();
     /** @brief Set state variable array `bv1` (two steps ago) to zero. */
     void set_zero_1();
     /** @brief Set RHS array `dbv` to zero (respecting `nonzero` flags). */
     void set_zero_d();
     /** @brief Set RHS array `dbv` to zero within the evolution domain (respecting `nonzero` flags). */
     void set_zero_d_exc();
     /** @brief Set Runge-Kutta sum array `bvr` to zero. */
     void set_zero_r();
     /** @brief Set boundary flag array `bflag` to zero. */
     void set_bflag_zero();
     /** @brief Set horizon flag array `hflag` to zero. */
     void set_hflag_zero();
 
     // --- UPDATE functions ---
 
     /** @brief Copy current state `bv` to previous state `bv0`. */
     void setv0();
     /** @brief Copy previous state `bv0` to state two steps ago `bv1`. */
     void set01();
     /** @brief Copy previous state `bv0` back to current state `bv`. (Used for RK steps). */
     void set0v();
     /**
      * @brief Add contribution to Runge-Kutta sum: bvr += dbv * dt.
      * @param dt Time step factor for this contribution.
      */
     void runge_kutta(double dt);
     /**
      * @brief Update state using Euler step: bv = bv0 + dbv * dt.
      * @param dt Time step size.
      */
     void new_bv(double dt);
     /** @brief Final update step for 4th order RK: bv = bv0 + bvr. */
     void new_bv4();
 
     // --- Check functions ---
 
     /**
      * @brief Calculate Hamiltonian and Momentum constraint violations.
      * @details Calculates constraints using `con` array values, computes average and maximum violations,
      *          updates `ham`, `hammax`, `mom`, `mommax`, `dGam`, `dGammax`, `lhm`, `lmm`.
      *          Excludes excised points and boundary points (`negb`) from averages.
      */
     void check_const();
     /**
      * @brief Calculate maximum Kretschmann invariant value.
      * @details Finds the maximum absolute value in `outv[0]` (excluding excised points)
      *          and updates `Kremax` and `lkm`.
      */
     void check_Kremax();
     /**
      * @brief Calculate maximum Weyl invariant value (?).
      * @details Finds the maximum absolute value in `outv[3]` (excluding excised points)
      *          and updates `Weylmax` and `lwm`.
      */
     void check_Weylmax();
     /**
      * @brief Check for apparent horizon formation and output data.
      * @param fout Output file stream for horizon data (`out_horizon_XX.dat`).
      * @details Searches for horizons based on outgoing null expansion (`outv[4]`) and compactness (`outv[9]/outv[5]`).
      *          Updates `hform` flag, potentially activates excision (`exc`), updates `horis` array,
      *          and writes horizon coordinate and areal radius to `fout`. Tracks multiple horizons.
      */
     void check_horizon(ofstream& fout);
     /**
      * @brief Check for necks (minima of areal radius) and output data.
      * @param fout Output file stream for neck data (`out_neck_XX.dat`).
      * @details Searches for minima based on `drR = outv[11]`. Updates `neck` array and writes neck coordinates to `fout`.
      *          Tracks multiple necks.
      */
     void check_neck(ofstream& fout);
 
     // --- Flag setting functions ---
 
     /**
      * @brief Set boundary flags `bflag` for excision.
      * @details Marks points inside `lli + exg` as excised (-1) and sets flags 1-6
      *          in the buffer region just outside the excision boundary.
      */
     void set_excflags();
 
     // --- Fluid functions (defined in cosmos_s_fluid.cpp) ---
 
     /**
      * @brief Calculate fluid pressure P from rest mass density rho.
      * @param rho Rest mass density.
      * @return Pressure P (using `fluidw`).
      */
     double pres(double rho);
     /**
      * @brief Calculate dP/drho.
      * @param rho Rest mass density.
      * @return dP/drho (which is `fluidw`).
      */
     double dpres(double rho);
     /**
      * @brief Calculate primitive rho and Lorentz factor Gamma from conserved E and momentum squared p2.
      * @param Ene Conserved energy density E.
      * @param p2 Squared momentum p^mu p_mu.
      * @param[out] rho Calculated primitive rest mass density.
      * @param[out] Gam Calculated Lorentz factor Gamma = -u^mu n_mu.
      */
     void get_rhoGam(double Ene, double p2,double& rho,double& Gam);
     /**
      * @brief Minmod slope limiter function.
      * @param a First slope estimate.
      * @param b Second slope estimate.
      * @return Limited slope using `b_minmod`.
      */
     double minmod(double a,double b);
     /**
      * @brief Sign function.
      * @param A Input value.
      * @return 1 if A>0, -1 if A<0, 0 if A=0.
      */
     double sign(double A);
     /**
      * @brief Convert conserved fluid variables (E, S_i, D) to primitive variables (rho, V^i, epsilon).
      * @details Reads conserved variables from `bv` array, calculates primitives, and stores them in `primv` array.
      */
     void dyntoprim();
 
     // --- Interpolation functions (defined in cosmos_s_ipol.cpp) ---
 
     /**
      * @brief Lagrange polynomial interpolation.
      * @param rr Target coordinate value.
      * @param xx Array of source grid point coordinates.
      * @param yy Array of source grid point variable values.
      * @param order Order of the interpolation polynomial.
      * @return Interpolated value at `rr`.
      */
     double ipol( double rr,double *xx,double *yy,int order );
     /**
      * @brief Interpolate state variable `bv[number]` to coordinate `zc` using Lagrange interpolation.
      * @param lc Lower grid index l near `zc`.
      * @param zc Target z-coordinate.
      * @param order Interpolation order.
      * @param number Variable index `i` for `bv`.
      * @return Interpolated value.
      * @note Assumes interpolation along z-axis (j=0, k=0).
      */
     double bv_ipol(int lc,double zc,int order,int number);
     /**
      * @brief Interpolate previous state variable `bv0[number]` to coordinate `zc` using Lagrange interpolation.
      * @param lc Lower grid index l near `zc`.
      * @param zc Target z-coordinate.
      * @param order Interpolation order.
      * @param number Variable index `i` for `bv0`.
      * @return Interpolated value.
      * @note Assumes interpolation along z-axis (j=0, k=0).
      */
     double bv0_ipol(int lc,double zc,int order,int number);
     /**
      * @brief Interpolate state variable `bv1[number]` (two steps ago) to coordinate `zc` using Lagrange interpolation.
      * @param lc Lower grid index l near `zc`.
      * @param zc Target z-coordinate.
      * @param order Interpolation order.
      * @param number Variable index `i` for `bv1`.
      * @return Interpolated value.
      * @note Assumes interpolation along z-axis (j=0, k=0).
      */
     double bv1_ipol(int lc,double zc,int order,int number);
     /**
      * @brief Calculate second derivatives `ddy[i]` for cubic spline interpolation along z-axis.
      * @param i Variable index `i` for `bv`.
      * @details Solves the tridiagonal system for the second derivatives using precomputed LU decomposition factors.
      *          Handles even/odd parity at the center.
      */
     void put_ddy(int i);
     /**
      * @brief Interpolate state variable `bv[number]` to coordinate `zc` using cubic spline interpolation.
      * @param lc Lower grid index l near `zc`.
      * @param zc Target z-coordinate.
      * @param number Variable index `i` for `bv`.
      * @return Interpolated value.
      * @note Requires `put_ddy(number)` to be called first. Assumes interpolation along z-axis (j=0, k=0).
      */
     double bv_spline_ipol(int lc,double zc,int number);
     /**
      * @brief Precompute LU decomposition factors for the spline tridiagonal matrix.
      * @details Calculates and stores `alp_e`, `bet_e`, `gam_e`, `alp_o`, `bet_o`, `gam_o`.
      */
     void set_luvec();
 
     // --- BSSN functions (defined in cosmos_s_bssn.cpp) ---
 
     /** @brief Calculate advection terms for BSSN variables using 4th order upwind/downwind stencils based on shift vector sign. */
     void BSSN_adv();
     /**
      * @brief Calculate the RHS for BSSN and gauge variables, and optionally fluid/scalar fields.
      * @param itype Stage number within the Runge-Kutta integration (1 to 4).
      * @details This is the core evolution function. It calculates spatial derivatives, computes source terms,
      *          and adds contributions to the `dbv` array. Includes fluid source terms and scalar field evolution if enabled.
      *          Also calculates constraints and curvature invariants if `itype == 1`.
      */
     void BSSN(int itype);
     /**
      * @brief Enforce algebraic constraints: det(tilde{gamma}) = 1 and tr(tilde{A}) = 0.
      * @details Rescales the conformal metric components and trace-free extrinsic curvature components.
      *          Adjusts `wa` (conformal factor) and `ek` (tr K) accordingly.
      */
     void enforce_const();
     /**
      * @brief Add Kreiss-Oliger dissipation to the RHS (`dbv`).
      * @param lup Upper grid index l up to which dissipation is applied.
      * @details Applies 6th-order dissipation based on the `KOep` parameter.
      */
     void KOdiss(int lup);
     /**
      * @brief Calculate numerical fluxes for fluid evolution using MUSCL-Hancock scheme and HLL solver.
      * @details Computes fluxes at cell interfaces (l-1/2) and stores them in `flux_x`, `flux_y`, `flux_z`.
      *          Uses primitive variables from `primv` and MUSCL reconstruction with minmod limiter.
      */
     void flux_fill();
     /**
      * @brief Apply simple extrapolation for RHS in the excision region buffer.
      * @details Copies `dbv` values from `lli+exg` to points just inside (`lli+exg-1` down to `lli+exg-6`).
      */
     void excision();
     /**
      * @brief Fill state variables `bv` inside the excision region using quadratic extrapolation.
      * @details Extrapolates values from `lli+exg-6` inwards based on values at `lli+exg-6`, `lli+exg-5`, `lli+exg-4`.
      *          Respects even/odd parity.
      */
     void excision_quadfunc();
 
     // --- Initial data functions (defined in cosmos_s_initial.cpp) ---
 
     /**
      * @brief Read simulation state from a checkpoint file.
      * @param fcontinue Input file stream opened to the checkpoint file (`out_all.dat`).
      * @details Reads time, dtp, dtpp, and all state variables (`bv`, `bv0`, `bv1`) from the file.
      */
     void initial_continue(ifstream& fcontinue);
     /**
      * @brief Inhomogeneous grid mapping function Z = f(z).
      * @param Z Computational coordinate (z in this code).
      * @return Physical coordinate (Z in this code, often denoted r).
      */
     double funcf(double Z);
     /**
      * @brief Inverse inhomogeneous grid mapping function z = f^{-1}(Z).
      * @param X Physical coordinate (Z in this code, often denoted r).
      * @return Computational coordinate (z in this code).
      * @note Uses an iterative solver.
      */
     double ifuncf(double X);
     /**
      * @brief First derivative of the grid mapping function dZ/dz.
      * @param Z Computational coordinate (z in this code).
      * @return Derivative value df/dz.
      */
     double df(double Z);
     /**
      * @brief Second derivative of the grid mapping function d^2Z/dz^2.
      * @param Z Computational coordinate (z in this code).
      * @return Second derivative value d^2f/dz^2.
      */
     double ddf(double Z);
     /**
      * @brief Third derivative of the grid mapping function d^3Z/dz^3.
      * @param Z Computational coordinate (z in this code).
      * @return Third derivative value d^3f/dz^3.
      */
     double dddf(double Z);
     /**
      * @brief Calculate and store flat-space geometric quantities for the inhomogeneous grid.
      * @details Computes `flat_df2z`, `coordZ`, `flat_Gamz`, `flat_dGamz`.
      */
     void set_flat();
     /**
      * @brief Calculate and store reference z-coordinates for off-axis points.
      * @details Computes `refz` array using `ifuncf`.
      */
     void set_refz();
     /**
      * @brief Calculate initial values for BSSN Gamma variables (tilde{Gamma}^i).
      * @details Computes `gamma0_x`, `gamma0_y`, `gamma0_z` based on initial metric derivatives and sets `bv[21]`, `bv[22]`, `bv[23]`.
      */
     void set_Gam();
     /**
      * @brief Set initial energy-momentum tensor components based on constraint equations.
      * @details Calculates Ricci tensor and extrinsic curvature terms to solve Hamiltonian and Momentum constraints
      *          for initial E, S_i, assuming they were initially zero or incorrect. Updates `bv[24]` to `bv[27]`.
      *          Also converts back to primitive variables and updates `bv[28]`.
      */
     void set_enemomini();
     /**
      * @brief Set various simulation parameters.
      * @param cfli CFL factor.
      * @param etaai Gauge parameter etaa.
      * @param etabi Gauge parameter etab.
      * @param etabbi Gauge parameter etabb.
      * @param lambdai Cosmological constant.
      * @param dt0i Initial base time step.
      * @param dtpi Initial previous time step.
      * @param dtppi Initial time step before previous.
      * @param ti Current time.
      * @param tinii Initial time.
      * @param Hbi Initial Hubble parameter.
      * @param KOepi Dissipation coefficient.
      * @param exgi Excision grid offset.
      * @param fluidwi Fluid EOS parameter w.
      * @param scalarmi Scalar field mass m.
      * @param kap_MUSCLi MUSCL kappa parameter.
      * @param b_minmodi Minmod b parameter.
      */
     void initial_params(double cfli,double etaai,double etabi,double etabbi,double lambdai,double dt0i,double dtpi,double dtppi,double ti,double tinii,double Hbi,double KOepi,int exgi,double fluidwi,double scalarmi,double kap_MUSCLi,double b_minmodi);
 
     // --- OUTPUT functions (defined in cosmos_s_output.cpp) ---
 
     /**
      * @brief Print 1D data along the z-axis (at specified j, k) to a file. Uses bv0.
      * @param fn Output file stream (`out_jkz.dat`).
      * @param j x-index to print.
      * @param k y-index to print.
      */
     void print_z(ofstream& fn, int j,int k);
     /**
      * @brief Print 1D data along the z-axis (at specified j, k) to a file. Uses bv.
      * @param fn Output file stream (`out_jkz.dat`).
      * @param j x-index to print.
      * @param k y-index to print.
      */
     void print_bz(ofstream& fn, int j,int k);
     /**
      * @brief Print max Kretschmann and Weyl invariants and their locations to a file.
      * @param fout Output file stream.
      */
     void print_Kremax(ofstream& fout);
     /**
      * @brief Print average and max constraint violations to a file.
      * @param fout Output file stream (`out_const.dat`).
      */
     void print_const(ofstream& fout);
     /**
      * @brief Print all state variables (bv, bv0, bv1) for checkpointing/restarting.
      * @param fout Output file stream (`out_all.dat`).
      * @details Writes time, dtp, dtpp, then dumps bv, bv0, bv1 arrays row by row.
      */
     void print_all(ofstream& fout);
     /**
      * @brief Print 1D data along the z-axis using spline interpolation to half-grid points.
      * @param fout Output file stream.
      */
     void print_ipolz(ofstream& fout);
 
     // --- Potential functions ---
 
     /**
      * @brief Calculate scalar field potential V(phi).
      * @param p Scalar field value phi.
      * @return Potential value 0.5 * (scalarm * p)^2.
      */
     double funcV(double p) { return(0.5*pow(scalarm*p,2)); }
     /**
      * @brief Calculate scalar field potential derivative dV/dphi.
      * @param p Scalar field value phi.
      * @return Derivative value scalarm^2 * p.
      */
     double funcdV(double p) { return(pow(scalarm,2)*p); }
 
 }; // End class Fmv0
 
 /**
  * @class Fmv
  * @brief Derived class representing the base grid (Layer 0).
  * @inherits Fmv0
  *
  * @details
  * Specializes Fmv0 for the base grid. Implements asymptotic boundary conditions
  * and specific initial data generation routines (e.g., long wavelength perturbations).
  */
 class Fmv : public Fmv0{
 
 public:
     /**
      * @brief Constructor for the base grid class Fmv.
      * @param tabs Number of buffer cells in z-direction.
      * @param tabsx Number of buffer cells in x,y-directions.
      * @param lupper Upper grid index (evolution domain) in z-direction.
      * @param xupper Max physical coordinate in x.
      * @param zupper Max physical coordinate in z (evolution domain).
      * @param am Amplitude for inhomogeneous grid mapping.
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature invariant evaluation.
      * @note Sets layer number `layn` to 0.
      */
     Fmv(int tabs,int tabsx,int lupper,
     double xupper,double zupper,double am,bool fld, bool scl, bool cuev) : Fmv0(tabs,tabsx,lupper,
     xupper,zupper,am,fld, scl, cuev){
         layn=0;
     }
 
  public:
 
     // --- BOUNDARY functions (defined in cosmos_s_boundary.cpp) ---
 
     /**
      * @brief Apply Sommerfeld-like outgoing wave boundary condition at outer boundary (z=lu[l]). Static version.
      * @param l Buffer layer index (1 to tab).
      * @param i Variable index.
      * @param bgv Background value for variable `i`.
      * @details Uses interpolation (`bv_ipol`) to estimate value at `z - dr`.
      */
     void asymcond(int l,int i,double bgv);
     /**
      * @brief Apply Sommerfeld-like outgoing wave boundary condition at outer boundary (z=lu[l]). Time-dependent version.
      * @param l Buffer layer index (1 to tab).
      * @param i Variable index.
      * @param bgv1 Background value for variable `i` at previous step (t-dtp).
      * @param bgv Background value for variable `i` at current step (t or intermediate RK step).
      * @param dtime Effective time elapsed since bv1 state (depends on RK stage).
      * @param itype RK stage type (1 to 4).
      * @details Uses interpolation (`bv1_ipol`) to estimate value at `z - dr(t)`.
      */
     void asymcond(int l,int i,double bgv1,double bgv,double dtime,int itype);
     /**
      * @brief Apply asymptotic boundary conditions to all variables in the outer buffer zones.
      * @param itype RK stage type (0 for initial, 1-4 for evolution steps).
      * @details Calls `asymcond` for each variable and buffer layer. Also handles symmetry at the inner boundary (l=ll[l])
      *          and triggers recalculation of interpolated values (`tmpipol`) via spline.
      */
     void boundary_asym(int itype);
 
     // --- Functions for specific initial conditions (defined in cosmos_s_initial.cpp) ---
 
     /**
      * @brief Final setup steps after initial data generation.
      * @details Calls boundary_asym(0), enforce_const(), set_Gam(), set_enemomini(), dyntoprim().
      */
     void set_initial_final();
     /**
      * @brief Load base grid data from a continuation file.
      * @param fcontinue Input file stream for the checkpoint file.
      * @details Calls `initial_continue`, `boundary_asym(0)`, and `dyntoprim`.
      */
     void base_initial_continue(ifstream& fcontinue);
     /**
      * @brief Set initial data for isocurvature long wavelength perturbations.
      * @param mu Amplitude of scalar field perturbation.
      * @param k Wavenumber of scalar field perturbation.
      * @param inr Inner radius of perturbation profile matching region.
      * @param L Outer radius of perturbation profile matching region.
      */
     void initial_iso_longwave(double mu,double k,double inr,double L);
     /**
      * @brief Set initial data for adiabatic long wavelength perturbations.
      * @param mu Amplitude of metric perturbation Psi.
      * @param k Wavenumber of metric perturbation Psi.
      * @param inr Inner radius of perturbation profile matching region.
      * @param L Outer radius of perturbation profile matching region.
      */
     void initial_longwave(double mu,double k,double inr,double L);
     /**
      * @brief Profile function for initial scalar field perturbation Phi(r).
      * @param r Physical radius.
      * @param mu Amplitude.
      * @param kk Wavenumber.
      * @param inr Inner matching radius.
      * @param L Outer matching radius.
      * @return Value of Phi(r).
      */
     double Phi(double r,double mu,double kk, double inr, double L);
     /** @brief Derivative dPhi/dr. @param r Radius. @param mu Amp. @param kk k. @param inr Inner r. @param L Outer r. @return dPhi/dr. */
     double dzPhi(double r,double mu,double kk, double inr, double L);
     /** @brief Second derivative d^2Phi/dr^2. @param r Radius. @param mu Amp. @param kk k. @param inr Inner r. @param L Outer r. @return d^2Phi/dr^2. */
     double ddzPhi(double r,double mu,double kk, double inr, double L);
     /**
      * @brief Profile function for initial metric perturbation Psi(r) = exp(psi(r)).
      * @param r Physical radius.
      * @param mu Amplitude of underlying psi perturbation.
      * @param kk Wavenumber.
      * @param inr Inner matching radius.
      * @param L Outer matching radius.
      * @return Value of Psi(r).
      */
     double Psi(double r,double mu,double kk, double inr, double L);
     /** @brief Derivative dPsi/dr. @param r Radius. @param mu Amp. @param kk k. @param inr Inner r. @param L Outer r. @return dPsi/dr. */
     double dzPsi(double r,double mu,double kk, double inr, double L);
     /** @brief Second derivative d^2Psi/dr^2. @param r Radius. @param mu Amp. @param kk k. @param inr Inner r. @param L Outer r. @return d^2Psi/dr^2. */
     double ddzPsi(double r,double mu,double kk, double inr, double L);
     /** @brief Third derivative d^3Psi/dr^3. @param r Radius. @param mu Amp. @param kk k. @param inr Inner r. @param L Outer r. @return d^3Psi/dr^3. */
     double dddzPsi(double r,double mu,double kk, double inr, double L);
 
 }; // End class Fmv
 
 /**
  * @class Fmv1
  * @brief Derived class representing a refinement grid (Layer > 0).
  * @inherits Fmv0
  *
  * @details
  * Specializes Fmv0 for a refinement level in the FMR hierarchy.
  * Contains pointers to the parent (`llay`) and child (`ulay`) layers.
  * Implements boundary conditions that involve interpolation from the parent layer
  * and restriction to the parent layer. Manages the time stepping synchronization
  * between layers.
  */
 class Fmv1 : public Fmv0{
 private:
     /** @brief Grid index on parent layer defining the boundary for this layer (?). */
     int lb; // TODO: Confirm exact meaning
     /** @brief Pointer to the parent (coarser) layer. */
     Fmv0* llay;
     /** @brief Pointer to the child (finer) layer, if it exists. */
     Fmv1* ulay;
 
 public:
     /**
      * @brief Constructor for a refinement layer class Fmv1.
      * @param tabs Number of buffer cells in z-direction.
      * @param tabsx Number of buffer cells in x,y-directions.
      * @param lupper Upper grid index (evolution domain) in z-direction for this layer.
      * @param xupper Max physical coordinate in x for this layer.
      * @param zupper Max physical coordinate in z (evolution domain) for this layer.
      * @param am Amplitude for inhomogeneous grid mapping.
      * @param fld Boolean flag to enable fluid evolution.
      * @param scl Boolean flag to enable scalar field evolution.
      * @param cuev Boolean flag to enable curvature invariant evaluation.
      * @param lolay Pointer to the parent (lower) layer object.
      * @note Sets layer number `layn` based on parent.
      */
     Fmv1(int tabs,int tabsx, int lupper,double xupper,double zupper,double am,bool fld, bool scl, bool cuev, Fmv0* lolay)
      : Fmv0(tabs,tabsx, lupper,xupper,zupper,am,fld, scl, cuev)
     {
         layn=lolay->get_layn()+1;
         llay=lolay;
         ulay=NULL; // Initialize ulay to NULL
         mrf=false; // Initially no child layer
         cout << "Fmv1 for layer #" << layn << " constructer done" << endl;
     }
 
 public:
     /** @brief Get the boundary index `lb`. @return lb. */
     int get_lb() const { return lb; }
 
     /** @brief Set the boundary index `lb`. @param p New value for lb. */
     void set_lb(int p) { lb=p; }
     /** @brief Set the pointer to the child (upper) layer. @param uplay Pointer to the child Fmv1 object. */
     void set_ulay(Fmv1* uplay) { ulay=uplay; }
 
     // --- FMR specific functions (defined in cosmos_s_fmr.cpp) ---
 
     /**
      * @brief Apply boundary conditions specific to FMR layers.
      * @details Handles symmetry at the inner boundary (l=ll[l]). Outer boundary values are set by `set_boundary`.
      *          Triggers recalculation of interpolated values (`tmpipol`) via spline.
      */
     void boundary_fmr();
     /**
      * @brief Set boundary values in the buffer zones by interpolating from the parent layer (`llay`).
      * @param btype Type of interpolation based on RK stage (1-5). Determines interpolation formula coefficients (aa, bb, cc).
      * @details Interpolates `llay->bv1`, `llay->bv0`, `llay->bv` based on `btype` and stores result in this layer's `bv` buffer zones.
      *          Calls `boundary_fmr` afterwards.
      */
     void set_boundary(int btype);
     /**
      * @brief Initialize this refinement layer's state by interpolating from the parent layer.
      * @details Sets up grid geometry (`set_flat`, `set_refz`), interpolates all state variables from `llay` to `bv`,
      *          calls `set_boundary(5)`, copies `bv` to `bv0`, and calculates initial primitive fluid variables (`dyntoprim`).
      */
     void set_fmr_initial();
     /**
      * @brief Perform two time steps on this layer and recursively call evolve on the child layer (`ulay`) if it exists.
      * @details Manages the 2:1 time stepping refinement. Calls `onestep` twice. If `mrf` is true, calls `ulay->evolve()`.
      *          Finally, calls `refine_llay` to update the parent layer with data from this layer.
      */
     void evolve();
     /**
      * @brief Update the parent layer (`llay`) data using data from this layer (restriction).
      * @details Copies values from this layer's grid points to the corresponding parent layer grid points.
      *          Effectively restricts data from fine to coarse grid in the overlapping region.
      */
     void refine_llay();
     /**
      * @brief Perform one full Runge-Kutta 4th order step on this layer.
      * @param btype Boundary interpolation type (passed to `set_boundary`).
      * @details Executes the 4 stages of RK4 (`BSSN(1)` to `BSSN(4)`), calling `set_boundary` after each stage.
      *          Updates time `t`. If `mrf` is true, recursively calls `ulay->evolve()` after the step.
      */
     void onestep(int btype);
     /**
      * @brief Load refinement layer data from a continuation file.
      * @param fcontinue Input file stream for the checkpoint file.
      * @details Calls `initial_continue`, `boundary_fmr`, and `dyntoprim`.
      */
     void fmr_initial_continue(ifstream& fcontinue);
 
 }; // End class Fmv1
 
 #endif // _COSMOS_S_H_
 