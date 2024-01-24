/**

for 2d code: -grid=multigrid
for 3d code: -grid=multigrid3D


qcc -grid=multigrid3D -lm -O3 convection.c -o cv.e

CC99='mpicc -std=c99' qcc -D_MPI=1  -D_NETCDF=1 -grid=multigrid3D -lm -lpnetcdf -O3 convection.c -o cv.e

HPC:
qcc -D_MPI=1  -D_NETCDF=1 -grid=multigrid3D -source convection.c
*/


double nu = 0.;
//double kappa = 1.e-5;
double N2b = 0.;  // background stratification
double B0 = 0.;   // surface buoyancy flux
double tau0 = 0.; // surface stress
double hm_i = 0.; // initial mixed layer thickness



int NMOUT = 5000;
int N0; // bug in restore
double tend = 1.;
double dtout = 1.;
double tout0 = 0;
char dpath[80]  = "./";
double tol_b = 0;
double tol_u = 0;
timer tel;
double tel_m = 1e10;
int npt = 0;
double nu0;
double b_noise = 1e-3;

// domain
int npz = 0;


//#include "grid/multigrid.h"
//#include "embed.h"
#include "../libs/extra.h"
#include "navier-stokes/centered.h"
#include "../libs/boussinesq.h"
#include "coriolis3d.h"
#include "auxiliar_input.h"
#if _NETCDF
#include "pnetcdf_bas.h"
#endif
#include "fractions.h"

//bug in netcdf
//int nl = 0;

int main(int argc,char* argv[]) {
  tel = timer_start();

  params = array_new();
  add_param ("N", &N, "int");
  add_param ("npz", &npz, "int");
  add_param ("L0", &L0, "double");
  add_param ("nu", &nu, "double");
  add_param ("b_noise", &b_noise, "double");
  add_param ("N2b", &N2b, "double");
  add_param ("B0", &B0, "double");
  add_param ("tau0", &tau0, "double");
  add_param ("hm_i", &hm_i, "double");
  add_param ("f0_x", &f0_x, "double");
  add_param ("f0_y", &f0_y, "double");
  add_param ("f0_z", &f0_z, "double");
  add_param ("alpha_H", &alpha_H, "double");
  add_param ("DT", &DT, "double");
  add_param ("CFL", &CFL, "double");
  add_param ("tend", &tend, "double");
  add_param ("dtout", &dtout, "double");
  add_param ("tout0", &tout0, "double");
  add_param ("NMOUT", &NMOUT, "int");

  // origin
  X0 = 0;
#if dimension > 2
  Y0 = 0;
  Z0 = -L0; // will be adjsuted in init event in case non cubic domain
#else
  Y0 = -L0;
#endif

  N0 = N; // bug in restore

  // Search for the configuration file with a given path or read params.in
  if (argc == 2)
    strcpy(file_param,argv[1]); // default: params.in

  read_params(file_param);
  create_outdir();
  backup_config(file_param);
  

  //TOLERANCE = 1e-4;
  periodic(right);
#if dimension > 2
  periodic(top);
  if (npz > 0){
#if _MPI
    int rank, comm_size;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &comm_size);
#endif
    dimensions(nz=npz);
  }
#endif

  u.x.gradient = minmod;
  u.y.gradient = minmod;
  u.z.gradient = minmod;
  run();
}

/**
## Initial conditions
Initial conditions correspond to a tanh profile for both buoyancy and velocity
*/



event init (t=0) {
#if _MPI
  Z0 = -L0*mpi_dims[2]/mpi_dims[0];
#endif


  if (!restore (file = "restart.bin")) {
    FILE * fp;
    if ((fp = fopen("b0.bas", "r"))) {
      input_matrixl ({b}, fp, ox=X0, oy=Y0);
      fclose(fp);
    } else {
      foreach() {
#if dimension > 2
        b[] = z < -hm_i ? N2b*(z + hm_i) : 0.;
#else
        b[] = y < -hm_i ? N2b*(y + hm_i) : 0.;
#endif
        b[] += b_noise*N2b*L0*noise();
      }
    }
    
    if ((fp = fopen("u0.bas", "r"))) {
      input_matrixl ((scalar *) {u}, fp, ox=X0, oy=Y0);
      fclose(fp);
    } else {
      foreach() {
        u.x[] = 0.;
        u.y[] = 0.;
      }
    }
  }
  N = N0; // restore function changes N


//#define WAVES (-1*0.5*(1-sin(8*2*pi*x/L0)) - y - 1.4987654321*L0/N)
#define WAVES (-0.5*(-0.5+sin(8*2*pi*x/L0-sqrt(tau0)*t/L0)) - y)
//#define WAVES (-1*0.5*(1-sin(8*2*pi*x/L0)) - y - 1.4987654321*L0/2)
/*   solid (cs, fs, WAVES); */


  /* double epsilon = 15*H_w*pow(sqrt(tau0),3.)/sq(+z0)/k_vk; */
  /* double l_charnock = k_vk*beta_charnock*tau0/g_grav; */
  /* nu0 = nu + Cd*pow(epsilon,1./3.)*pow(l_charnock, 4./3.); */
  /* fprintf(stdout,"nu0 = %e\n",nu0); */
 

#if dimension > 2
  /* b[front] = neumann(q0/nu); */
//  u.t[front] = neumann(tau0/nu);
  b[back] = neumann(-N2b);
#else
//  b[top] = neumann(q0/nu);
  b[bottom] = neumann(-N2b);
  /* b[embed] = neumann(q0/kappa); */
  
//  u.r[top] = WAVES < 0 ? dirichlet(0) : neumann(tau0/nu);
//  u.t[top] = neumann(tau0/nu);
//  u.t[top] = WAVES < 0 ? dirichlet(0) : neumann(tau0/nu);
  /* u.t[embed] = tau0/nu*Delta + u.x[0,-1]; */
  /* u.n[embed] = dirichlet(0.); */
#endif

  boundary (all);

#if _NETCDF
  char file_nc[80];
  sprintf (file_nc, "%svars.nc", dpath); 
  create_nc({b,u}, file_nc);
#endif
}


/**
   Surface forcing (buoyancy flux)
*/
// trying this formulation instead of specification of neumann BC so that it is
// independent of nu
event surface_fluxes (i++) {
#if dimension > 2
  foreach_boundary(front) {
#else
  foreach_boundary(top) {
#endif
    b[]   += dt*B0/Delta;
//    b[]   += dt*q0/Delta;
//    u.x[] += dt*tau0/Delta;
  }
}

/**
   Surface forcing (momentum flux)
*/

event acceleration (i++) {

#if dimension > 2
  foreach_boundary(front) {
#else
  foreach_boundary(top) {
#endif
    u.x[] += dt*tau0/Delta;
  }

}

/* scalar waves[]; */

/* event moving_cylinder (i++) { */
/*   coord vc = {0.,0.}; // the velocity of the cylinder */
/*   fraction (waves, -WAVES); */

/*   foreach() */
/*     foreach_dimension() */
/*       u.x[] = waves[]*vc.x + (1. - waves[])*u.x[]; */
/*   boundary ((scalar *){u}); */
/* } */

event output (t = tout0; t <= tend+1e-10;  t += dtout) {

  double ket = 0., vd = 0.;
  foreach(reduction(+:ket) reduction(+:vd)) {
      ket += 0.5*(sq(u.x[]) + sq(u.y[]))*sq(Delta);
      //      viscous dissipation
      vd += sq(Delta)*(sq(u.x[1] - u.x[-1]) +
                          sq(u.x[0,1] - u.x[0,-1]) +
                          sq(u.y[1] - u.y[-1]) +
                          sq(u.y[0,1] - u.y[0,-1])
                          )/sq(2.*Delta);
  }
  vd *= nu;
  

  fprintf (stdout,"i = %d, dt = %g, t = %g, ke = %.2e, vd = %.2e\n",
           i, dt, t, ket, vd);

  if (N < NMOUT){
#if _NETCDF
    write_nc();
#else
    char name[100];
    sprintf (name,"%silast.dat", dpath);
    FILE * fp = fopen (name, "w");
    fprintf(fp,"%d",i);
    fclose(fp);
  
    sprintf (name,"%sb%09d.bas", dpath, i);
    fp = fopen (name, "w");
    output_matrixl ({b}, fp, linear = true);
    fclose(fp);
  
    sprintf (name,"%su%09d.bas", dpath, i);
    fp = fopen (name, "w");
    output_matrixl ((scalar *) {u}, fp, linear = true);
    fclose(fp);
  
#if TREE
    if (tol_b || tol_u) {
      scalar l[];
      foreach()
        l[] = level;
      boundary({l});
      sprintf (name,"%slevel%09d.bas", dpath, i);
      fp = fopen (name, "w");
      output_matrixl ({l}, fp, linear = true);
      fclose(fp);
    }
#endif
#endif
  } else {
    char name[100];
    sprintf (name, "%soutput-%09d", dpath, i);
    fprintf(stdout,name);
    dump (name);
  }

}

event snapshot(i++){
  double elapsed = timer_elapsed(tel);

  if ( elapsed > tel_m) {
    char name[100];
    sprintf (name, "%ssnapshot-%09d", dpath, i);
    dump (name);
    return 1;
  }
}
