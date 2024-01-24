/**
   Netcdf interface for basilisk

   TODO finish level implementation
   TODO for non level fields, no need to have a z dimension
*/

#include <stdio.h>
#include <string.h>
#include <pnetcdf.h>

#define NDIMS 4
#define Z_NAME "z"
#define Y_NAME "y"
#define X_NAME "x"
#define REC_NAME "time"
#define LVL_NAME "level"

#if LAYERS == 0
int nl = 1;
#endif


/* /\* For the units attributes. *\/ */
/* #define UNITS "units" */
/* #define PRES_UNITS "hPa" */
/* #define TEMP_UNITS "celsius" */
/* #define MAX_ATT_LEN 80 */

/* Handle errors by printing an error message and exiting with a
 * non-zero status. */
int nc_err;

static void handle_error(int status, int lineno)
{
    fprintf(stderr, "Error at line %d: %s\n", lineno, ncmpi_strerror(status));
    MPI_Abort(MPI_COMM_WORLD, 1);
}
/* IDs for the netCDF file, dimensions, and variables. */
int ncid;

// User global variables
scalar * scalar_list_nc;
char file_nc[80];


// temporary
int nc_varid[1000];
char * nc_varname[1000];
int nvarout = 0;
int nc_rec = -1;
int it_nc = -1;
int rec_varid;

struct OutputNetcdf {
  scalar * scalar_list;
  char * file_nc;
  int nl;
  int it;
};


void create_nc(scalar * list_out, char* file_out){
   /* LOCAL IDs for the netCDF file, dimensions, and variables. */
   int x_dimid, y_dimid, z_dimid, rec_dimid;
   int z_varid, y_varid, x_varid;
   int dimids[NDIMS];
   
   // make it global variable
   sprintf (file_nc,"%s", file_out);
   scalar_list_nc = list_copy(list_out);


   /* Create the file. */
   if ((nc_err = ncmpi_create(MPI_COMM_WORLD, file_nc,
                              NC_CLOBBER, MPI_INFO_NULL,&ncid)))
     handle_error(nc_err, __LINE__);

   /* Define the dimensions. The record dimension is defined to have
    * unlimited length - it can grow as needed. In this example it is
    * the time dimension.*/
#if _MPI
   int npx = mpi_dims[0];
   int npy = mpi_dims[1];
#else
   int npx = 1;
   int npy = 1;
#endif
   int Nloc = (1 << depth());
   int Nx = Nloc*npx;
   int Ny = Nloc*npy;

   if ((nc_err = ncmpi_def_dim(ncid, REC_NAME, NC_UNLIMITED, &rec_dimid)))
     handle_error(nc_err, __LINE__);

#if dimension > 2
#if _MPI
   int npz = mpi_dims[2];
#else
   int npz = 1;
#endif
   int Nz = Nloc*npz;
   if ((nc_err = ncmpi_def_dim(ncid, Z_NAME, Nz, &z_dimid)))
     handle_error(nc_err, __LINE__);
#else
   if ((nc_err = ncmpi_def_dim(ncid, LVL_NAME, nl, &z_dimid)))
     handle_error(nc_err, __LINE__);
#endif
   if ((nc_err = ncmpi_def_dim(ncid, Y_NAME, Ny, &y_dimid)))
     handle_error(nc_err, __LINE__);
   if ((nc_err = ncmpi_def_dim(ncid, X_NAME, Nx, &x_dimid)))
     handle_error(nc_err, __LINE__);

   /* Define the coordinate variables. We will only define coordinate
      variables for lat and lon.  Ordinarily we would need to provide
      an array of dimension IDs for each variable's dimensions, but
      since coordinate variables only have one dimension, we can
      simply provide the address of that dimension ID (&y_dimid) and
      similarly for (&x_dimid). */

   if ((nc_err = ncmpi_def_var(ncid, REC_NAME, NC_FLOAT, 1, &rec_dimid,
			    &rec_varid)))
     handle_error(nc_err, __LINE__);
#if dimension > 2
   if ((nc_err = ncmpi_def_var(ncid, Z_NAME, NC_FLOAT, 1, &z_dimid,
			    &z_varid)))
     handle_error(nc_err, __LINE__);
#endif
   if ((nc_err = ncmpi_def_var(ncid, Y_NAME, NC_FLOAT, 1, &y_dimid,
			    &y_varid)))
     handle_error(nc_err, __LINE__);
   if ((nc_err = ncmpi_def_var(ncid, X_NAME, NC_FLOAT, 1, &x_dimid,
			    &x_varid)))
     handle_error(nc_err, __LINE__);

   /* The dimids array is used to pass the dimids of the dimensions of
      the netCDF variables. Both of the netCDF variables we are
      creating share the same four dimensions. In C, the
      unlimited dimension must come first on the list of dimids. */
   dimids[0] = rec_dimid;
   dimids[1] = z_dimid;
   dimids[2] = y_dimid;
   dimids[3] = x_dimid;

   /* Define the netCDF variables */
//   char * str1 = "empty_var";
   for (scalar s in scalar_list_nc){
//     if (strcmp(str1,s.name) != 0) {
       if ((nc_err = ncmpi_def_var(ncid, s.name, NC_FLOAT, NDIMS,
                                   dimids, &nc_varid[nvarout])))
         handle_error(nc_err, __LINE__);
//       nc_varname[nvarout] = strdup(s.name);
       nvarout += 1;
//       str1 = strdup(s.name);
//     }
   }
   
   /* /\* Assign units attributes to the netCDF variables. *\/ */
   /* if ((nc_err = nc_put_att_text(ncid, pres_varid, UNITS,  */
   /*      			 strlen(PRES_UNITS), PRES_UNITS))) */
   /*    ERR(nc_err); */
   /* if ((nc_err = nc_put_att_text(ncid, temp_varid, UNITS,  */
   /*      			 strlen(TEMP_UNITS), TEMP_UNITS))) */
   /*    ERR(nc_err); */

   /* End define mode. */
   if ((nc_err = ncmpi_enddef(ncid)))
     handle_error(nc_err, __LINE__);

   /*  write coordinates*/
   float yc[Ny], xc[Nx];
   double Delta = L0*1.0/Nx;
   for (int i = 0; i < Nx; i++)
     xc[i] = X0 + (i + 0.5)*Delta;
   
   for (int i = 0; i < Ny; i++)
      yc[i] = Y0 + (i + 0.5)*Delta;

#if dimension > 2
   float zc[Nz];
   for (int i = 0; i < Nz; i++)
     zc[i] = Z0 + (i + 0.5)*Delta;
   if ((nc_err = ncmpi_put_var_float_all(ncid, z_varid, &zc[0])))
     handle_error(nc_err, __LINE__);
#endif

   if ((nc_err = ncmpi_put_var_float_all(ncid, y_varid, &yc[0])))
     handle_error(nc_err, __LINE__);
   if ((nc_err = ncmpi_put_var_float_all(ncid, x_varid, &xc[0])))
     handle_error(nc_err, __LINE__);

   /* Close the file. */
   if ((nc_err = ncmpi_close(ncid)))
     handle_error(nc_err, __LINE__);
   fprintf(stdout,"*** SUCCESS creating example file %s!\n", file_nc);
}


void write_nc(struct OutputNetcdf p) {

  int Nloc = (1 << depth());

  /* open file. */
  if ((nc_err = ncmpi_open(MPI_COMM_WORLD, 
                           file_nc, NC_WRITE, MPI_INFO_NULL, &ncid)))
    handle_error(nc_err, __LINE__);
  
  /* if (it_nc != p.it){ */
  /*   it_nc = p.it; */
  /*   nc_rec += 1; */
  /* } */

  // write time
  nc_rec += 1;
  float loctime[1];
  loctime[0] = t;
  

  MPI_Offset startt[1], countt[1];
  startt[0] = nc_rec; //time
  countt[0] = 1;

  if ((nc_err = ncmpi_put_vara_float_all(ncid, rec_varid, startt, countt,
                                         &loctime[0])))
    handle_error(nc_err, __LINE__);
  
#if dimension > 2
  float * field = (float *)malloc(Nloc*Nloc*Nloc*sizeof(float));
#else
//  float ** field = matrix_new (Nloc, Nloc, sizeof(float));
  float * field = (float *)malloc(Nloc*Nloc*sizeof(float));
#endif
  
  /* The start and count arrays will tell the netCDF library where to
     write our data. */
  MPI_Offset start[NDIMS], count[NDIMS];
  
  
  /* These settings tell netcdf to write one timestep of data. (The
     setting of start[0] inside the loop below tells netCDF which
     timestep to write.) */
  start[0] = nc_rec; //time
#if dimension > 2
#if _MPI
  start[1] = mpi_coords[2]*Nloc;     //z
#else
  start[1] = 0;
#endif  // MPI
#else
  start[1] = 0;     //level
#endif
#if _MPI
  start[2] = mpi_coords[1]*Nloc;      //y
  start[3] = mpi_coords[0]*Nloc;      //x
#else
  start[2] = 0;      //y
  start[3] = 0;      //x
#endif
  
  count[0] = 1;
#if dimension > 2
  count[1] = Nloc;     //z
#else
  count[1] = nl;
#endif
  count[2] = Nloc;
  count[3] = Nloc;
  
  
   int nv = -1;
 /* char * str1; */
  for (scalar s in scalar_list_nc){
    nv += 1;
    foreach(noauto){
      int i = point.i - GHOSTS;
      int j = point.j - GHOSTS;
#if dimension > 2
      int k = point.k - GHOSTS;
      field[Nloc*Nloc*k + Nloc*j + i] = s[];
//      fprintf(stdout,"val= %f !\n", field[k][j][i]);
#else
      field[Nloc*j + i] = s[];
//      field[j][i] = s[];
#endif
    }
  
/*     int nv; */
/*     for (nv = 0; nv < nvarout; nv ++) */
/*       if (strcmp(s.name, nc_varname[nv]) == 0) { */
/* //        start[1] += 1; */
/*         break; */
/*       } */
    
    if ((nc_err = ncmpi_put_vara_float_all(ncid, nc_varid[nv], start, count,
                                           &field[0])))
      //                                           &field[0][0])))
      handle_error(nc_err, __LINE__);
    
    /* if (start[1] == p.nl - 1) */
    /*   start[1] = -1; */
  }
  free(field);

   /* Close the file. */
   if ((nc_err = ncmpi_close(ncid)))
     handle_error(nc_err, __LINE__);
}


/**
   Read NC
 */

void read_nc(scalar * list_in, char* file_in){

  int i, ret;
  int ncfile, ndims, nvars, ngatts, unlimited;
  int var_ndims, var_natts;
  nc_type type;
  char varname[NC_MAX_NAME+1];
  int *dimids=NULL;


  int Nloc = (1 << depth());
#if dimension > 2
  float * field = (float *)malloc(Nloc*Nloc*Nloc*sizeof(float));
#else
//  float ** field = matrix_new (Nloc, Nloc, sizeof(float));
  float * field = (float *)malloc(Nloc*Nloc*nl*sizeof(float));
#endif
 

//  float ** field = matrix_new (Nloc, Nloc, sizeof(float));


  ret = ncmpi_open(MPI_COMM_WORLD, file_in, NC_NOWRITE, MPI_INFO_NULL,
                   &ncfile);

  if (ret != NC_NOERR) handle_error(ret, __LINE__);
    /* reader knows nothing about dataset, but we can interrogate with query
     * routines: ncmpi_inq tells us how many of each kind of "thing"
     * (dimension, variable, attribute) we will find in the file  */

    /* no communication needed after ncmpi_open: all processors have a cached
     * view of the metadata once ncmpi_open returns */

    ret = ncmpi_inq(ncfile, &ndims, &nvars, &ngatts, &unlimited);
    if (ret != NC_NOERR) handle_error(ret, __LINE__);


  MPI_Offset start[NDIMS], count[NDIMS];
  // default (no MPI)
  start[0] = 0; //time
  start[1] = 0; //level
  start[2] = 0;
  start[3] = 0;
#if _MPI
#if dimension > 2
  start[1] = mpi_coords[2]*Nloc; 
#endif
  start[2] = mpi_coords[1]*Nloc;      //y
  start[3] = mpi_coords[0]*Nloc;      //x
#endif  

  count[0] = 1;
#if dimension > 2
  count[1] = Nloc;
#else
  count[1] = nl; // TODO finish level implementation
#endif
  count[2] = Nloc;
  count[3] = Nloc;




  for (scalar s in list_in){
    for(i=0; i<nvars; i++) {


      ret = ncmpi_inq_var(ncfile, i, varname, &type, &var_ndims, dimids,
                          &var_natts);
      if (ret != NC_NOERR) handle_error(ret, __LINE__);


      if (strcmp(varname,s.name) == 0) {
        fprintf(stdout,"Reading variable  %s!\n", s.name);

          if ((nc_err = ncmpi_get_vara_float_all(ncfile, i, start, count,
                                                 &field[0])))
            handle_error(nc_err, __LINE__);

          foreach(noauto){
            int i = point.i - GHOSTS;
            int j = point.j - GHOSTS;
#if dimension > 2
            int k = point.k - GHOSTS;
             s[] = field[Nloc*Nloc*k + Nloc*j + i];
#else
            s[] = field[Nloc*j + i];
#endif
//            s[] = field[j][i];

          }


        }

    }
  }


//  matrix_free (field);
  free(field);

  ret = ncmpi_close(ncfile);
  if (ret != NC_NOERR) handle_error(ret, __LINE__);

}

event cleanup (t = end)
{
  free(scalar_list_nc);
}
