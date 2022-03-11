#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

//routine to initialize positions on lattice
int init_posv(long double **coord, long double d, int n, int npart){

      //set initial positions
      int count=0, i, j, k, a, b;
      for(i=0;i<n;i++){
         for(j=0;j<n;j++){
            for(k=0;k<n;k++){
                coord[count][0]=d*i;
                coord[count][1]=d*j;
                coord[count][2]=d*k;
                count+=1;
            }
         }
      }

      return 0;

}

//routine to calculate LJ potential for whole system
long double ljpot(int npart, long double box, long double **coord, long double rcut){
       long double en=0, rvec[3], ecut=4*(pow((1/rcut), 12)-pow((1/rcut), 6)), r2i, r6i;
       int a, b, c;

       //loop through pairs
       for(a=0;a<npart-1;a++){
          for(b=a+1;b<npart;b++){
              long double r=0;
              for(c=0; c<3; c++){
	          //calculate distance vector
	          rvec[c] = coord[a][c]-coord[b][c];
	          //apply MIC PBC
	          rvec[c] = rvec[c]-box*round(rvec[c]/box);
                  //square of distance vector
                  r = r+pow(rvec[c], 2);
	      }

	      //calculate distance (scalar)
	      r = sqrt(r);

	      //check cutoff
	      if(r<rcut){

	      //calculate energy
              r2i=1/(pow(r, 2));
              r6i=pow(r2i, 3);
              en=en+4*r6i*(r6i-1)-ecut;

	      }
          }
       }

       return en;
}

//routine to calculate LJ potential energy difference
long double delta_en(int npart, long double box, long double **coord, long double rcut, int sel_part){
       long double en=0, rvec[3], ecut=4*(pow((1/rcut), 12)-pow((1/rcut), 6)), r2i, r6i;
       int a, b, c;

       //loop through pairs
       for(a=sel_part;a<sel_part+1;a++){
          for(b=0;b<npart;b++){
              long double r=0;
	      if(b!=a){
              for(c=0; c<3; c++){
	          //calculate distance vector
	          rvec[c] = coord[a][c]-coord[b][c];
	          //apply MIC PBC
	          rvec[c] = rvec[c]-box*round(rvec[c]/box);
                  //square of distance vector
                  r = r+pow(rvec[c], 2);
	      }

	      //calculate distance (scalar)
	      r = sqrt(r);

	      //check cutoff
	      if(r<rcut){

	      //calculate energy
              r2i=1/(pow(r, 2));
              r6i=pow(r2i, 3);
              en=en+4*r6i*(r6i-1)-ecut;

              }
	      }
          }
       }

       return en;
}

//routine for displacing particle (Metropolis MC move)
long double mcmove(long double npart, long double box, long double **coord, long double beta, long double rcut, long double dx){
      int sel_part, c;
      long double old_coord[3], den;

      //randomly select particle
      sel_part = int(rand()/(RAND_MAX+1.)*npart);

      //calculate energy
      long double old_en = delta_en(npart, box, coord, rcut, sel_part);
      for(c=0;c<3;c++){
          old_coord[c] = coord[sel_part][c];
      }

      //displace particle by random number
      for(c=0;c<3;c++){
          long double dxyz = (rand()/(RAND_MAX+1.)-0.5)*dx;
          coord[sel_part][c] = coord[sel_part][c]+(dxyz); 
      }

      //calculate new energy
      long double new_en = delta_en(npart, box, coord, rcut, sel_part);

      //check acceptance
      long double factor= exp(-beta*(new_en-old_en));
      long double rd_nb = (rand()/(RAND_MAX+1.));
      if((rd_nb)>factor){
         for(c=0;c<3;c++){
	      coord[sel_part][c] = old_coord[c]; //not accepted
	      den = 0;
         }
      }
      else{
	      den = (new_en-old_en);
      }

      return den;
}

//routine for converting properties
long double prop_conv(int npart, long double **coord, long double box, long double **coordwrap, long double epot){
      //Ar parameters and other
      long double sigma=3.4e-10, epsilon=1.6568e-21;

      int a, b;

      //wrap and convert units on coordinates
      for (a=0;a<npart;a++){
	   for(b=0;b<3;b++){
               coordwrap[a][b]=coord[a][b]*sigma-box*sigma*floor(coord[a][b]/(box));
	   }
      }

      epot = epot*epsilon;

      return epot;
}

//routine for calculating RDF
int calc_rdf(int npart, long double **coord, long double box, long double *g){

      int a, b, c;
      long double rvec[3];
      long double dr=box/100; //spacing for RDF
      long double hist[100] = { 0 }; //histogram for RDF

      //loop through pairs
      for(a=0;a<npart-1;a++){
         for(b=a+1;b<npart;b++){
             long double r=0;
             for(c=0; c<3; c++){
                 //calculate distance vector
                 rvec[c] = coord[a][c]-coord[b][c];
                 //apply MIC PBC
                 rvec[c] = rvec[c]-box*round(rvec[c]/box);
                 //square of distance vector
                 r = r+pow(rvec[c], 2);
             }

             //calculate distance (scalar)
             r = sqrt(r);

             //add to RDF
             int index = round(r/dr);
             hist[index]= hist[index]+1;

         }
      }

      //calculate RDF
      for(a=0; a<100; a++){
         float N = (float) npart;
         g[a] = (pow(box, 3)/(N*4*M_PI*dr))*(hist[a]*(2/N))*(1/pow((dr*a), 2));
      }

      return 0;

}
