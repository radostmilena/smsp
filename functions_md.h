#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
using namespace std;
using namespace Eigen;

//useful structures
struct result{
        long double r, en, sumv, sumv2;

};

struct prop_conv{
	long double en, temp, ekin, nstep_ps;
};

//routine to initialize positions and generate random velocities
int init_posv(long double **coord, long double **v, long double d, long double &vsum, long double &v2sum, int n, int npart){

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

      //assign random velocities
      for(a=0;a<npart;a++){
         for(b=0;b<3;b++){
          v[a][b]=(rand()/(RAND_MAX+1.)-0.5);
          vsum=vsum+v[a][b];
          long double curv=v[a][b];
          v2sum=v2sum+pow(curv, 2);
          }
      }
      
      return 0;

}

//routine for velocity rescaling
int vel_resc(long double **v, long double vsum, long double v2sum, int npart, long double temp){

      int a, b;
      long double fs=sqrt((3*temp)/v2sum);
      for(a=0;a<npart;a++){
          for(b=0;b<3;b++){
            v[a][b]=(v[a][b]-vsum)*fs;
          }
      }

      return 0;

}

//routine to calculate LJ potential
result ljpot(int npart, long double box, long double **coord, long double rcut, long double **f, long double *g){
       result rs;
       rs.en=0;
       int a, b, c;
       long double rvec[3];
       long double ecut=4*(pow((1/rcut), 12)-pow((1/rcut), 6));
       long double r2i, r6i, ff;
       long double dr=box/100; //spacing for RDF
       long double hist[100] = { 0 }; //histogram for RDF

       //set forces to zero
       for(int a = 0; a < npart; a++){
	  for(b=0; b<3;b++){
              f[a][b] = 0;
	  }
       }

       //loop through all pairs
       for(a=0;a<npart-1;a++){
          for(b=a+1;b<npart;b++){
              rs.r=0;
              for(c=0; c<3; c++){
	          //calculate distance vector
	          rvec[c] = coord[a][c]-coord[b][c];
	          //apply MIC PBC
	          rvec[c] = rvec[c]-box*round(rvec[c]/box);
                  //square of distance vector
                  rs.r = rs.r+pow(rvec[c], 2);
	      }

	      //calculate distance (scalar)
	      rs.r = sqrt(rs.r);

	      //add to RDF
	      int index = round(rs.r/dr);
	      hist[index]= hist[index]+1;

	      //check cutoff
	      if(rs.r<rcut){

	      //calculate force
              r2i=1/(pow(rs.r, 2));
              r6i=pow(r2i, 3);
              ff=48*r2i*r6i*(r6i-0.5);

	      for(c=0; c<3; c++){
	         //update force on particles a and b
                 f[a][c]=f[a][c]+rvec[c]*ff;
                 f[b][c]=f[b][c]-rvec[c]*ff;
	      }

	      //calculate energy
              rs.en=rs.en+4*r6i*(r6i-1)-ecut;
              }
          }
       }

       //calculate RDF
       for(a=0; a<100; a++){
	  float N = (float) npart;
          g[a] = (pow(box, 3)/(N*4*M_PI*dr))*(hist[a]*(2/N))*(1/pow((dr*a), 2));
       }

       return rs;
}

//routine for time integration
result t_int(int npart, long double **f, long double **coord, long double **v, long double dt, long double box, long double rcut, long double *g){
         result rs;
	 rs.sumv=0;
	 rs.sumv2=0;
         int a, b;
         long double *prev_f[npart];
         for(int a = 0; a < npart; a++){
             prev_f[a] = new long double[3];
	 }

          //calculate new positions
          for(a=0;a<npart;a++){
             for(b=0; b<3; b++){
                 prev_f[a][b]=f[a][b];
                 coord[a][b]=coord[a][b]+v[a][b]*dt+((prev_f[a][b])/2)*pow(dt, 2);
	     }
	  }

         //calculate new forces
         result res2=ljpot(npart, box, coord, rcut, f, g);
	 rs.en=res2.en;

         //calculate new velocities
         for(a=0;a<npart;a++){
             for(b=0; b<3; b++){
                 v[a][b]=v[a][b]+((f[a][b]+prev_f[a][b])/2)*dt;
                 rs.sumv=rs.sumv+v[a][b];
                 rs.sumv2=rs.sumv2+pow(v[a][b], 2);
	     }
          }

         return rs;
}

//routine for calculating and converting properties
prop_conv calcprop_conv(int npart, long double sumv2, long double **v, long double nstep, long double **coord, long double en, long double box, long double temp, long double ekin, long double **vcur, long double **coordwrap){
      prop_conv pc;

      //Ar parameters and other
      long double sigma=3.4e-10, m=39.95*1.6747e-27, epsilon=1.6568e-21;
      long double k_B=1.3806504e-23;

      int a, b;

      //wrap and convert units on coordinates and velocities
      for (a=0;a<npart;a++){
	   for(b=0;b<3;b++){
               coordwrap[a][b]=coord[a][b]*sigma-box*sigma*floor(coord[a][b]/(box));
               vcur[a][b]=v[a][b]*sqrt(epsilon/m);
	   }
      }

      //convert units on other properties
      pc.temp=temp*(epsilon/k_B);
      pc.ekin=ekin*epsilon;
      pc.en=en*epsilon;

      //convert units on timestep
      pc.nstep_ps=nstep*(sqrt(epsilon/m)*(1/sigma))*1e-12;

      return pc;

}
