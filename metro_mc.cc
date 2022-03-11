#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <time.h>
#include <fstream>
#include "functions_mc.h"
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {

      ofstream file1, file2, file3;

      //Ar parameters
      long double sigma=3.4e-10, m=39.95*1.6747e-27, epsilon=1.6568e-21;
      long double rho=1374, temp=94.4, k_B=1.3806504e-23, rcut=8.0e-10/sigma;
      temp=temp*(k_B/epsilon);
      long double beta = (1/(temp));

      //initialize some variables
      int a, b;
      long double nstep=0, ntot=5000000, box, ekin, dx=0.20;
      long double g[100] = { 0 };

      //initialize system
      int n = 10;
      int npart = pow(n, 3);
      long double d=cbrt(m/rho)/sigma;
      box = n*d;
      long double *coord[npart], *coordwrap[npart];
      for(int a = 0; a < npart; a++){
          coord[a] = new long double[3];
          coordwrap[a] = new long double[3];
      }

      int init = init_posv(coord, d, n, npart);

      ekin=1.5*npart*k_B*temp*(epsilon/k_B); //in joule

      //MC loop
      time_t start,end;
      time (&start);

      file1.open("energies.txt");
      file2.open("traj.txt");
      file3.open("rdf.txt");

      long double epot = ljpot(npart, box, coord, rcut);

      while(nstep<=ntot){

	  //MC move
	  long double delta_en = mcmove(npart, box, coord, beta, rcut, dx); 

	  epot = epot + delta_en;

	  //convert units, wrap
	  long double epot_conv = prop_conv(npart, coord, box, coordwrap, epot);

          //print energies
	  file1 << nstep << ' ' << ekin << ' ' << epot_conv << ' ' << ekin+epot_conv << endl;

          //print coordinates, wrap if specified
	  if (delta_en != 0){
	  if(strcmp(argv[1], "-wrap") == 0){
		  file2 << npart << endl;
		  file2 << 's' << 't' << 'e' << 'p' << '=' << nstep << endl;
		  for(a=0; a<npart; a++){
	              file2 << 'A' << ' ' << coordwrap[a][0]*1e10 << ' ' << coordwrap[a][1]*1e10 << ' ' << coordwrap[a][2]*1e10 << endl;
		  }
	  }
	  else{
		  file2 << npart << endl;
		  file2 << 's' << 't' << 'e' << 'p' << '=' << nstep << endl;
		  for(a=0; a<npart; a++){
	              file2 << 'A' << ' ' << coord[a][0] << ' ' << coord[a][1] << ' ' << coord[a][2] << endl;
		  }
          }
	  }

	  //calculate and print RDF
	  if(int(nstep)%1000==0){
		  int rdf = calc_rdf(npart, coord, box, g);
	      for(a=0; a<100; a++){
	          file3 << a*(box/100)*sigma << ' ' << g[a] << endl;
	      }
	  }

	  //increase step
	  nstep+=1;

      }  

      //monitor performance
      time (&end);
      double dif = difftime (end,start);
      printf("Elapsed time is %.2lf seconds. \n", dif );

      file1.close();
      file2.close();
      file3.close();

      return 0;
}
