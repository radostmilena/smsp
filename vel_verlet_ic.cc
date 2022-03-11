#include <iostream>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <time.h>
#include <fstream>
#include "functions_md.h"
using namespace std;
using namespace Eigen;

int main(int argc, char* argv[]) {

      ofstream file1, file2, file3;

      //Ar parameters
      long double sigma=3.4e-10, m=39.95*1.6747e-27, epsilon=1.6568e-21;
      long double rho=1374, temp=94.4, k_B=1.3806504e-23, rcut=8.0e-10/sigma;
      temp=temp*(k_B/epsilon);

      //initialize some variables
      int a, b;
      long double nstep=0, dt, ntot, box, vsum=0, v2sum=0, ekin, itemp;
      long double g[100] = { 0 };

      //read in timestep
      cin >> dt;
      ntot=20/dt;
      dt=dt*(sqrt(m/epsilon)*sigma)*1e12;

      //initialize system
      int n = 10;
      int npart = pow(n, 3);
      long double d=cbrt(m/rho)/sigma;
      box = n*d;
      long double *coord[npart], *v[npart], *f[npart], *vcur[npart], *coordwrap[npart];
      for(int a = 0; a < npart; a++){
          coord[a] = new long double[3];
          v[a] = new long double[3];
          f[a] = new long double[3];
          vcur[a] = new long double[3];
          coordwrap[a] = new long double[3];
      }

      int init = init_posv(coord, v, d, vsum, v2sum, n, npart);

      //velocity rescaling
      vsum=vsum/npart;
      v2sum=v2sum/npart;
      int vresc = vel_resc(v, vsum, v2sum, npart, temp);
      
      //MD loop
      time_t start,end;
      time (&start);

      file1.open("energies.txt");
      file2.open("traj.txt");
      file3.open("rdf.txt");

      while(nstep<=ntot){

          //calculate force and energy
          result res=ljpot(npart, box, coord, rcut, f, g);

	  //propagate coord and v
          result res2=t_int(npart, f, coord, v, dt, box, rcut, g);

	  //calculate temperature and total energy in red. units
          itemp=(res2.sumv2)/(3*npart);
          ekin=0.5*res2.sumv2;

	  //convert units, wrap
	  prop_conv pco=calcprop_conv(npart, res2.sumv2, v, nstep, coord, res2.en, box, itemp, ekin, vcur, coordwrap);

          //print energies
	  file1 << pco.nstep_ps*dt << ' ' << pco.temp << ' ' << pco.ekin << ' ' << pco.en << ' ' << pco.ekin+pco.en << endl;

          //print coordinates, velocities and wrap if specified
	  if(strcmp(argv[1], "-wrap") == 0){
		  for(a=0; a<npart; a++){
	              file2 << coordwrap[a][0] << ' ' << coordwrap[a][1] << ' ' << coordwrap[a][2] << ' ' << vcur[a][0] << ' ' << vcur[a][1] << ' ' << vcur[a][2] << endl;
		  }
	  }
	  else{
		  for(a=0; a<npart; a++){
	              file2 << coord[a][0] << ' ' << coord[a][1] << ' ' << coord[a][2] << ' ' << vcur[a][0] << ' ' << vcur[a][1] << ' ' << vcur[a][2] << endl;
		  }
          }

	  //print RDF
	  for(a=0; a<100; a++){
	     file3 << a*(box/100)*sigma << ' ' << g[a] << endl;
	  }

	  //increase timestep
	  nstep+=1;

	  //velocity rescaling
	  res2.sumv2 = (res2.sumv2/npart);
	  vsum = 0;
	  if(int(nstep)%50==0 && nstep < 1000){
             int v_resc = vel_resc(v, vsum, res2.sumv2, npart, temp);
	  }
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
