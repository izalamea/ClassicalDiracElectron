#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <signal.h> 
#include <string.h>

const int Nparticles=1;
const int Ndist=1000;

int main(){
  double r, rmin=3.74, rmax=8.54;
  double x,y,z;
  int it, dist[Ndist];
  char name[100];
  FILE *IfPtr;
 
  sprintf(name,"script_%d_Trajectories_.gp",Nparticles);
  if ( ( IfPtr = fopen( name, "w" ) ) == NULL )
    printf( "File could not be opened. %s\n",name);

      fprintf(IfPtr, "set nokey\nsp ",it);
  //     fprintf(IfPtr, "sp ",it);
  for(it=0;it<Nparticles-1;it++){
    fprintf(IfPtr, "\"./Trajectories/CE_trajectory_RK_%d_.dat\" using 4:5:6 w l, ",it);
  }

  fprintf(IfPtr, "\"./Trajectories/CE_trajectory_RK_%d_.dat\" using 4:5:6 w l\n pause 10.",Nparticles-1);

  fclose(IfPtr);


  for(it=0;it<Ndist;it++)
    dist[it]=0;

  sprintf(name,"out1");
  if ( ( IfPtr = fopen( name, "r" ) ) == NULL )
    printf( "File could not be opened. %s\n",name);

  while ( !feof( IfPtr ) ){
    fscanf( IfPtr, "%lf %lf %lf %lf %lf\n",&z,&z,&x,&y,&z);

    r=sqrt(x*x+y*y);
    
    it=floor((r-rmin)/(rmax-rmin)*Ndist);
    
    dist[it]++;

  }



  fclose(IfPtr);

  z=0.;
  for(it=0;it<Ndist;it++){
    z+=1.*dist[it];
    printf("%e %d %e\n",rmin+it*(rmax-rmin)/Ndist,dist[it],z);
  }

  return 0;
}

//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
