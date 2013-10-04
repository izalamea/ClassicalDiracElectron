/**

This code solves the time evolution for a Clasical Dirac Electron.

Time integrations uses 4th order Runge-Kutta

The main theory behind is in the paper:

A. O. Barut, Nino Zanghi
Phys. Rev. Lett. 52, 2009–2012 (1984)
Classical Model of the Dirac Electron

**/
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

const int run_name=2;
const double lambda=1.;
const double q=1.;
const double QC=.0;
const double EX=0.00;
const double EZ=0.0;
const double BZ=-.9;
const double dt0=1.e-4;
const double T=30.*M_PI;

void XDOT(double *zr,double *zi, double *xdot);
void ZDOT(double *zr,double *zi, double *p, double *zrdot, double *zidot);
void PDOT(double **Fmunu, double *xdot, double *pdot);
void ConstantEB(double **Fmunu);
void OneOverr2E(double **Fmunu, double *x);

void assingAtoB(double *A, double *B);//B=A 
double *allocate_1(int m);
double **allocate_2(int m, int n);
void freemn(double **vn, int m);

int main()
{
    double **Fmunu;
    double *x0, *x1, *zr1, *zr0, *zi1, *zi0, *p1, *p0;
    double *x2, *x3, *zr2, *zr3, *zi2, *zi3, *p2, *p3;
    double  *xdot, *zrdot, *zidot, *pdot;
    double *xdot1, *zrdot1, *zidot1, *pdot1;
    double *xdot2, *zrdot2, *zidot2, *pdot2;
    double *xdot3, *zrdot3, *zidot3, *pdot3;
    
    double alpha, beta,t=0., H=0.,dt=dt0,v0, theta[7];
    char name[400];
    int i=0;
    FILE *IfPtr;
    
    sprintf(name,"CE_trajectory_RK_%d_.dat",run_name);
    if ( ( IfPtr = fopen( name, "w" ) ) == NULL )
        printf( "File could not be opened. %s\n",name);

    
    Fmunu=allocate_2(4,4);
    x0=allocate_1(4);
    x1=allocate_1(4);
    x2=allocate_1(4);
    x3=allocate_1(4);
    
    zr0=allocate_1(4);
    zr1=allocate_1(4);
    zr2=allocate_1(4);
    zr3=allocate_1(4);
    
    zi0=allocate_1(4);
    zi1=allocate_1(4);
    zi2=allocate_1(4);
    zi3=allocate_1(4);
    
    p0=allocate_1(4);
    p1=allocate_1(4);
    p2=allocate_1(4);
    p3=allocate_1(4);

    xdot=allocate_1(4);
    xdot1=allocate_1(4);
    xdot2=allocate_1(4);
    xdot3=allocate_1(4);
    
    zrdot=allocate_1(4);
    zrdot1=allocate_1(4);
    zrdot2=allocate_1(4);    
    zrdot3=allocate_1(4);
    
    zidot=allocate_1(4);
    zidot1=allocate_1(4);
    zidot2=allocate_1(4);
    zidot3=allocate_1(4);
    
    pdot=allocate_1(4);
    pdot1=allocate_1(4);
    pdot2=allocate_1(4);
    pdot3=allocate_1(4);
    
    //Initialize EM filed and charge state.
    ConstantEB(Fmunu);
    
    x0[0]=0.;
    x0[1]=0.; 
    x0[2]=0.;
    x0[3]=0.;

    theta[0]=0.25*M_PI;

    theta[1]=0.25*M_PI;
    theta[2]=0.5*M_PI;
    theta[3]=0.5*M_PI;

    theta[4]=0.5*M_PI;
    theta[5]=0.25*M_PI;
    theta[6]=0.*M_PI;
    /**
    zr0[0]=cos(theta[0])*cos(theta[1]);
    zr0[1]=cos(theta[0])*sin(theta[1])*cos(theta[2]); 
    zr0[2]=cos(theta[0])*sin(theta[1])*sin(theta[2])*cos(theta[3]);
    zr0[3]=cos(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[3]);
    
    zi0[0]=sin(theta[0])*cos(theta[4]);
    zi0[1]=sin(theta[0])*sin(theta[4])*cos(theta[5]);
    zi0[2]=sin(theta[0])*sin(theta[4])*sin(theta[5])*cos(theta[6]);
    zi0[3]=sin(theta[0])*sin(theta[4])*sin(theta[5])*sin(theta[6]);

    /**/
    zr0[0]=1.;
    zr0[1]=0.; 
    zr0[2]=0.;
    zr0[3]=-zr0[0];
    
    theta[0]=-1.*M_PI;
    zi0[0]=0.;
    zi0[1]=0.*cos(theta[0]); 
    zi0[2]=0.*sin(theta[0]);
    zi0[3]=0.;
    /**/
    alpha=1.*0.5*M_PI;
    beta=0.*2.*M_PI;
    p0[0]=0.;

    p0[1]=p0[0]*cos(alpha); 
    p0[2]=p0[0]*sin(alpha)*sin(beta);
    p0[3]=p0[0]*sin(alpha)*cos(beta);
    p0[0]=sqrt(1.+p0[1]*p0[1]+p0[2]*p0[2]+p0[3]*p0[3]);


    while(t<T && p0[0]*p0[0]-p0[1]*p0[1]-p0[2]*p0[2]-p0[3]*p0[3]>-0.1){
        assingAtoB(zr0,zr1);
        assingAtoB(zi0,zi1);
        assingAtoB(p0,p1);
        assingAtoB(x0,x1);

        assingAtoB(zr0,zr2);
        assingAtoB(zi0,zi2);
        assingAtoB(p0,p2);
        assingAtoB(x0,x2);

        assingAtoB(zr0,zr3);
        assingAtoB(zi0,zi3);
        assingAtoB(p0,p3);
        assingAtoB(x0,x3);

        XDOT(zr0,zi0,xdot);
        ZDOT(zr0,zi0,p0,zrdot,zidot);/**
	OneOverr2E(Fmunu,x0);/**/
        PDOT(Fmunu,xdot,pdot);

	v0=sqrt(xdot[1]*xdot[1]+xdot[2]*xdot[2]+xdot[3]*xdot[3]);

	//	if(v0*dt>1.e-3)
	//	  dt=1.e-3/v0;
	//	else dt=dt0;

        H=-p0[0]*xdot[0]+p0[1]*xdot[1]+p0[2]*xdot[2]+p0[3]*xdot[3];
        alpha=zr0[0]*zr0[0]+zr0[1]*zr0[1]-zr0[2]*zr0[2]-zr0[3]*zr0[3];
        alpha+=zi0[0]*zi0[0]+zi0[1]*zi0[1]-zi0[2]*zi0[2]-zi0[3]*zi0[3];
        
        if(i%20==0){
            fprintf(IfPtr,"%e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e %e\n",t,H,x0[0],x0[1],x0[2],x0[3],p0[0],p0[1],p0[2],p0[3],zr0[0],zr0[1],zr0[2],zr0[3],zi0[0],zi0[1],zi0[2],zi0[3],xdot[0],xdot[1],xdot[2],xdot[3],alpha);
            printf("%e %e %e %e %e\n",t,x0[0],x0[1],x0[2],x0[3]);
        }
/**/        
        x1[0]+=0.5*dt*xdot[0];
        x1[1]+=0.5*dt*xdot[1]; 
        x1[2]+=0.5*dt*xdot[2];
        x1[3]+=0.5*dt*xdot[3];
        
        zr1[0]+=0.5*dt*zrdot[0];
        zr1[1]+=0.5*dt*zrdot[1]; 
        zr1[2]+=0.5*dt*zrdot[2];
        zr1[3]+=0.5*dt*zrdot[3];
        
        zi1[0]+=0.5*dt*zidot[0];
        zi1[1]+=0.5*dt*zidot[1]; 
        zi1[2]+=0.5*dt*zidot[2];
        zi1[3]+=0.5*dt*zidot[3];
        
        p1[0]+=0.5*dt*pdot[0];
        p1[1]+=0.5*dt*pdot[1]; 
        p1[2]+=0.5*dt*pdot[2];
        p1[3]+=0.5*dt*pdot[3];

        
        XDOT(zr1,zi1,xdot1);
        ZDOT(zr1,zi1,p1,zrdot1,zidot1);/**
	OneOverr2E(Fmunu,x1);/**/
        PDOT(Fmunu,xdot1,pdot1);

        
        x2[0]+=0.5*dt*xdot1[0];
        x2[1]+=0.5*dt*xdot1[1]; 
        x2[2]+=0.5*dt*xdot1[2];
        x2[3]+=0.5*dt*xdot1[3];
        
        zr2[0]+=0.5*dt*zrdot1[0];
        zr2[1]+=0.5*dt*zrdot1[1]; 
        zr2[2]+=0.5*dt*zrdot1[2];
        zr2[3]+=0.5*dt*zrdot1[3];
        
        zi2[0]+=0.5*dt*zidot1[0];
        zi2[1]+=0.5*dt*zidot1[1]; 
        zi2[2]+=0.5*dt*zidot1[2];
        zi2[3]+=0.5*dt*zidot1[3];
        
        p2[0]+=0.5*dt*pdot1[0];
        p2[1]+=0.5*dt*pdot1[1]; 
        p2[2]+=0.5*dt*pdot1[2];
        p2[3]+=0.5*dt*pdot1[3];        

        XDOT(zr2,zi2,xdot2);
        ZDOT(zr2,zi2,p2,zrdot2,zidot2);/**
	OneOverr2E(Fmunu,x2);/**/
        PDOT(Fmunu,xdot2,pdot2);

        x3[0]+=dt*xdot2[0];
        x3[1]+=dt*xdot2[1]; 
        x3[2]+=dt*xdot2[2];
        x3[3]+=dt*xdot2[3];
        
        zr3[0]+=dt*zrdot2[0];
        zr3[1]+=dt*zrdot2[1]; 
        zr3[2]+=dt*zrdot2[2];
        zr3[3]+=dt*zrdot2[3];
        
        zi3[0]+=dt*zidot2[0];
        zi3[1]+=dt*zidot2[1]; 
        zi3[2]+=dt*zidot2[2];
        zi3[3]+=dt*zidot2[3];
        
        p3[0]+=dt*pdot2[0];
        p3[1]+=dt*pdot2[1]; 
        p3[2]+=dt*pdot2[2];
        p3[3]+=dt*pdot2[3];        
        
        XDOT(zr3,zi3,xdot3);
        ZDOT(zr3,zi3,p3,zrdot3,zidot3);/**
        OneOverr2E(Fmunu,x3);/**/
        PDOT(Fmunu,xdot3,pdot3);
        
        /**/        
        x0[0]+=dt*(xdot[0]+2.*xdot1[0]+2.*xdot2[0]+xdot3[0])/6.;
        x0[1]+=dt*(xdot[1]+2.*xdot1[1]+2.*xdot2[1]+xdot3[1])/6.; 
        x0[2]+=dt*(xdot[2]+2.*xdot1[2]+2.*xdot2[2]+xdot3[2])/6.;
        x0[3]+=dt*(xdot[3]+2.*xdot1[3]+2.*xdot2[3]+xdot3[3])/6.;
        
        zr0[0]+=dt*(zrdot[0]+2.*zrdot1[0]+2.*zrdot2[0]+zrdot3[0])/6.;
        zr0[1]+=dt*(zrdot[1]+2.*zrdot1[1]+2.*zrdot2[1]+zrdot3[1])/6.; 
        zr0[2]+=dt*(zrdot[2]+2.*zrdot1[2]+2.*zrdot2[2]+zrdot3[2])/6.;
        zr0[3]+=dt*(zrdot[3]+2.*zrdot1[3]+2.*zrdot2[3]+zrdot3[3])/6.;

        zi0[0]+=dt*(zidot[0]+2.*zidot1[0]+2.*zidot2[0]+zidot3[0])/6.;
        zi0[1]+=dt*(zidot[1]+2.*zidot1[1]+2.*zidot2[1]+zidot3[1])/6.; 
        zi0[2]+=dt*(zidot[2]+2.*zidot1[2]+2.*zidot2[2]+zidot3[2])/6.;
        zi0[3]+=dt*(zidot[3]+2.*zidot1[3]+2.*zidot2[3]+zidot3[3])/6.;
                
        p0[0]+=dt*(pdot[0]+2.*pdot1[0]+2.*pdot2[0]+pdot3[0])/6.;
        p0[1]+=dt*(pdot[1]+2.*pdot1[1]+2.*pdot2[1]+pdot3[1])/6.; 
        p0[2]+=dt*(pdot[2]+2.*pdot1[2]+2.*pdot2[2]+pdot3[2])/6.;
        p0[3]+=dt*(pdot[3]+2.*pdot1[3]+2.*pdot2[3]+pdot3[3])/6.;
        
        
        t+=dt;
        i++;
    }
    
    freemn(Fmunu,4);
    free(Fmunu);
    free(x0);
    free(x1);    
    free(x2);
    free(x3);

    free(zr0);
    free(zr1);
    free(zr2);
    free(zr3);
    
    free(zi0);
    free(zi1);    
    free(zi2);
    free(zi3);    
    
    free(p0);
    free(p1);
    free(p2);
    free(p3);    
    
    free(xdot);
    free(xdot1);
    free(xdot2);
    free(xdot3);
    
    free(pdot);
    free(pdot1);
    free(pdot2);
    free(pdot3);
    
    free(zrdot);
    free(zrdot1);
    free(zrdot2);
    free(zrdot3);
    
    free(zidot);
    free(zidot1);    
    free(zidot2);
    free(zidot3);    
    
    fclose(IfPtr);
    
    return 0;
}

/*************************************************************************/

double *allocate_1(int m) 
{
double *a;

a = (double *) malloc((unsigned) (m)*sizeof(double));

return a;
}

/************************************************************************************************/

double **allocate_2(int m, int n) 
{
    int i;
    double **a;
    
    a = (double **) malloc((unsigned) (m)*sizeof(double*));
    
    for (i=0; i<m; i++) 
	{
        a[i] = (double *) malloc((unsigned) n*sizeof(double));
	}
    return a;
}

/*************************************************************************/


void freemn(double **vn, int m)
{
    int i;
    
    for(i=0;i<m;i++)
        free(vn[i]);
    return;
}

/*************************************************************************/

void assingAtoB(double *A, double *B)//B=A 
{
    B[0]=A[0];
    B[1]=A[1];
    B[2]=A[2];
    B[3]=A[3];
    
    return;
}

/************************************************************************************************/


void XDOT(double *zr,double *zi, double *xdot)
{
    
    xdot[0]=zr[0]*zr[0]+zr[1]*zr[1]+zr[2]*zr[2]+zr[3]*zr[3]+zi[0]*zi[0]+zi[1]*zi[1]+zi[2]*zi[2]+zi[3]*zi[3];

    xdot[1]=2.*( zr[0]*zr[3] +zi[0]*zi[3] +zr[1]*zr[2] +zi[1]*zi[2]);

    xdot[2]=2.*( zr[0]*zi[3] -zi[0]*zr[3] -zr[1]*zi[2] +zi[1]*zr[2]);

    xdot[3]=2.*( zr[0]*zr[2] +zi[0]*zi[2] -zr[1]*zr[3] -zi[1]*zi[3]);
    
    return;
}

/*************************************************************************/


void ZDOT(double *zr,double *zi, double *p, double *zrdot, double *zidot)
{
    
    zrdot[0]=(-p[0]*zi[0] + p[1]*zi[3] - p[2]*zr[3] + p[3]*zi[2])/lambda;

    zrdot[1]=(-p[0]*zi[1] + p[1]*zi[2] + p[2]*zr[2] - p[3]*zi[3])/lambda;
    
    zrdot[2]=( p[0]*zi[2] - p[1]*zi[1] + p[2]*zr[1] - p[3]*zi[0])/lambda;
    
    zrdot[3]=( p[0]*zi[3] - p[1]*zi[0] - p[2]*zr[0] + p[3]*zi[1])/lambda;
    
    
    zidot[0]=(-p[0]*zr[0] + p[1]*zr[3] + p[2]*zi[3] + p[3]*zr[2])/(-lambda);
    
    zidot[1]=(-p[0]*zr[1] + p[1]*zr[2] - p[2]*zi[2] - p[3]*zr[3])/(-lambda);
    
    zidot[2]=( p[0]*zr[2] - p[1]*zr[1] - p[2]*zi[1] - p[3]*zr[0])/(-lambda);
    
    zidot[3]=( p[0]*zr[3] - p[1]*zr[0] + p[2]*zi[0] + p[3]*zr[1])/(-lambda);
    
    return;
}

/*************************************************************************/


void PDOT(double **Fmunu, double *xdot, double *pdot)
{
    pdot[0]=q*( Fmunu[0][0]*xdot[0] + Fmunu[0][1]*xdot[1] + Fmunu[0][2]*xdot[2] + Fmunu[0][3]*xdot[3]);
    
    pdot[1]=-q*( Fmunu[1][0]*xdot[0] + Fmunu[1][1]*xdot[1] + Fmunu[1][2]*xdot[2] + Fmunu[1][3]*xdot[3]);

    pdot[2]=-q*( Fmunu[2][0]*xdot[0] + Fmunu[2][1]*xdot[1] + Fmunu[2][2]*xdot[2] + Fmunu[2][3]*xdot[3]);
    
    pdot[3]=-q*( Fmunu[3][0]*xdot[0] + Fmunu[3][1]*xdot[1] + Fmunu[3][2]*xdot[2] + Fmunu[3][3]*xdot[3]);
    
    return;
}


/*************************************************************************/


void ConstantEB(double **Fmunu)
{
    Fmunu[0][0]=0.;   Fmunu[0][1]=EX;   Fmunu[0][2]=0.;   Fmunu[0][3]=EZ;
    
    Fmunu[1][0]=-EX;   Fmunu[1][1]=0.;   Fmunu[1][2]=-BZ;   Fmunu[1][3]=0.;
    
    Fmunu[2][0]=0.;   Fmunu[2][1]=BZ;   Fmunu[2][2]=0.;   Fmunu[2][3]=0.;
    
    Fmunu[3][0]=-EZ;   Fmunu[3][1]=0.;   Fmunu[3][2]=0.;   Fmunu[3][3]=0.;
    
    return;
}

/*************************************************************************/


void OneOverr2E(double **Fmunu, double *x)
{   double r3=pow(x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+1.e-2,-1.5);
    double Ex=QC*x[1]*r3, Ey=QC*x[2]*r3, Ez=QC*x[3]*r3;
    
    Fmunu[0][0]=0.;   Fmunu[0][1]=Ex;   Fmunu[0][2]=Ey;   Fmunu[0][3]=Ez;
    
    Fmunu[1][0]=-Ex;   Fmunu[1][1]=0.;   Fmunu[1][2]=0.;   Fmunu[1][3]=0.;
    
    Fmunu[2][0]=-Ey;   Fmunu[2][1]=0.;   Fmunu[2][2]=0.;   Fmunu[2][3]=0.;
    
    Fmunu[3][0]=-Ez;   Fmunu[3][1]=0.;   Fmunu[3][2]=0.;   Fmunu[3][3]=0.;
    
    return;
}

/*************************************************************************/

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
//
