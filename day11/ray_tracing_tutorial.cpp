#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>

using namespace std;
using std::ifstream;

#define SHADOW (0)
#define DISK (1)

//edit 0: which type of problem?
//#define whichstop SHADOW
#define whichstop DISK

#define PI 3.14159265358979323846
#define N 6
//resolution
#define m 256
int M=m*m;

#define t0 0
#define r0 1000.0
//edit 1: field of view
#define fov 0.05
//#define fov 0.025
//edit 2: inclination
double theta0=75.*PI/180.0;
//double theta0=60.*PI/180.0;
//double theta0=20.*PI/180.0;
#define phi0 0
//black hole spin
#define a 0.9
double a2 = a*a;
double RHor=1+sqrt(1-a2)+0.00001;
#define PI 3.14159265358979323846
#define RDisc 20.0

double L;
double kappa;
double RMstable;

void initial(double *y0, double *ydot0, double x, double y)
{
      y0[0]=r0;
      y0[1]=theta0;
      y0[2]=phi0;
      y0[3]=t0;
      y0[4]=-cos(y)*cos(x);
      y0[5]=sin(y)/r0;
      double rdot0=y0[4];
      double thetadot0=y0[5];
      
      double r2=r0*r0;
      double costheta2=cos(theta0)*cos(theta0);
      double sintheta2=sin(theta0)*sin(theta0);
      double sum=r2+a2*costheta2;
      double delta=r2-2.0*r0+a2;

      y0[4]= rdot0*sum/delta;
      y0[5]= thetadot0*sum;

      ydot0[0] = rdot0;
      ydot0[1] = thetadot0;
      ydot0[2] = cos(y)*sin(x)/(r0*sin(theta0));
      
      double phidot0=ydot0[2];
      double energy2=(sum-2.0*r0)*(rdot0*rdot0/delta+thetadot0*thetadot0)+delta*sintheta2*phidot0*phidot0;

      double energy = sqrt(energy2);
    
      L=(sum*delta*phidot0-2.0*a*r0*energy)*sintheta2/(sum-2.0*r0);


      y0[4]=y0[4]/energy;
      y0[5]=y0[5]/energy;     
      L=L/energy;
      kappa=y0[5]*y0[5]+a2*sintheta2+L*L/sintheta2;
} 



void geodesic(double *y, double *dydlamda)
{
  double r, theta, phi, t, pr, ptheta;
  
  r=y[0];
  theta= y[1];
  phi = y[2];
  t = y[3];
  pr = y[4];
  ptheta = y[5];

  //if(r>1000.0) r=0.001;
  double r2=r*r;
  double twor=2.0*r;
  double cos2=cos(theta)*cos(theta);
  double sin2=sin(theta)*sin(theta);
  double sum=r2+a2*cos2;
  double delta=r2-twor+a2;
  double sd=sum*delta;

  dydlamda[0] = pr*delta/sum;
  dydlamda[1] = ptheta/sum;
  dydlamda[2] = (twor*a+(sum-twor)*L/sin2)/sd;
  dydlamda[3] = 1.0+(twor*(r2+a2)-twor*a*L)/sd;
  dydlamda[4] = ((r-1.0)*(-kappa)+twor*(r2+a2)-2.0*a*L)/sd-2.0*pr*pr*(r-1.0)/sum;
  dydlamda[5] = sin(theta)*cos(theta)*(L*L/(sin2*sin2)-a2)/sum;

}

void rkck (double *y, double *dydx, double h, double *yout, double *yerr)
{
  static const double b21 = 0.2, b31 = 3.0/40.0, b32 = 9.0/40.0, b41 = 0.3, b42 = -0.9,
    b43 = 1.2, b51 = -11.0/54.0, b52 = 2.5, b53 = -70.0/27.0,
    b54 = 35.0/27.0, b61 = 1631.0/55296.0, b62 = 175.0/512.0,
    b63 = 575.0/13824.0, b64 = 44275.0/110592.0, b65 = 253.0/4096.0,
    c1 = 37.0/378.0, c3 = 250.0/621.0, c4 = 125.0/594.0, c6 = 512.0/1771.0,
    dc1 = c1-2825.0/27648.0, dc3 = c3-18575.0/48384.0,
    dc4 = c4-13525.0/55296.0, dc5 = -277.00/14336.0, dc6 = c6-0.25;

  
  int i;
  int n=N;

  double ak2[6], ak3[6], ak4[6], ak5[6], ak6[6], ytemp[6];


  for (i=0;i<n;i++) ytemp[i]=y[i]+b21*h*dydx[i];
  geodesic(ytemp, ak2);
  
  
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
  geodesic(ytemp, ak3);
 
  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
  geodesic(ytemp, ak4);

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  geodesic(ytemp, ak5);  

  for (i=0;i<n;i++) ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  geodesic(ytemp, ak6); 

  for (i=0;i<n;i++) yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
 
  for (i=0;i<n;i++) yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);

}

double SIGN(double x)
{
  double val;
  if (x>=0.) val=1.;
  else val=-1.;
  return val;
}

double RMS(double spin)
{
  double Z1=1.+pow((1.-spin*spin),1./3.)*(pow(1.+spin,1./3.)+pow(1.-spin,1./3.));
  double Z2=sqrt(3.*spin*spin+Z1*Z1);
  double rms = 3.+Z2-SIGN(spin)*sqrt((3.-Z1)*(3.+Z1+2.*Z2));
  return rms;
}

double MAX(double x, double y)
{
    if (x >= y)
      return x;
    else 
      return y;
}

double MIN(double x, double y)
{
    if (x <= y)
      return x;
    else
      return y;
}

void nrerror(const string error_text)
{
  cerr << "Numerical Recipes run-time error..." << endl;
  cerr << error_text << endl;
  cerr << "...now exiting to system..." << endl;
  exit(1);
}

#include <cmath>

void rkqs (double *y, double *dydx, double &x, double htry, double eps, double *yscal, double &hdid, double &hnext)
{
  const double SAFETY=0.9, PGROW=-0.2, PSHRNK=-0.25, ERRCON=1.89e-4;
  
  int i;

  double errmax, h, htemp, xnew;

  int n=N;
  h=htry;
  double yerr[6], ytemp[6];


  for (;;) 
    {   
      rkck(y, dydx, h, ytemp, yerr);

      errmax=0.0;
      for(i=0;i<n;i++) errmax = MAX(errmax, fabs(yerr[i]/yscal[i]));
      //      cout << "errmax " << errmax <<"\n";
      //      cout << "rkqs errstep: " << htry << "\t" << yerr[0] << "\n";

      errmax/=eps;
      if(errmax<=1.0) break;

      htemp=SAFETY*h*pow(errmax, PSHRNK);
      h=(h>=0.0 ? MAX(htemp,0.1*h) : MIN(htemp, 0.1*h));
    
      xnew=x+h;
      if(xnew==x) nrerror("stepsize underflow in rkqs");
    }

  if (errmax>ERRCON) hnext=SAFETY*h*pow(errmax, PGROW);
  else hnext=5.0*h;

  x+=(hdid=h);

  for (i=0;i<n;i++) y[i]=ytemp[i];

}

void binarysearch (double *y, double *dydx, double &hsmall, double &hbig)
{
  int q=0;
  double side;
  if (y[1]>PI/2) side=1.0;
  else if (y[1]<PI/2) side=-1.0;
  else side=0.0;

  geodesic(y,dydx);

   while(y[0]>RHor && y[0]<1001 && side!=0.0)
   {
	if (hbig < hsmall)
	{
		cout <<"left bigger than right" <<"\n";
       		break;
	}

	int i;
	int n = N;

	double yout[6], yerr[6]; 
 
	double difference = hbig-hsmall;

        double hmid = difference/2.0+hsmall;

        rkck(y, dydx, hmid, yout, yerr);

	q++;

        if (q>20000) break;

        if (difference<0.0000001)
	{
	  //cout << "Done!"<< "\n";
	  hbig=hmid;
          for(i=0;i<n;i++) 
	    {
	     y[i]=yout[i];
	    }
	  // if (y[0]>RDisc) y[0]=10000001;
	  if (side>0)y[1]=PI/2.0-0.0000001;
	  else y[1]=PI/2.0+0.0000001;

	     //printf("y[%d] is %e\n", i, y[i]); 
	      //cout<< "y[" << i <<"] is "<<y[i]<<endl;
	  break;
	} 

	if (side*(yout[1]-PI/2.0)>0)
	{
	  hsmall=hmid;
	}

	else
	{ 
	  hbig=hmid;
	}
   }
}

int main()
{
  int n = N;
  double y[6], dydx[6], yscal[6];

  double SIDE;
  double ylaststep[6];
  double xlaststep=0;

  std::ofstream datafile;
  datafile.open("position.dat");
  
  datafile<<m<<" "<<a<<endl;
  // marginally stable orbit
  RMstable = RMS(a);
  //cout << "rms: " << RMstable << endl;
  
  for(int j=0; j<M; j++)

    {

      int order=0;

      double htry=0.5, eps=1e-11, hdid=0.0, hnext=0.0;
      int i=0, q=0, t_prev=0;
      double intensity3=0.;
      const double TINY=1.0e-8;
      double x=0.0;
 
      // create 2D grid camera
      double range = fov/double(m-1);
      double p, t;
      t = -(j/m-m/2+0.5)*range;

      if(j/m != m) p = (j%m-m/2+0.5)*range;
      else p = (0.5-m/2)*range;

      double alpha = t*r0;
      double beta = p*r0;

      //cout <<"p is "<< p <<"; t is " << t << endl;

      // initial position and momentum for each ray throug the grid
      initial(y, dydx, p, t);
      
      while (y[0]<1001&&y[0]>0)     // go till the particle falls into the BH or fly faraway
	{
	  for (i=0;i<n; i++) ylaststep[i]=y[i];
	  xlaststep=x; 
       
      // calculate the geodesic until reaching the black hole, the disk, or infinity 	
	  geodesic(y, dydx);
       
	  for (i=0; i<n; i++) yscal[i]=fabs(y[i])+fabs(dydx[i]*htry)+TINY;
       
	  if (y[1]>PI/2) SIDE= 1.0;
	  else if (y[1]<PI/2) SIDE = -1.0;
	  else SIDE = 0.0;
	

	  rkqs(y, dydx, x, htry, eps, yscal, hdid, hnext);
	  //	  cout << "hnext htry: " << htry << "\t" << hnext << "\n";
	  q++;
	  if (whichstop==DISK)
	    {
	      if(y[0]<RHor)
		{
		  y[0]=100001;
		  break;
		}
	      else if ((y[1]-PI/2)*SIDE<0)
		{
		  if(y[0]<(RDisc+4)&&y[0]>=(RMstable-0.5))
		    {
		      for (i=0;i<n;i++) y[i]=ylaststep[i];
		      x=xlaststep;
	      	
		      double hbig=hdid;
		      double hsmall=0.0;
		      binarysearch(y, dydx, hsmall, hbig);

		      if (y[0]<=RDisc&&y[0]>=RMstable)
			{
			  break;
			}
		      else order++;
		    }
		  else order++;
		}
	      
	      else if(SIDE==0.0&&y[0]<=RDisc&&y[0]>=RMstable)
		{
		  break;
		}
	      else if (q>1000000) 
		{
		  break;
		}
        
	      htry=hnext;
	    }
	
      // stop when we reach the horizon, large radius, or take too long
	  else if (whichstop==SHADOW)
	    {
	      double stheta0=sin(theta0);
	      double VPL=0.0;
	      if (fabs(cos(y[1])) < 0.98) VPL=0.2;
	      double D=(y[0]*y[0]-2.*y[0]+a*a);
	      double AR=(y[0]*y[0]+a*a)*(y[0]*y[0]+a*a)-a*a*D*stheta0*stheta0;
	      double RHO=y[0]*y[0]+a*a*cos(theta0);
	      double OM=2.*a*y[0]/AR;
	      double ENU=sqrt(D*RHO/AR);
	      double EPSI=sin(theta0)*sqrt(AR/RHO);
	      double Omega=ENU/EPSI*VPL+OM;
	      double GAM=1./sqrt(1.-VPL*VPL);
	      double G=ENU/GAM/(1.-L*Omega);
	      double j=1./y[0]/y[0]*exp(-cos(y[1])*cos(y[1]));
	      intensity3+=(y[3]-t_prev)*j*G*G;
	      //	      cout << "intensity3: " << intensity3 << " " << y[3] << " " << t_prev << " " << j << " " << G << endl;
	      if (y[0] < RHor)
		{
		  break;
		}
	      else if (q>1000000)
		{
		  break;
		}
	      htry=hnext;
	      t_prev=y[3];
	    }
	}
      double ut, uphi, xi,pdotu;
      double g, intensity, intensity2;

      if (y[0]>RDisc)
	{
	  g=100; 
	  intensity=0;
	  intensity2=0;
	}

      else
	{
	  // GR Keplerian speed for thin disk
	  //	  double yr=MAX(y[0],3.);
	  double yr=y[0];
	  double r2=yr*yr;
	  double r32=yr*sqrt(yr);
	  double rp=1.;
	  if (yr >= RMstable)
	    {
	      xi=sqrt((r2*yr-3.0*r2+2.0*a*r32));
	      ut=(a+r32)/xi;
	      uphi=1.0/xi;
	      pdotu=-ut+L*uphi;
     	      g=-1.0/pdotu;
	    }
	  
	  double g2 =g*g;
	  // cout << "g: " << g << " "<<ut << " "<< uphi << endl;
	  // assume intensity ~ 1/r^2 * g^4
	  intensity=fabs(g2*g/r2);
	  intensity2=fabs(g2*g/rp);
	}
	  // output for the ray through each grid: r, phi, order of the ring, redshift, intensity
      //      cout << "intensity3 write: " << intensity3 << endl;
      datafile << alpha << " " << beta << " " << y[0] << " " <<y[2] << " " << order << " " << g <<" "<<intensity<<" "<<intensity2<<" "<<y[3]<<" "<<intensity3<<endl;
      }

  datafile.close();
  return 0;
}

