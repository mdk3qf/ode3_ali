///
/// Starter template for first baseball problem
/// Solve for the initial speed of the pitch given the initial parameters
/// xend : distance to home plate [18.5] m
/// z0 : height of release of ball [1.4] m
/// theta0 : angle of release above horizontal [1] degree
///
///  Do not change the interface for running the program
///  Fill in the value of vPitch in the print statement with your solution
///  at the end of main()
///

#include "RKn.hpp"
#include "TROOT.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TFile.h"
#include "TStyle.h"
#include "TGClient.h"
#include "TF1.h"
#include "TCanvas.h"
#include <iostream>
#include <cstdio>
#include <cstdlib>

using namespace std;

struct Params {
  double g;   // acceleration [m/s^2]
  double m;   // mass of object [kg], nb proj. In vacuum funcs do not depend on the mass
  double d;   // m diameter of ball
  double b;   // b,c params for air resistance
  double c;
};

int main(int argc, char **argv){

  // examples of parameters
  Params pars;
  pars.g=9.81;
  pars.m=0.145;    
  pars.d=0.075;   
  pars.b=1.6e-4;  
  pars.c=0.25;
  void *p_par = (void*) &pars;

  double xend=18.5;       // meters to plate
  double z0=1.4;             // height of release [m]
  double theta0=1;         // angle of velocity at release (degrees)
                                      // convert to radians before using!
  bool showPlot=false;    // keep this flag false by default
  
  // allow changing the parameters from the command line
  int c;
  while ((c = getopt (argc, argv, "x:z:t:p")) != -1)
    switch (c) {
    case 'x':
      xend = atof(optarg);
      break;
    case 'z':
      z0 = atof(optarg);
      break;
    case 't':
      theta0 = atof(optarg);
      break;
    case 'p':
      showPlot=true;
      break;
    case '?':
      fprintf (stderr, "Unknown option `%c'.\n", optopt);
    }
  TApplication theApp("App", &argc, argv); // init ROOT App for displays


  double vPitch = 0;   // m/s of pitch needed to land in strike zone at 0.9 meters
  // write code to solve for vPitch here


  // y[0] = x position (horizontal distance)
  // y[1] = vx (velocity in x direction)
  // y[2] = z position (height)
  // y[3] = vz (velocity in z direction)
  // ----------------------------------------------------------------
  // Define ODE functions for baseball motion with drag F = bv + cv^2
  
  // dx/dt = vx
  auto f_x = [](double t, const vector<double> &y, void *params) -> double {
    (void) t;
    (void) params;
    return y[1];  // vx
  };
  
  // dvx/dt = -drag_x/m
  auto f_vx = [](double t, const vector<double> &y, void *params) -> double {
    (void) t;
    Params *p = (Params*)params;
    double vx = y[1];
    double vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    if (v < 1e-10) return 0;
    double drag_coeff = p->b * p->d + p->c * p->d * p->d * v;
    return -drag_coeff * vx / p->m;
  };
  
  // dz/dt = vz
  auto f_z = [](double t, const vector<double> &y, void *params) -> double {
    (void) t;
    (void) params;
    return y[3];  // vz
  };
  
  // dvz/dt = -g - drag_z/m
  auto f_vz = [](double t, const vector<double> &y, void *params) -> double {
    (void) t;
    Params *p = (Params*)params;
    double vx = y[1];
    double vz = y[3];
    double v = sqrt(vx*vx + vz*vz);
    if (v < 1e-10) return -p->g;
    double drag_coeff = p->b * p->d + p->c * p->d * p->d * v;
    return -p->g - drag_coeff * vz / p->m;
  };
  
  // Stopping condition: stop when x >= xend
  auto f_stop = [](double t, const vector<double> &y, void *params) -> double {
    (void) t;
    (void) params;
    if (y[0] >= 18.5) return 1;
    return 0;
  };

  // define derivatives
  vector<pfunc_t> v_fun(4);
  v_fun[0] = +f_x;
  v_fun[1] = +f_vx;
  v_fun[2] = +f_z;
  v_fun[3] = +f_vz;
  
  pfunc_t stop_func = +f_stop;
  
  // Convert angle to radians
  double theta_rad = theta0 * M_PI / 180.0;
  
  // Binary search for correct initial velocity
  double v_low = 20.0;   // m/s (increased lower bound)
  double v_high = 50.0;  // m/s
  double tolerance = 0.01; // m/s
  double z_target = 0.9;  // m
  
  while (v_high - v_low > tolerance) {
    double v_mid = (v_low + v_high) / 2.0;
    
    // Set initial conditions
    vector<double> y(4);
    y[0] = 0.0;                      // x0
    y[1] = v_mid * cos(theta_rad);   // vx0
    y[2] = z0;                       // z0
    y[3] = v_mid * sin(theta_rad);   // vz0
    
    double t = 0.0;
    double tmax = 5.0;  // seconds
    int nsteps = 1000;

    // Use RK4SolveN to evolve time, updates y values according to the input functions and parameters
    auto tgN = RK4SolveN(v_fun, y, nsteps, t, tmax, p_par, stop_func);
    
    double z_final = y[2];
    
    // If ball lands too high, lower velocity
    // If ball lands too low (or hits ground), increase velocity
    if (z_final > z_target) {
      v_high = v_mid;  // Ball too high, reduce velocity
    } else {
      v_low = v_mid;   // Ball too low, increase velocity
    }
  }

  // At t_final: y = [18.5, vx_final, 0.9, vz_final]
  vPitch = (v_low + v_high)/2.0;

  // ----------------------------------------------------------------


  // do not change these lines
  printf("********************************\n");
  printf("(xend,z0,theta0) = (%lf,%lf,%lf)\n",xend,z0,theta0);
  printf("v_pitch = %lf m/s\n",vPitch);
  printf("********************************\n");

  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}