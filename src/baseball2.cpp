///
/// Starter template for second baseball problem
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


int main(int argc, char **argv){

  // we have 6 initial conditions for this problem
  // y[0] = y[2] = y[4] = 0;  // init x,y,z
  // y[1] = v0*cos(theta0);   // vx  "x is line towards the plate
  // y[3] = 0;                // vy  "y" is measured as left/right divergence from line to plate
  // y[5] = v0*sin(theta0);   // vz  "z" is vertival measure
  vector<double> y0(6);

  bool showPlot=false;
  // pitches
  // slider ip=0
  // curve ip=1
  // screwball ip=2
  // fast ip=3
  int ip=1;    // default pitch
  int c;
  while ((c = getopt (argc, argv, "p:n")) != -1)
    switch (c) {
    case 'p':
      ip = atoi(optarg);
      break;
    case 'n':
      showPlot=false;
      break;
    }

  TString title;
  if (ip==0){
    cout << "Setting up initial conditions for slider" << endl;
    //SetupSlider(y0);
  }
  else if (ip==1){
    cout << "Setting up initial conditions for curveball" << endl;
    //SetupCurve(y0);
  }
  else if (ip==2){
    cout << "Setting up initial conditions for screwball" << endl;
    //SetupScrewball(y0);
  }
  else {
    cout << "Setting up initial conditions for fastball" << endl;
    //SetupFastball(y0);
  }

  TApplication theApp("App", &argc, argv); // init ROOT App for displays

  double xend=60;   // feet
  double yend=0;    // tbd
  double zend=0;    // tbd
  double vxend=0;
  double vyend=0;
  double vzend=0;

  // ---------------------------------------------------------------------------
  // WRITE CODE HERE

  // 3D model of baseball with air resistance and magnus effect
      // ω = ω (0, sin φ, cos φ) where φ is angle relative to z axis
      // F_mag_x = m * B * ω * (vz * sin φ - vy * cos φ)
      // F_mag_y = m * B * ω * (vx * cos φ)
      // F_mag_z =  -m * B * ω * (vx * sin φ)
      // F_drag_i =  -m * F * v *v_i
  
    // For a baseball approximate:
      // B = S/m = 4.1×10−4
      // F = 0.0039 + (0.0058)/[1+exp(0.2(v-35))]

  // Define parameters struct
  struct Params {
    double g;      // gravity [m/s^2]
    double B;      // Magnus coefficient
    double omega;  // spin rate [rad/s]
    double phi;    // spin axis angle relative to z axis [rad]
  };

  Params pars; // SI units
  pars.g = 9.81;  
  pars.B = 4.1e-4;
  
  double v0, omega_rpm, phi;

  // pitch parameters
  if (ip == 0) {              // Slider
      v0 = 38.0;              // m/s  (≈85 mph)
      omega_rpm = 1700.0;
      phi = 0.0;              // spin axis along z
  }
  else if (ip == 1) {         // Curveball
      v0 = 33.5;              // m/s  (≈75 mph)
      omega_rpm = 2000.0;
      phi = M_PI/4.0;         // spin axis 45 degrees from z
  }
  else if (ip == 2) {         // Screwball
      v0 = 31.3;              // m/s
      omega_rpm = 1800.0;
      phi = 3.0 * M_PI/4.0;   // spin axis 135 degrees from z
  }
  else {                      // Fastball
      v0 = 40.2;              // m/s (≈90 mph)
      omega_rpm = 1500.0;
      phi = M_PI/2.0;         // spin axis 90 degrees from z
  }

  // convert rpm to rad/s
  pars.omega = omega_rpm * 2.0 * M_PI / 60.0;
  
  // set spin axis angle in pars
  pars.phi = phi;
  
  // Initial conditions
  double theta = 0.0;               // release angle
  y0[0] = 0.0;                      // x0
  y0[1] = v0 * cos(theta);          // vx0
  y0[2] = 0.0;                      // y0
  y0[3] = 0.0;                      // vy0
  y0[4] = 1.8;                      // z0 (6 feet in meters)
  y0[5] = v0 * sin(theta);          // vz0
  
  void *p_par = (void*)&pars;
      

  // Define ODE functions
  auto f_x = [](double t, const vector<double> &y, void *params) -> double {
    return y[1]; // dx/dt = vx
  };


  auto f_vx = [](double t, const vector<double> &y, void *params) -> double {
    Params *p = (Params*)params;
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-10) return 0;
    
    // Drag coefficient F(v)
    double F = 0.0039 + 0.0058 / (1.0 + exp(0.2 * (v - 35.0)));
    
    // Drag force
    double drag_x = -F * v * vx;
    
    // Magnus force in x direction
    double mag_x = p->B * p->omega * (vz * sin(p->phi) - vy * cos(p->phi));
    
    return drag_x + mag_x;
  };

  auto f_y = [](double t, const vector<double> &y, void *params) -> double {
    return y[3]; // dy/dt = vy
  };


  auto f_vy = [](double t, const vector<double> &y, void *params) -> double {
    Params *p = (Params*)params;
    double vx = y[1];
    double vy = y[3];
    double vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-10) return 0;
    
    // Drag coefficient F(v)
    double F = 0.0039 + 0.0058 / (1.0 + exp(0.2 * (v - 35.0)));
    
    // Drag force
    double drag_y = -F * v * vy;
    
    // Magnus force in y direction
    double mag_y = p->B * p->omega * vx * cos(p->phi);
    
    return drag_y + mag_y;
  };


  auto f_z = [](double t, const vector<double> &y, void *params) -> double {
    return y[5]; // dz/dt = vz
  };


  auto f_vz = [](double t, const vector<double> &y, void *params) -> double {
    Params *p = (Params*)params;
    double vx = y[1], vy = y[3], vz = y[5];
    double v = sqrt(vx*vx + vy*vy + vz*vz);
    if (v < 1e-10) return -p->g;
    
    // Drag coefficient F(v)
    double F = 0.0039 + 0.0058 / (1.0 + exp(0.2 * (v - 35.0)));
    
    // Drag force
    double drag_z = -F * v * vz;
    
    // Magnus force in z direction
    double mag_z = -p->B * p->omega * vx * sin(p->phi);
    
    return drag_z + mag_z - p->g;
  };


  // Stopping condition: stop when x >= 60 feet (18.288 m)
  auto f_stop = [](double t, const vector<double> &y, void *params) -> double {
    (void) t; (void) params;
    if (y[0] >= 18.288) return 1;  // 60 feet in meters
    return 0;
  };
  
  
  vector<pfunc_t> v_fun = {f_x, f_vx, f_y, f_vy, f_z, f_vz};
  
  pfunc_t stop_func = +f_stop;

  // Now we have everything to solve the trajectory
  vector<double> y = y0;
  double t = 0.0;
  double tmax = 4.0;  // seconds
  int nsteps = 2000;
  
  auto tgN = RK4SolveN(v_fun, y, nsteps, t, tmax, p_par, stop_func);
  
  // Extract final values and convert to feet
  xend = y[0] * 3.28084;
  yend = y[2] * 3.28084;
  zend = y[4] * 3.28084;
  vxend = y[1] * 3.28084;
  vyend = y[3] * 3.28084;
  vzend = y[5] * 3.28084;

  // ---------------------------------------------------------------------------


  // to compare to the plots in Fitzpatrick, output your results in **feet**
  // do not change these lines
  printf("********************************\n");
  printf("Coordinates when x=60 feet\n");
  printf("(x,y,x) = (%lf,%lf,%lf)\n",xend,yend,zend);
  printf("(vx,vy,vz) = (%lf,%lf,%lf)\n",vxend,vyend,vzend);
  printf("********************************\n");

  // plot the trajectory.  See Fitzpatrick for plot details
  if (showPlot){
    cout << "Press ^c to exit" << endl;
    theApp.SetIdleTimer(30,".q");  // set up a failsafe timer to end the program  
    theApp.Run();
  }
  
  return 0;
}

