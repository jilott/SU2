#include<math.h>
#include<iostream>
#include<iomanip>

using namespace std;

void printJ(double J[][2]);
void printx(double *x);

int main() {

  int k;
  int kmax = 1000;
  double fk[2], dx[2], J[2][2];
  double re_theta,  f_lambda,  lambda, dre_dlamb;
  double re_theta0, f_lambda0;
  double newton_tol = 1e-10;
  double normres, factor;

  /*-- Test values -- */
  double tu    = 0.03;
  double rho   = 1.03838;
  double mu    = 1.72261e-05;
  double U     = 66.3579;
  double du_ds = 192.376;

  bool debug=false;

  /*-- Initial guess -- */
  f_lambda = 1.; 
  if (tu <= 1.3)
    re_theta = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
  else
    re_theta = (331.5*pow(tu-0.5658,-0.671))*f_lambda;
   
  /*-- Newton-Raphson loop --*/
  /*-- x={re_theta, f_lambda} --*/
  k = 0;
  while (1) {
    k++;
    if (k>=kmax)
      break;

    /*-- Evaluate residual for present guess --*/
    if (tu <= 1.3)
      re_theta0 = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
    else
      re_theta0 = (331.5*pow(tu-0.5658,-0.671))*f_lambda;

    re_theta0 = max(re_theta0,20.);

    lambda = pow(re_theta,2) * mu/(rho*pow(U,2))*du_ds;
    lambda = min(max(lambda,-0.1),0.1);

    if (lambda <= 0.0)
      f_lambda0 = 1 - (-12.986*lambda - 123.66*pow(lambda,2) - 405.689*pow(lambda,3))*exp(-pow(tu/1.5,1.5));
    else
      f_lambda0 = 1 + 0.275*(1-exp(-35.0*lambda))*exp(-tu/0.5);

    fk[0] = re_theta - re_theta0;
    fk[1] = f_lambda - f_lambda0;

    /*-- Flip sign of fk for RHS of Newton's method --*/
    fk[0] *= -1.; fk[1] *= -1.;

    /*-- Check for convergence --*/
    normres = sqrt(fk[0]*fk[0] + fk[1]*fk[1]);
    cout << normres << " " << re_theta << " " << f_lambda << " " << lambda << endl;
    if (normres < newton_tol) {
      cout << "Solution converged." << endl;
      cout << "x: {" << re_theta << " " << f_lambda <<  "}" << endl;
      break;
    }

    if (debug) {
    // cout << "Initial re_theta: " << re_theta << endl;
     cout << endl << "Initial fk: " << endl;
     printx(fk);
    }

    /*-- Evaluate the Jacobian of the nonlinear correlation system d(fk)/dx --*/
    J[0][0] = 1.;
    if (tu <= 1.3)
      J[0][1] = (1173.51 - 589.428*tu + 0.2196/pow(tu,2));
    else
      J[0][1] = (331.5*pow(tu-0.5658,-0.671));
    J[0][1] *= -1;

    dre_dlamb = 2*re_theta * mu/(rho*pow(U,2))*du_ds;
    if (lambda <= 0) 
      J[1][0] = 0. - (-12.986 - 2*123.66*lambda - 3*405.689*pow(lambda,2))*dre_dlamb*exp(-pow(tu/1.5,1.5));
    else
      J[1][0] = 0. + 0.275*(0.0+35.0*exp(-35.0*lambda))*dre_dlamb*exp(-tu/0.5);
    J[1][0] *= -1;
    J[1][1] = 1.;
    

    if (debug) {
      cout << "-----------------------------" << endl;
      cout << "Initial matrix: " << endl;
      printJ(J);
    }

    /*-- Forward solve --*/
    factor = J[1][0]/J[0][0];
    J[1][0] = 0.0;
    J[1][1] -= J[0][1]*factor;
    fk[1]   -= fk[0]*factor;

    if (debug) {
    //cout << endl << "Zero first column" << endl;
    //printJ(J);
    }

    if (debug) {
    // cout << endl << "Zero second column" << endl;
    // printJ(J);

    // cout << endl << "fk after forward solve " << endl;
    // printx(fk);
    }

    /*-- Back substitution --*/
    dx[1] = fk[1]/J[1][1];
    dx[0] = (fk[0]-J[0][1]*dx[1])/J[0][0];

    if (debug) {
     cout << endl << "dx after update: " << endl;
     printx(dx);

     cout << "-----------------------------" << endl;
    }

    re_theta += dx[0];
    f_lambda += dx[1];

    //break;
  }

}

void printx(double *x) {
  int k;

  for(k=0; k<2; k++)
    cout << x[k] << endl;
}

void printJ(double J[][2]) {

  int j,k;

  cout << scientific;
  for(j=0; j<2; j++) {
    for(k=0; k<2; k++) {
      cout << setw(10) << J[j][k] << "  ";
    }
    cout << endl;
  }
}
