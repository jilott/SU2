#include<math.h>
#include<iostream>
#include<iomanip>

using namespace std;

/*-- Test values -- */
double tu    = 0.03;
double rho   = 1.03838;
double mu    = 1.72261e-05;
double U     = 66.3579;
double du_ds = 1.92376;

double corr_func(double lambda) {

  double re_thetaRHS, re_thetaLHS;
  double f_lambda, f;

  /*-- compute f --*/
  re_thetaLHS = sqrt(rho*pow(U,2)/(mu*du_ds)*lambda);

  if (lambda <= 0.0)
    f_lambda = 1 - (-12.986*lambda - 123.66*pow(lambda,2) - 405.689*pow(lambda,3))*exp(-pow(tu/1.5,1.5));
  else
    f_lambda = 1 + 0.275*(1-exp(-35.0*lambda))*exp(-tu/0.5);

  if (tu <= 1.3)
    re_thetaRHS = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
  else
    re_thetaRHS = (331.5*pow(tu-0.5658,-0.671))*f_lambda;

  //cout << "LHS, RHS: " << re_thetaLHS << ", " << re_thetaRHS << endl;
  f = re_thetaLHS-re_thetaRHS;

  return f;
}

int main() {

  double lambda_a, lambda_b, lambda_c, lambda;
  double f_a, f_b, f_c;

  double bracket_tol = 1e-10;

  double re_theta, f_lambda;
  double lambda_check;

  bool debug=false;

  /*-- Initial bracket for lambda --*/
  if (du_ds>=0) {
    lambda_a = 0.0;
    lambda_b = 0.1;
  } else {
    lambda_a = -0.1;
    lambda_b =  0.0;
  }

  f_a = corr_func(lambda_a);
  f_b = corr_func(lambda_b);
  if (f_a*f_b > 0) {
    cout << "ERROR: f_a and f_b have same sign!" << endl;
    lambda = 0.1*copysign(1.0,du_ds);
  } else {

    /*-- Begin method of false position to solve for lambda --*/
    while(1) {

      /*-- Mathews and Fink, p.57 eq. 22 --*/
      lambda_c = lambda_b - (f_b*(lambda_b-lambda_a))/(f_b-f_a);
      f_c = corr_func(lambda_c);

      /*-- Monitor convergence --*/
      cout << lambda_a << " " << lambda_c << " " << lambda_b << "\t" << f_a << " " << f_c << " " << f_b << endl;

      if (f_c==0.0) {
        lambda = lambda_c;
        cout << "Exact solution found: lambda=" << lambda << endl;
        break;
      }
      if (lambda_c-lambda_a <= bracket_tol) {
        lambda = 0.5*(lambda_a+lambda_c);
        cout << "Bracket tolerance reached. lambda = " << lambda << endl;
        break;
      } else if (lambda_b-lambda_c <= bracket_tol) {
        lambda = 0.5*(lambda_b+lambda_c);
        cout << "Bracket tolerance reached. lambda = " << lambda << endl;
        break;
      }

      /*-- Modify bracket for next iteration --*/
      if (f_a*f_c <=0) {
        lambda_b = lambda_c;
        f_b      = f_c;
      } else {
        lambda_a = lambda_c;
        f_a      = lambda_a;
      }


    }
  }


  if (lambda <= 0.0)
    f_lambda = 1 - (-12.986*lambda - 123.66*pow(lambda,2) - 405.689*pow(lambda,3))*exp(-pow(tu/1.5,1.5));
  else
    f_lambda = 1 + 0.275*(1-exp(-35.0*lambda))*exp(-tu/0.5);

  if (tu <= 1.3)
    re_theta = (1173.51 - 589.428*tu + 0.2196/pow(tu,2))*f_lambda;
  else
    re_theta = (331.5*pow(tu-0.5658,-0.671))*f_lambda;

  lambda_check = pow(re_theta,2)*mu/(rho*pow(U,2))*du_ds;


  cout <<  "Solution: " << endl;
  cout << lambda << " " << lambda_check << " " << f_lambda << " " << re_theta << endl;

}
