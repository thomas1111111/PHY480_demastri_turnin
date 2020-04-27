//  Name: Thomas DeMastri     email: demastri@msu.edu
//
//  filename: integ.cpp
//
//  This puprose of this program is to test the accuracy
//  of the integration routines we wrote.
// 
//  This program is based off of integ_test.cpp from session_03
//
//  Revision History:
//  Created: 2/27/20 - Added integrand and function of integral
//                     Chose to use cos(x) as integrand since
//                     sin(x) is the analytical definite integral
//                     and there are routines to calculate this with
//                     high precision.
//  
//  Version 2:
//
//  Got my Simpson's Rule and Milne's rule to work. Adding the GSL Routine now.
//  Updated my integrand to A*cos(k*x+b) since that is a more interesting problem.
//
//  Version 3:
//
//  I had a lot of trouble getting the gsl integration routine to work,
//  since it requires a certain format for the integrand function to be in
//  that didn't play well with what I had already written.
//
//  Version 4:
//  I've decided it's easier to use the adaptive gsl routine.
//  I've adapted the other integration rules to work with the integrand
//  that gsl is expecting. 
//  Things seem to working now!
//
//
//DISCUSSION FOR 1 b.) and 1 c.):
//	There are three distinct sections on the error plot: For low values of N, both the simpson's and the milne's rule don't do a very good job
//	of estimating the value of the integral.
//	There is a middle region where the relative error in the integral begins to decrease like 1/N^a for both rules. As we'd expect of the Simpson's Rule,
//	it decreases like 1/N^4, and the Milne's rule picks up two extra orders of magntitude, decreasing like 1/N^6. 
//	Then, there is a roundoff region, where further increases in N do not increase precision.
//	It appears that the optimal amount of steps for the Milne's Rule occurs around N=10^5.
//	
//	The optimal step size is predicted analytically to occur at 10^-2, since the error ought to fall off at 1/N^6 and we're trying to satisfy
//	the equation 1/(N^6) = sqrt(N)*(10^-16) --> solved by N^13/2 = 10^16, N is about 200 or so.
//	I think the reason for this large discrpeancy is that the length of my integral is large; the error only starts to reasonably fall after
//	N is on the order of 10^3. After this, The error does hit the bottom two orders of magnitude from this tipping point. 
//
//	I think an adaptive routine would be able to take this into account. 
//
//DISCUSSION FOR 1 d.):
//	The method works because we cheat a little bit. We can analytically evaluate the integral 1/sqrt(x) very easily, 
//	and we can can handle the zero case of (1/(1+x) -1)/sqrt(x) much more easily than the zero case of (1/((1+x)*sqrt(x)))
//	In the former, we know that at x=0, we have a problem of 0/0, but in the limit the numerator dies off much more quickly so
//	the integrand approaches zero. It is much easier to correct for this than finding someway to deal with the gangly infinity occuring in 
//	the latter integrand.
//	I'm more impressed that the GSL routine was able to get an answer even without the tricks. The adaptive routine is no joke.
//
//We need to include some packages:
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
using namespace std;

//We'll include the routines we wrote as well:
#include "integ_methods.h"

//And declare our functions:
double test_integrand (double x, void *p);
double analytic_integrand (double x, void *p);
double tricky_integrand (double x, void *p);
double analytic_trick (double x, void *p);
double tricky_integrand_changed (double x, void *p);



int
main ()
{
  //We'll start by putting some bounds on our problem:
  double xmin = 1.;
  double xmax = 1000.;

  //And set up the parameters for our integrand:
  double A = 1.0;		//decay constant
  void *params = &A;

  //Open up a data file to write to:
  ofstream my_out;
  my_out.open ("integ.dat");
  my_out << setw (16) << "log(N)" << setw (16) << "log(Simp Rel Err)" <<
    setw (16) << "log(Miln Rel Err)" << setw (16) << "log(GSL Rel Err)" <<
    endl;

  double simps;
  double milne;
  double exact_answer = analytic_integrand (xmax, params) - analytic_integrand (xmin, params);
  double gsl_cop = integrate_gsl (*test_integrand, xmin, xmax, 1e-12, 1e-12, A);

  for (long N = 10; N < 1.1e7; N = N * 2)
    {

      simps = simpsons_rule (*test_integrand, N, xmin, xmax, A);
      milne = milnes_rule (*test_integrand, N, xmin, xmax, A);
  
      my_out << setw (16) << log10(N) << setw (16) <<
	log10(fabs ((simps - exact_answer) / exact_answer)) << setw (16) <<
	log10(fabs ((milne - exact_answer) / exact_answer)) << setw (16) <<
	log10(fabs ((gsl_cop - exact_answer) / exact_answer)) << endl;

	cout << N << endl;

    }

  my_out.close ();

  // This part estimates the integrand: 1/((1+x)*sqrt(x) from 0 to 2
  xmin = 0.;
  xmax = 2.;
  
  double B = 1.;

  long new_N = 40000.;
  simps = simpsons_rule(*tricky_integrand, new_N, xmin, xmax, B);
  milne = milnes_rule(*tricky_integrand, new_N, xmin, xmax, B); 
  gsl_cop = integrate_gsl (*tricky_integrand, xmin, xmax, 1e-12, 1e-12, B);
  exact_answer = analytic_trick(xmax,&B) - analytic_trick(xmin,&B);
   
  cout << endl;
  cout << "Estimations of 1/((1+x)*sqrt(x)) from 0 to 2 without any tricks:" << endl;
  cout << "Simpson's Rule: "<< simps  << endl;
  cout << "Milne's Rule: " << milne << endl;
  cout << "GSL Rule: " << gsl_cop << endl;;
  cout << "Exact Answer: " << exact_answer << endl;
  cout << endl;

  // I tried the trick outlined in by equation 3 in the handout.
  // The 2*sqrt(2)'s hanging around on the end of these is the value of the analytical integral 1/sqrt(x) evaluated
  // from 0 to 2. 
  
  
  simps = simpsons_rule(*tricky_integrand_changed, new_N, xmin, xmax, B)+sqrt(2.)*2;
  milne = milnes_rule(*tricky_integrand_changed, new_N, xmin, xmax, B)+2.*sqrt(2);
  gsl_cop = integrate_gsl (*tricky_integrand_changed, xmin, xmax, 1e-12, 1e-12, B)+sqrt(2.)*2.;
  exact_answer = analytic_trick(xmax,&B) - analytic_trick(xmin,&B);

  cout << endl;
  cout << "Estimations of 1/((1+x)*sqrt(x)) from 0 to 2 with tricks:" << endl;
  cout << "Simpson's Rule: "<< simps  << endl;
  cout << "Milne's Rule: " << milne << endl;
  cout << "GSL Rule: " << gsl_cop << endl;;
  cout << "Exact Answer: " << exact_answer << endl;
  cout << endl;

}



double
test_integrand (double x, void *p)
{
  //declaring our function to integrate over   

  double A = *(double *) p;
  return ((exp ((-A) * x)) * x);


}

double
tricky_integrand (double x, void *p)
{
  double B = *(double *) p;
  return ((1/(B+x))/sqrt(x));
}

double
tricky_integrand_changed (double x, void *p)
{
  double B = *(double *) p;
  if (isnan(((1/(B+x) - 1/B)/sqrt(x))) == true)
  {
	  return 0;
  }
  else
  {
	  return (((1/(B+x) - 1/B)/sqrt(x)));
  }
}


double analytic_trick (double x, void *p)
{
  double B = *(double *) p;
  return (2*atan(sqrt(x/B))/sqrt(B));
}

double
analytic_integrand (double x, void *p)
{
  //this is the exact analytical answer to our test integrand,
  //included here for error analysis

  double A = *(double *) p;

  return (-(exp ((-A) * x)) * (x + 1));
}
