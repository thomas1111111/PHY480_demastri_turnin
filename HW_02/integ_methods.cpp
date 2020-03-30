// Name: Thomas DeMastri      email: demastri@msu.edu
//
// filename: integ_methods.cpp
//  
// The purpose of this script is to create methods that integrate
// over a user specified integrand using Simpson's rule and Milne's rule.
// This program is based off of integ_routines.cpp from session 3.
//
// Revision History:
//   Created: 2/27/20 - Simpson's rule method added, updated to double percision
//			from integ_routines
//		      - Milne's rule method added, adapted from session three notes
//		      - After a quick search, I've elected to use a non-adaptive 
//		        integration routine from GSL. It doesn't seem 'fair' to compare
//		        an adaptive routine to the one's written here.

#include <cmath>
#include "integ_methods.h"
#include <gsl/gsl_integration.h>

double 
simpsons_rule (int N, double min, double max, double (*integrand) (double x, double A, double k, double b))
{
   
	//This method integrates a function using simpson's rule	
   
   double x, A, k, b;

   double h = ((max - min)/float(N - 1));  //calculate our mesh size
   double integral = 0; //initiate integral at 0

   for (int n=2; n<N; n+=2)                // Simpson's rule requires we loop over odds and evens,
   {		   // this loop does the odds
     x = min + h * double(n-1);
     integral += (4./3.)* h * integrand(x,A,k,b);
   }
   for (int n=3; n<N; n+=2)                // This loop goes over the evens
   {
     x = min + h * double(n-1);
     integral += (2./3.)* h * integrand(x,A,k,b);
   }
   // We need to remember to add in our endpoints:
   integral +=  (h/3.) * (integrand(max,A,k,b) + integrand(min,A,k,b));

   return (integral); // return the final calculated value
}


double 
milnes_rule (int N, double min, double max, double (*integrand) (double x, double A, double k, double b) )
{
  
       	//This method integrates a function using milne's rule

   double x, A, k, b;
	
   double h = ((max - min)/float(N - 1));  //calculate our mesh size
   double integral = 0; //initiate integral at 0

   for (int n=2; n<N; n+=4)                // Milne's rule requires seperate calculations
   {                                             // this loop does grabs the contributions from
    
     x = min + h * double(n-1);           // the elementary weights for 64/45
     integral += (64./45.)* h * integrand(x,A,k,b);

     x = min + h * double(n+1);
     integral += (64./45.) * h * integrand(x,A,k,b);
   }
   for (int n=3; n<N; n+=4)                // This loop goes over the contributions from the
   {
   
     x = min + h * double(n-1);     // other two elementary weights
     integral += (24/45.)* h * integrand(x,A,k,b);

     x = min + h * double(n+1);
     integral += (28/45.)* h * integrand(x,A,k,b); // we have to do twice the first elementary weight here since the intervals are meshed together
   }
   // We need to remember to add in our endpoints:
   integral +=  (14./45.)* h * (integrand(min,A,k,b) + integrand(x,A,k,b));

   return (integral); // return the final calculated value
}

double 
integrate_gsl(double (integrand) (double x, void *p), integ_params params, double min, double max, double eps_abs, double eps_rel)
{

	//This method uses a built in GSL routine to integrate our given integrand:

    gsl_function F;
    F.function = integrand;   
    F.params = reinterpret_cast<void *>(&params); //recasting our parameters as void parameters

    double result, error;
    size_t neval;

    int func = gsl_integration_qng(&F, min, max, eps_abs,eps_rel,&result,&error,&neval);

    if (func)
    {
	return(-1.);
    }
    else
    {
	return(result);
    }

    return (-2.);

}
