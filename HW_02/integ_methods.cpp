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
//		      - After a long time struggling to get this routine to work, I've given up
//		      - I'm attempting the problem with an adaptive routine.

#include <cmath>
#include "integ_methods.h"
#include <gsl/gsl_integration.h>

double 
simpsons_rule (double (*integrand)(double x, void *p), long N, double min, double max, double A)
{
   
	//This method integrates a function using simpson's rule	
   
   double x;
   void *params = &A;

   double h = ((max - min)/double(N - 1));  //calculate our mesh size
   double integral = 0; //initiate integral at 0

   for (long n=2; n<N; n+=2)                // Simpson's rule requires we loop over odds and evens,
   {		   // this loop does the odds
     x = min + h * double(n-1);
     integral += (4./3.)* h * integrand(x,params);
   }
   for (long n=3; n<N; n+=2)                // This loop goes over the evens
   {
     x = min + h * double(n-1);
     integral += (2./3.)* h * integrand(x,params);
   }
   // We need to remember to add in our endpoints:
   integral +=  (h/3.) * (integrand(max,params) + integrand(min,params));

   return (integral); // return the final calculated value
}


double 
milnes_rule (double (*integrand)(double x, void *p),long N, double min, double max, double A)
{
  
       	//This method integrates a function using milne's rule

   double x;
   void *params =&A;

   double h = ((max - min)/double(N - 1));  //calculate our mesh size
   double integral = 0; //initiate integral at 0

   for (long n=2; n<N; n+=4)                // Milne's rule requires seperate calculations
   {                                             // this loop does grabs the contributions from
    
     x = min + h * double(n-1);           // the elementary weights for 64/45
     integral += (64./45.)* h * integrand(x,params);

     x = min + h * double(n+1);
     integral += (64./45.) * h * integrand(x,params);
   }
   for (long n=3; n<N; n+=4)                // This loop goes over the contributions from the
   {
   
     x = min + h * double(n-1);     // other two elementary weights
     integral += (24/45.)* h * integrand(x,params);

     x = min + h * double(n+1);
     integral += (28/45.)* h * integrand(x,params); // we have to do twice the first elementary weight here since the intervals are meshed together
   }
   // We need to remember to add in our endpoints:
   integral +=  (14./45.)* h * (integrand(min,params) + integrand(max,params));

   return (integral); // return the final calculated value
}

double 
integrate_gsl(double (integrand) (double x, void *p), double min, double max, double eps_abs, double eps_rel, double A)
{

	//This method uses a built in GSL routine to integrate our given integrand:

    void *params = &A;
    double result;
    double error;

    gsl_integration_workspace *work_ptr = gsl_integration_workspace_alloc (2000);
	
    gsl_function My_func;

    My_func.function = integrand;
    My_func.params = params;

    gsl_integration_qags (&My_func, min, max, eps_abs, eps_rel, 1000, work_ptr, &result, &error);

    gsl_integration_workspace_free (work_ptr);

    return result;
}
