//  Name: Thomas DeMastri      email: demastri@msu.edu
//
//  filename: integ_methods.h
//
//  This is a header file for integ_methods.cpp
//  It declares the two functions I need for the HW
//  This script is based off of integ_routines.h in session 3
//
//  Revision history:
//  Created 2/27/20 - Added simpsons_rule and milnes_rule declaration.
//


extern double simpsons_rule (int N, double min, double max, double (*integrand) (double x, double A, double k, double b) );

extern double milnes_rule (int N, double min, double max, double (*integrand) (double x, double A, double k, double b) );

extern double integrand_gsl (double x, void *p);

struct //created solely to satisty the GSL routine
integ_params
{
   double A; //Amplitude
   double k; //wave number
   double b; //phase
};

extern double integrate_gsl(double (integrand) (double x, void *p), integ_params, double min, double max, double eps_abs, double eps_rel);

