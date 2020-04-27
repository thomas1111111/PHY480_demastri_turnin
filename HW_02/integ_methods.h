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


extern double simpsons_rule (double (*integrand)(double x, void *p), long N, double min, double max, double A);

extern double milnes_rule (double (*integrand)(double x, void *p), long N, double min, double max,double A);

extern double integrate_gsl(double (*integrand) (double x, void *p), double min, double max, double eps_abs, double eps_rel, double A);

