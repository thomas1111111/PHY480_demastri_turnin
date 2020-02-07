//filename: sum_up_down.cpp
//
//author: Thomas DeMastri email: demastri@msu.edu
//
//Revision History:
//
//Created 2/6/2020, 4:56 P.M.
//
//Changes:
//
// 2/6/2020:
//--created an s_up and s_down function, comparing the methods specified in 
//  Landau-Paez problem 3.4 in float precision (~1e-8)
//--wrote the main function which calculates s_up and s_down values for several 
//  values of N increasing exponentially (powers of 10), and stores log(error) in 
//  a .dat file
//
// 2/7/2020;
//
//	--had to store N in a long in order to fit 1e10
//	--had to change n from a float to an integer so that things would actually work.
//
//	Answers to part c):
//	There are two main regions of this graph, a linear portion between 100 and 1e7,
//	and a flat portion at 1e8. The error follows a power law with exponent 1.2 (extracted from
//	my linear fit), but after 1e8 the values of sum_up and sum_down shouldn't change. This is because
//	the machine precision of floats is 1e-8, so 1/1e8 will not change anything in the sum.
//
//	The upward sum starts at a larger value and thus has to truncate digits off in 1e-7 and 1e-8 earlier.
//	The downward sum starts at these smaller values and is thus more precise. 
//
//***************************************************************//

#include <iomanip>
#include <fstream>
using namespace std;		//all of the packages I think I'll need
#include <iostream>
#include <cmath>


float
s_up (long N)
{
  //For s_up, we start at 1 and work up to N
  float sum = 0.0;		//Our starting value for the sum is 0.0 (all floats by choice of the problem)
  for (long n=1;n<=N;n++)
    {
	sum += float(1)/ n;
    }
  cout << N << endl;
  return (sum);

}


float
s_down (long N)
{
  //For s_down, we start at N and work down to 1
  float sum = 0.0;		//Our starting value for the sum is 0.0 (all floats by choice of the problem)
  for (int n= N; n>0; n--)
    {
      sum += float(1)/ n;
    }

  cout << N << endl;

  return (sum);

}


int
main ()
{
  //This program will calculate sum_up and sum_down for various values of N
  //I've chosen to start N at a 'low' value of 10 and scale N by powers of 2
  //so that things are evenly spaced on the log-log plot
  
  //let's declare some things we're going to use later:
  float sum_up;
  float sum_down;
  float diff;
  
  
  //Let's open up a file to run to
  ofstream sumout ("sum_up_down.dat");

  sumout << "          N       log(N)     log(diff)" << endl;

  for(long N = 100; N <= 1e10; N*=10)
  {  
      sum_up = s_up(N);	//calculate our sum up and down values
      sum_down = s_down (N);

      diff = fabs (sum_up - sum_down) / (0.5 * (fabs (sum_up) + fabs (sum_down)));	//the difference divided by the average, as the problem specifies


      //write out the N, log(N), and the difference to our data file
      sumout << fixed << setprecision (8) << setw (11) << N << " "
	<< scientific << setprecision (8)
	<< setw (11) << log10 (N) << " " << setw (11) << log10 (diff) << endl;
    }
  
  cout << "done";

  return (0);
}
