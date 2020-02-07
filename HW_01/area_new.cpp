//  file: area_new.cpp
//
//  Name: Thomas DeMastri; email: demastri@msu.edu
//  
//  compile with g++ -o area_new.x area_new.cpp
//
//  Revision history
//
//  2/6/20:
//      - Copied over from the file area.cpp in the session_01 package
//      - let the user specify the precision with a cin statement (item one on to-do list)
//      - imported cmath, to do the second item on to-do list. Deleted the constant value of pi defined at the start.
//      - implemented an inline square function with the help of activity 5 materials, used it in my area calculation
//      - wrote a function for doing the area calculation.
//
//      - For item number 5, i think having a user input is much easier that editing a file for reading in the radius and precision
//      - It's now writing to a file area_new.dat with the radius, precision, and area labeled/provided.
//      
//      - For item 6, it now checks that radius and precision are positive numbers and makes sure that precision is an integer. 
//  To do:
//   1. output the answer with higher precision (more digits)
//   2. use a "predefined" value of pi or generate it with atan
//   3. define an inline square function
//   4. split the calculation off into a function (subroutine)
//   5. output to a file (and/or input from a file)
//   6. add checks of the input (e.g., for non-positive radii)
//   7. rewrite using a Circle class
//
//*********************************************************************// 





// include files
#include <iostream>		// this has the cout, cin definitions
#include <iomanip>		// included for setprecision
using namespace std;		// if omitted, then need std::cout, std::cin 
#include <fstream>
#include <cmath>

//*********************************************************************//

inline double			//function which takes in a double and squares it
sqr (double r)
{
  return r * r;
};

double
area_calc (double r)
{
// splitting the calulation into a subroutine:

  return M_PI * sqr (r);

}

int
main ()
{
  double radius;		// every variable is declared as int or double or...
  double precision;		// delcare precision as an int

  cout << "Enter the radius of a circle: " << endl;	// ask for radius
  cin >> radius;
  cout << "Enter the precision you want your answer" << endl;
  cin >> precision;

  double tol = 1e-14;

  if ((radius > 0) && (precision > 0) && (fabs (precision - int (precision)) < tol))	//check that we have numbers that are valid
    {
      ofstream aout ("area_new.dat");

      double area = area_calc (radius);	// standard area formula

      aout << setprecision (precision) << "radius = " << radius <<
	",  area = " << area << ", precision = " << precision << endl;

      cout << "Results outputted to area_new.dat" << endl;

    }
  else
    {
      cout <<
	"invlaid input, please enter radius and precision larger than 0, and an integer precision"
	<< endl;
    }
  return 0;			// "0" for successful completion
}

/*********************************************************************/
