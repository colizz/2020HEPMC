//  Breit-Wigner
//   mk=10., gk=1.0;
//   int_0^200 sqrt(x)/((x-mk*mk)*(x-mk*mk)+mk*mk*gk*gk);
//   theta: atan(-10.) to atan(10.)
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>
using namespace std;

double func(double x);

int main()
{
   srand(unsigned(time(0)));  // seed the randomizer
   int i, n;
   double crude_mc, x, sum_sigma, fx, variance;
   cout << "Read in the number of Monte-Carlo samples" << endl;
   cin >> n;
   crude_mc = sum_sigma=0. ;
   // evaluate the integral with the a crude Monte-Carlo method
   for ( i = 1; i <= n; i++){
      x=rand()/(double)RAND_MAX;
      x=200.*x; // 0-200
      fx=func(x);
      crude_mc += 200.*fx;
      sum_sigma += 200.*200.*fx*fx;
   }
   crude_mc = crude_mc/((double) n );
   sum_sigma = sum_sigma/((double) n );
   variance=sum_sigma-crude_mc*crude_mc;
// above seems to be just for V[f]
   variance=variance/double(n);

   cout << " variance= " << sqrt(variance) << " Integral = "
        << crude_mc<< endl;
} // end of main program

// this function defines the function to integrate
double func(double x)
{
   double value;
   double mk=10., gk=1.0;
   value = sqrt(x)/((x-mk*mk)*(x-mk*mk)+mk*mk*gk*gk);
   return value;
} // end of function to evaluate

