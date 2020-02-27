//  int_0^1 exp(-x/2) dx =2(1-exp(-0.5))
//  importance sampling: g(x)=1-x/2.0
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
   double crude_mc, tx, x, sum_sigma, fx, variance;
   cout << "Read in the number of Monte-Carlo samples" << endl;
   cin >> n;
   crude_mc = sum_sigma=0. ;
   // evaluate the integral with the a crude Monte-Carlo method
   for ( i = 1; i <= n; i++){
      tx=3.0/4.0*rand()/(double)RAND_MAX;
      //0=<tx<=3/4;

      x=2.*(1-sqrt(1.0-tx));
      fx=3.0/4.0*func(x);
      crude_mc += fx;
      sum_sigma += fx*fx;
   }
   crude_mc = crude_mc/((double) n );
   sum_sigma = sum_sigma/((double) n );
   variance=sum_sigma-crude_mc*crude_mc;
// above seems to be just for V[f]
   variance=variance/double(n);

   cout << " Deviation= " << sqrt(variance) << " Integral = "
        << crude_mc<< endl;
   cout << " Exact= " <<2.*(1.-exp(-0.5))<< endl;
} // end of main program

// this function defines the function to integrate
double func(double x)
{
   double value;
   value = exp(-x/2.0)/(1.-x/2.0);
   return value;
} // end of function to evaluate

