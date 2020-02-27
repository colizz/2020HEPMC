void random()
{
//https://root.cern.ch/doc/master/classTRandom1.html
//https://root.cern.ch/doc/master/TRandom1_8cxx_source.html
//The algorithm for this random engine has been taken from the 
//original implementation in FORTRAN by Fred James as part of 
//CLHEP.
//The initialisation is carried out using a Multiplicative 
//Congruential generator using formula constants of L'Ecuyer as
// described in "F.James, Comp. Phys. Comm. 60 (1990) 329-344". 

//https://root.cern.ch/doc/master/classTRandom2.html
//https://root.cern.ch/doc/master/TRandom2_8cxx_source.html
//Random number generator class based on the maximally 
//quidistributed combined Tausworthe generator by L'Ecuyer.
//The period of the generator is 2**88 (about 10**26) and it 
//uses only 3 words for the state.

//http://root.cern.ch/root/html/TRandom3.html
//https://root.cern.ch/doc/master/TRandom3_8cxx_source.html
//TRandom3, is based on the "Mersenne Twister generator", and is the recommended one, since it has good random proprieties (period of about 10**6000 ) and it is fast
        // create random number generator
 gRandom = new TRandom2(0);



 double mean=0.0; 
 double variance=0.0;
 double chi2=0.0;
 double f[10];
 for(int i=0; i<10; i++) {
  f[i]=0.;
 }
 int NN;
 cout << "Number =";
 cin >> NN;

 
 for(int i=1; i<NN; i++) {
  
   double fx=gRandom->Uniform(0,1);

   if(fx<0.1) {
     f[0]+=1.0; 
   }
   else if(fx<0.2) {
     f[1]+=1.0; 
   }
   else if(fx<0.3) {
     f[2]+=1.0; 
   }
   else if(fx<0.4) {
     f[3]+=1.0; 
   }
   else if(fx<0.5) {
     f[4]+=1.0; 
   }
   else if(fx<0.6) {
     f[5]+=1.0; 
   }
   else if(fx<0.7) {
     f[6]+=1.0; 
   }
   else if(fx<0.8) {
     f[7]+=1.0; 
   }
   else if(fx<0.9) {
     f[8]+=1.0; 
   }
   else{
     f[9]+=1.0; 
   }
   mean+=fx;
   variance+=fx*fx;
 }
   mean=mean/double(NN);
   variance=sqrt(variance/double(NN)-mean*mean);

 for(int i=0; i<10; i++) {
  cout<<i<<" "<<f[i]/double(NN)<<endl;
  chi2+=(f[i]-double(NN)/10.0)*(f[i]-double(NN)/10.0)/(double(NN)/10.0);
 }
  chi2=chi2/9.0;

 cout<<"mean= "<<mean<<endl;  
 cout<<"variance= "<<variance<<endl;  
 cout<<"chi2= "<<chi2<<endl;  



}
