#include "TF1.h"
#include "TMath.h"

void random3b()
{
// f(x)=1/sqrt(2*pi)*exp(-x^2/2), 0<x<inf
// h(x)=exp(-x)

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,500);
        gRandom = new TRandom3(0);
        gRandom2 = new TRandom3(1);

        double pi  = 3.14159265358;
        int n=5000000; int n1=0;
        double u,v;
        double L=sqrt(exp(1.0)/2.0/pi);

        TH1F *h = new TH1F("h","demo",100,0,3.);

        for (int i = 0; i < n; ++i) {
           double k1=gRandom->Uniform(0,1);
           double k2=gRandom2->Uniform(0,1);
           double x=-TMath::Log(k1);
           v=L*exp(-(x-1.)*(x-1.)/2.0)/L;
           if(k2<=v) {h->Fill(x,L); n1++;}
        }

        h->SetLineColor(kRed);
// dN/dX -> density
        h->Scale(1.0/0.03/double(n));
        h->Draw();

        cout<<double(n1)/double(n)<<" "<<1.0/L/2.0<<" "<<1.0/sqrt(2.0*pi)<<" "<<1.0/sqrt(2.0*pi)*exp(-0.5)<<endl;

        c1->SaveAs("sample_3b_x.png");
}
