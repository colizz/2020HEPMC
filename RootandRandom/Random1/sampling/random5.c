#include "TF1.h"
#include "TMath.h"

void random5()
{
//Maxwell Distri.
// f(x)=2*sqrt(x)/sqrt(pi)*exp(-x), 0<x<inf
// f(x)=H(x)*g(x)
//    g(x)=2/3*exp(-2x/3)
//    H(x)=3/sqrt(pi)*sqrt(x)*exp(-x/3)
//     L=max(H(x))=sqrt(27/2/pi/e)

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,500);
        gRandom = new TRandom3(0);
        gRandom2 = new TRandom3(1);

        double pi  = 3.14159265358;
        int n=1000000; int n1=0;
        double u,v;
        double L=sqrt(27./2./pi/exp(1.0));

        TH1F *h = new TH1F("h","demo",100,0,5.);

        for (int i = 0; i < n; ++i) {
           double k1=gRandom->Uniform(0,1);
           double k2=gRandom2->Uniform(0,1);
           double x=-3./2.*TMath::Log(k1);
           if(k2*k2<=(-exp(1.0)*k1*TMath::Log(k1))) {h->Fill(x); n1++;}
        }
        h->SetLineColor(kRed);
// dN/dX -> density
        h->Scale(1.0/0.05/double(n1));
        h->Draw();

        cout<<double(n1)/double(n)<<" "<<1.0/L<<" "<<2./sqrt(pi)/exp(1.0)<<" "<<2./sqrt(pi)*sqrt(2.)*exp(-2.0)<<endl;

        c1->SaveAs("sample_5_x.png");
}
