#include "TF1.h"
#include "TMath.h"

void random4()
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
        gRandom3 = new TRandom3(2);
        gRandom4 = new TRandom3(3);

        double pi  = 3.14159265358;
        int n=1000000; int n1=0;
        double ld=7.;
        double p1, p2, p3;
        p1=3./ld;
        p2=3./ld;
        p3=1./ld;
        double x;

        TH1F *h = new TH1F("h","demo",100,0,3.);

        for (int i = 0; i < n; ++i) {
           double k1=gRandom->Uniform(0,1);
           double k2=gRandom2->Uniform(0,1);
           double k3=gRandom3->Uniform(0,1);
           double k4=gRandom4->Uniform(0,1);
           x=k2;
           if(k1>p1 && k1<(p1+p2)) {
            if(k3>x){x=k3;}
           } 
           if(k1>(p1+p2)) {
            if(k3>x){x=k3;}
            if(k4>x){x=k4;}
           } 
           double y=(2.-1.)*x+1.;
           h->Fill(y); 
        }
        h->SetLineColor(kRed);
        h->Scale(1.0/0.03/double(n));
        h->Draw();
        cout<<3./7.<<" "<<3.*4./7.<<endl;
        c1->SaveAs("sample_4_x.png");
}
