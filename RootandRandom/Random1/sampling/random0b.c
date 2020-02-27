#include "TF1.h"
#include "TMath.h"

void random0b()
{
// a*exp(-ax)dx=d(1-exp(-a*y))
// a=2.0 here

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,700);
        gRandom = new TRandom3(0);
        int n=2000000; 

        TH1F *h = new TH1F("h","demo",1000,0,1.);

        for (int i = 0; i < n; ++i) {
              double y=gRandom->Uniform(0,1);
              double x=-1./2.*TMath::Log(1.0-y);
          h->Fill(x);
        }

        h->SetLineColor(kRed);
// dN/dX -> density

        h->Scale(1.0/0.001/double(n));

        h->Draw();
        cout<<2.0*exp(-2.)<<endl;

        c1->SaveAs("sample_0b_x.png");
}
