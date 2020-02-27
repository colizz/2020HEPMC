#include "TF1.h"
#include "TMath.h"

void random3()
{
// xdx=1/2d(x^2)

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,700);
        gRandom = new TRandom3(0);
        gRandom2 = new TRandom3(1);

        int n=200000; 

        TH1F *h = new TH1F("h","demo",100,0,1.);
        TH1F *h0 = new TH1F("h0","demo",100,0,1.);

        for (int i = 0; i < n; ++i) {
           double x=gRandom->Uniform(0,1);
           double y=gRandom2->Uniform(0,1);
           if(y<=x) {h->Fill(x);}
           h0->Fill(x, x);
        }

        h->SetLineColor(kRed);
// dN/dX -> density
        h->Scale(1.0/0.01/double(n));
        h0->Scale(1.0/0.01/double(n));
        h->Draw();
        h0->Draw("same");

        c1->SaveAs("sample_3_x.png");
}
