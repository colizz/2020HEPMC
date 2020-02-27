#include "TF1.h"
#include "TMath.h"

void random0()
{
// xdx=1/2d(x^2)

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,700);
        gRandom = new TRandom3(0);
        int n=2000000; 

        TH1F *h = new TH1F("h","demo",100,0,1.);
        TH1F *h0 = new TH1F("h0","demo",100,0,1.);

        for (int i = 0; i < n; ++i) {
              double y=gRandom->Uniform(0,1);
              double x=sqrt(2.*y);  
           h->Fill(x);
           h0->Fill(y, y);
        }

        h->SetLineColor(kRed);
// dN/dX -> density
        h->Scale(1.0/0.01/double(n));
        h0->Scale(1.0/0.01/double(n));

        h->Draw();
        h0->Draw("same");

//        TF1 *f1=new TF1("f1","x",0,1);
//        f1->Draw("same");


        c1->SaveAs("sample_x.png");
}
