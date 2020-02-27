#include "TF1.h"
#include "TMath.h"

void random2()
{
// f(x,y)=1/(2*pi)*exp[-(x^2+y^2)/2]
// note \int_{-inf}^{inf} exp(-x^2/2)dx=sqrt(2*pi)

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,700);

        gRandom = new TRandom3(0);
        gRandom2 = new TRandom3(1);

        double pi  = 3.14159265358;
        int n=2000000; 
        double u,v, x, y;

        TH2F *h = new TH2F("h","demo",20,-1,1.,20,-1,1.);

        for (int i = 0; i < n; ++i) {
           u=gRandom->Uniform(0,1);
           v=gRandom2->Uniform(0,1);
           x=sqrt(-2.0*TMath::Log(u))*TMath::Cos(2.*pi*v);
           y=sqrt(-2.0*TMath::Log(u))*TMath::Sin(2.*pi*v);
           h->Fill(x, y);
        }

        h->SetLineColor(kRed);
        h->Scale(1.0/(0.1*0.1)/double(n));
        h->Draw("colz");

/*
        TH1D *hx = h->ProjectionX();
          hx->Scale(1.0/(0.1)/double(n));
        hx->Draw();
*/
//        TH1D *hy = h->ProjectionY();
//        hy->Scale(1.0/(0.1)/double(n));
//        hy->Draw();

        cout<<1.0/sqrt(2.0*pi)<<endl;

        c1->SaveAs("sample2_x.png");
}
