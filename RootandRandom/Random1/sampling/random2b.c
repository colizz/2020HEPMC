#include "TF1.h"
#include "TMath.h"

void random2b()
{
// f(x,y)=1/(2*pi)*exp[-(x^2+y^2)/2]
// note \int_{-inf}^{inf} exp(-x^2/2)dx=sqrt(2*pi)
//Marsaglia Method

        TCanvas *c1= new TCanvas("c1", "demo", 10,10,700,700);

        gRandom = new TRandom3(0);
        gRandom2 = new TRandom3(1);

        double pi  = 3.14159265358;
        int n=8000000; int n1=0;
        double u,v,x,y,w,z;

        TH2F *h = new TH2F("h","demo",20,-1,1.,20,-1,1.);

        for (int i = 0; i < n; ++i) {
           u=gRandom->Uniform(0,1);
           v=gRandom2->Uniform(0,1);
           w=(2.*u-1.)*(2.*u-1.)+(2.*v-1.)*(2.*v-1.);   
           if(w>1) continue;
           z=sqrt(-2.*TMath::Log(w)/w);    
           x=(2.*u-1.)*z;
           y=(2.*v-1.)*z;
           h->Fill(x, y);
           n1++;   
        }

        h->SetLineColor(kRed);
        h->Scale(1.0/(0.1*0.1)/double(n));
        h->Draw("colz");


//        TH1D *hx = h->ProjectionX();
//          hx->Scale(1.0/(0.1)/double(n1));
//        hx->Draw();

//        TH1D *hy = h->ProjectionY();
//        hy->Scale(1.0/(0.1)/double(n1));
//        hy->Draw();

        cout<<double(n1)/double(n)<<" "<<pi/4.0<<" "<<1.0/sqrt(2.0*pi)<<endl;

        c1->SaveAs("sample2b_x.png");
}
