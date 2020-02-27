#include <iostream>
#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "TTree.h"
#include "TVirtualPad.h"
#include "TApplication.h"
#include "TFile.h"

using namespace Pythia8;

int main(int argc, char* argv[]) {

	if (argc != 3) {
		cerr << " Unexpected number of command-line arguments. \n You are"
		<< " expected to provide one input and one output file name. \n"
		<< " Program stopped! " << endl;
		return 1;
	}

	// Create file
	TFile* outFile = new TFile(argv[2], "RECREATE");

	// Book Tree.
	int scale;
	int event, pid, status;
	double pt, m, px, py, pz, gm, theta, phi, e;

	TTree *h = new TTree("h","h");
	h->Branch("event",&event);
	h->Branch("pid",&pid);
	h->Branch("status",&status);
	h->Branch("pt",&pt);
	h->Branch("m",&m);
 	h->Branch("px",&px);
	h->Branch("py",&py);
 	h->Branch("pz",&pz);
 	h->Branch("gm",&gm);
	h->Branch("theta",&theta);
 	h->Branch("phi",&phi);
	h->Branch("e",&e);
	h->Branch("scale",&scale);	

  	// specify an input file
	HepMC::IO_GenEvent ascii_in(argv[1],std::ios::in);

  	// get the first event
  	HepMC::GenEvent* evt = ascii_in.read_next_event();
	
	int scale_h, scale_w, scale_z, scale_gamma, middle, is_over, is_write, copy_time=0;
  	int pN=0, vN, EveN=0, vbarcode=-1;
  	int w=0, z=0, gamma=0, q=0, tau=0, mu=0, elec=0, gluon=0;
  	// loop until we run out of events

        HepMC::GenEvent::vertex_iterator tempt_v;
        HepMC::GenVertex::particles_out_const_iterator tempt_p;
        while ( evt ) {
  		EveN++;
		is_over = 0;
		is_write =0;
		//if(EveN%1000==0) cout<<"event_number= "<<EveN<<endl;
		//vertex barcode 
		scale_h = 0;
		scale_w = 0;
		scale_z = 0;
		scale_gamma = 0;
		middle = 0;
		copy_time = 0;
		cout<<"event= -----------------"<<EveN<<endl;
  		for ( HepMC::GenEvent::vertex_iterator v = evt->vertices_begin();v != evt->vertices_end(); ++v ){
			if(is_over==1 && scale_z ==0 && scale_w ==0) continue; //h > l gamma quark gluon
			if(is_over==1 && (scale_z ==2 || scale_w ==2)) continue; 
			if(is_over==2 && (scale_z == 1 || scale_w ==1)){
				if(middle==0){
					is_write=0;
					v =tempt_v;
					vbarcode = (*v)->barcode();
				}
				middle++;
			}		
	
			if(scale_h == 1){
				if(vbarcode!=(*v)->barcode()) continue;
			}
                	for ( HepMC::GenVertex::particles_out_const_iterator p = (*v)->particles_out_const_begin(); p!=(*v)->particles_out_const_end(); ++p ){
				if(middle==1) {middle++; is_over = 0; cout<<"middle"<<middle<<endl;p = tempt_p;}
				if(abs((*p)->pdg_id())>50) continue;
				if((*p)->pdg_id()==25){
					scale_h = 1;
					//h income vertex barcode
					vbarcode	=	(*p)->end_vertex()->barcode();
        				cout<<"h--- "<<vbarcode<<endl;
					if((*p)->status()==22)	{
						scale	=	0;
						event	=	EveN;
                   		     		status  =       (*p)->status();
						pid     =       (*p)->pdg_id();
						pt      =       (*p)->momentum().perp();
                        			px      =       (*p)->momentum().px();
                       		 		py      =       (*p)->momentum().py();
                       				pz      =       (*p)->momentum().pz();
                       				m       =       (*p)->momentum().m();
                       				gm      =       (*p)->generated_mass();
                       				theta   =       (*p)->polarization().theta();
						phi     =       (*p)->polarization().phi();
						e	=	(*p)->momentum().e();
    						h->Fill();               
                		
					}	
				break;
				}
				
				//for z boson
				else if(abs((*p)->pdg_id())==23 && scale_h ==1){	     
//					scale_z++;
	                        	vbarcode	=	(*p)->end_vertex()->barcode();
					cout<<"z "<<vbarcode<<endl;
					copy_time++;
//					if(is_over == 4) continue;
					if(copy_time==1){
						if(scale_gamma==1) scale_z++;
//						scale_z++;
						else {
							tempt_v = v;
							tempt_p = p+1;
							cout<<"2 z vbarcode "<<(*p)->end_vertex()->barcode()<<" "<<(*tempt_p)->end_vertex()->barcode()<<endl;
						}
					}
cout<<"is_write "<<is_write<<endl;
				
                        		if(scale_h==1 && is_write ==0)  {
						scale	=	2502;
						z	=	z + 1;
                        	        	event   =       EveN;
                                		status  =       (*p)->status();
                                		pid     =       (*p)->pdg_id();
                               		 	pt      =       (*p)->momentum().perp();
                                		px      =       (*p)->momentum().px();
                                		py      =       (*p)->momentum().py();
                   		             	pz      =       (*p)->momentum().pz();
                                		m       =       (*p)->momentum().m();
                            		    	gm      =       (*p)->generated_mass();
                                		theta   =       (*p)->polarization().theta();
                       		         	phi     =       (*p)->polarization().phi();
                                		e       =       (*p)->momentum().e();
                                		h->Fill();
                                        cout<<"z "<<vbarcode<<endl;
						scale_z++;						
						is_write = 1;
						
					}

					break;
				}
				
				//for w boson
                                else if(abs((*p)->pdg_id())==24 && scale_h ==1){
//					scale_w++;
                                        vbarcode        =       (*p)->end_vertex()->barcode();
                                        cout<<"w "<<vbarcode<<endl;
					copy_time++;
//					cout<<copy_time<<endl;
//					if(is_over == 4) continue;
                                        if(copy_time==1){
                                                tempt_v = v;
//						scale_w++;
                                                tempt_p = p+1;
//						cout<<"2 vbarcode "<<(*p)->end_vertex()->barcode()<<" "<<(*tempt_p)->end_vertex()->barcode()<<endl;
                                        }
//cout<<"is_write "<<is_write<<endl;
                                        if(scale_h==1 && is_write==0)  {
						scale	=	2502;
                                                w       =       w + 1;
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();
//                                        cout<<"w "<<vbarcode<<endl;
						scale_w++;
						is_write = 1;

                                        }
					break;

                                }

				//for gamma
                                else if((*p)->pdg_id()==22 && scale_h ==1){
                                        //vbarcode        =       (*p)->end_vertex()->barcode();
//                                        cout<<"gamma "<<endl;
					scale_gamma++;
					if(scale_z==1) scale_z = 2;
                                        if(scale_h==1)  {
						scale	=	2502;
                                                gamma       =       gamma + 1;
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();
                                        cout<<"gamma "<<endl;

						is_over = 1;
					

                                        }

                                }

				//for quark
                                else if(abs((*p)->pdg_id())<=5 && scale_h ==1){
                                        vbarcode        =       (*p)->end_vertex()->barcode();
                                        cout<<"q "<<vbarcode<<endl;
                                        if(scale_h==1)  {
						scale = 2502;
						if(scale_w!=0) { scale = 250202402;}
						if(scale_z!=0) { scale = 250202302;}
                                                q       =       q + 1;
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();
                                        	cout<<"scale_z="<<scale_z<<" , scale_w="<<scale_w<<" , quark_pid= "<<pid<<endl;

						is_over = 1;
						if(scale_w == 1 || scale_z ==1) is_over = 2;
                                                if(scale_w == 2 || scale_z ==2) is_over = 1;


//cout<<"is_over "<<is_over<<" "<<scale_w<<" "<<scale_z<<endl;
                                        }

                                }

				//for charge lepton
                                else if((abs((*p)->pdg_id())==15||abs((*p)->pdg_id())==13||abs((*p)->pdg_id())==11) && scale_h ==1){
//                                        vbarcode        =       (*p)->end_vertex()->barcode();
//                                        cout<<"tau "<<vbarcode<<endl;
                                        if(scale_h==1)  {
                                                scale = 2502;
                                                if(scale_w!=0) {scale = 250202402;}
                                                if(scale_z!=0) {scale = 250202302;}
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();

                                                cout<<"scale_z="<<scale_z<<" , scale_w="<<scale_w<<" , charge_lepton_pid= "<<pid<<endl;

  
                                                is_over = 1;
                                                if(scale_w == 1 || scale_z ==1) is_over = 2;
                                                if(scale_w == 2 || scale_z ==2) is_over = 1;

                                        }

                                }


                                //for neutrios
                                else if((abs((*p)->pdg_id())==12||abs((*p)->pdg_id())==14||abs((*p)->pdg_id())==16) && scale_h ==1){
                                        //vbarcode        =       (*p)->end_vertex()->barcode();
//                                        cout<<"elec "<<endl;
                                        if(scale_h==1)  {
                                                scale = 2502;
                                                if(scale_w!=0) {scale = 250202402;}
                                                if(scale_z!=0) {scale = 250202302;}

                                                //elec       =       elec + 1;
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();

                                                cout<<"scale_z="<<scale_z<<" , scale_w="<<scale_w<<" , neutrios_pid= "<<pid<<endl;


                                                is_over = 1;
                                                if(scale_w == 1 || scale_z ==1) is_over = 2;
                                                if(scale_w == 2 || scale_z ==2) is_over = 1;


                                        }

                                }

 
                                //for gluon
                                else if(abs((*p)->pdg_id())==21 && scale_h ==1){
//                                        vbarcode        =       (*p)->end_vertex()->barcode();
//                                        cout<<"gluon "<<vbarcode<<endl;
                                        if(scale_h==1)  {
						scale 	= 2502;
                                                if(scale_w!=0) {scale = 250202402;}
                                                if(scale_z!=0) {scale = 250202302;}

                                                gluon       =       gluon + 1;
                                                event   =       EveN;
                                                status  =       (*p)->status();
                                                pid     =       (*p)->pdg_id();
                                                pt      =       (*p)->momentum().perp();
                                                px      =       (*p)->momentum().px();
                                                py      =       (*p)->momentum().py();
                                                pz      =       (*p)->momentum().pz();
                                                m       =       (*p)->momentum().m();
                                                gm      =       (*p)->generated_mass();
                                                theta   =       (*p)->polarization().theta();
                                                phi     =       (*p)->polarization().phi();
                                                e       =       (*p)->momentum().e();
                                                h->Fill();

                                                cout<<"scale_z="<<scale_z<<" , scale_w="<<scale_w<<" , gluon_pid= "<<pid<<endl;

	                                        cout<<"gluon "<<endl;

                                                is_over = 1;

                                                if(scale_w == 1 || scale_z ==1) is_over = 2;
                                                if(scale_w == 2 || scale_z ==2) is_over = 1;
                                        }

                                }



				else { continue;}

			}
	  

  		}

  	delete evt;
  	ascii_in >> evt;
	}

	cout<<"h:w:z:q:gamma:tau:gluon:mu:elec =  "<<EveN<<" : "<<w<<" : "<<z<<" : "<<q<<" : "<<gamma<<" : "<<tau<<" : "<<gluon<<" : "<<mu<<" : "<<elec<<endl;

  	h->Write();
  	//delete outFile;
  	return 0;
}

