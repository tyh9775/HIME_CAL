#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>

void h_cmb(TH1D *h, TFile *f,std::string hname, double scl,int num){

    TH1D *h_temp = (TH1D*)f->Get(hname.c_str())->Clone();
    double N_ent = h->GetEntries();
    double N_ent_temp = h_temp->GetEntries();
    h_temp->Scale(num); //undo prob normalization
    double width = h_temp->GetXaxis()->GetBinWidth(1);
    h_temp->Scale(width); //undo density normalization
    h_temp->Scale(scl); 
    h->Add(h_temp);
    h->SetEntries(N_ent_temp+N_ent);
    h_temp->Delete();
    return;
}

void ha_cmb(TH1D *h, TFile *f,std::string hname, double scl,int num){

    TH1D *h_temp = (TH1D*)f->Get(hname.c_str())->Clone();
    double N_ent = h->GetEntries();
    double N_ent_temp = h_temp->GetEntries();
    h_temp->Scale(num); //undo prob normalization
    double width = h_temp->GetXaxis()->GetBinWidth(1);
    h_temp->Scale(width*width); //undo area density normalization
    h_temp->Scale(scl); 
    h->Add(h_temp);
    h->SetEntries(N_ent_temp+N_ent);
    h_temp->Delete();
    return;
}

void h2_cmb(TH2D *h, TFile *f,std::string hname, double scl,int num){
    TH2D *h_temp = (TH2D*)f->Get(hname.c_str())->Clone();
    double N_ent = h->GetEntries();
    double N_ent_temp = h_temp->GetEntries();
    h_temp->Scale(num); //undo prob normalization
    double width_x = h_temp->GetXaxis()->GetBinWidth(1);
    double width_y = h_temp->GetYaxis()->GetBinWidth(1);
    h_temp->Scale(width_x*width_y); //undo density normalization
    h_temp->Scale(scl); 
    h->Add(h_temp);
    h->SetEntries(N_ent_temp+N_ent);
    h_temp->Delete();
    return;
}


void h_den_write(TH1D *h, TFile *f, int total_events){
    f->cd();
    double n = h->GetEntries();
    h->Scale(1.0, "width");
    h->Scale(1.0/total_events);
    //h->SetEntries(n);
    h->Write();
    return;
}

void ha_den_write(TH1D *h, TFile *f, int total_events){
    f->cd();
    double n = h->GetEntries();
    h->Scale(1.0, "width");
    h->Scale(1.0, "width");
    h->Scale(1.0/total_events);
    //h->SetEntries(n);
    h->Write();
    return;
}

void h2_den_write(TH2D *h, TFile *f, int total_events){
    f->cd();
    double n = h->GetEntries();
    h->Scale(1.0, "width");
    h->Scale(1.0/total_events);
    //h->SetEntries(n);
    h->Write();
    return;
}

void sb_combiner_v2(
    const std::vector<std::string> &inputs = {"input1.root,input2.root"},
    const std::string &ofile = "output.root",
    bool useVetoWall = true,
    const std::vector<int> &num_ent = {1,1},
    double scaler = 1.23
){
    double scale_factor = 1.0/scaler;
    TFile *fo1 = new TFile(ofile.c_str(),"RECREATE");
    
    std::string nov_run = "runs34n"; //the run that needs to be scaled

    int tot_ent=0;
    for (size_t i=0;i<num_ent.size();i++){
        tot_ent += num_ent[i];
    }
    int hitPatternBins = nbins_hitpattern;
    TH2D *hHitPattern = new TH2D("hHitPattern", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern = new TH2D("hBarHitPattern", "", hitPatternBins, -1000, 1000, 72, 0, 72);
    hHitPattern->Sumw2();
    hBarHitPattern->Sumw2();

    TH2D *hHitPattern_rnd = (TH2D*)hHitPattern->Clone("hHitPattern_rnd");
    TH2D *hBarHitPattern_rnd = (TH2D*)hBarHitPattern->Clone("hBarHitPattern_rnd");


    TH2D *hHitPatternLayer[3];

	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i] =
			new TH2D(Form("hHitPatternLayer%d", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer[i]->Sumw2();
	}
    
	auto hHitPattern_exc = new TH2D("hHitPattern_exc", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern_exc = new TH2D("hBarHitPattern_exc", "", hitPatternBins, -1000, 1000, 72, 0, 72);
	hHitPattern_exc->Sumw2();
    hBarHitPattern_exc->Sumw2();
    
    TH2D *hHitPatternLayer_exc[3];
    TH2D *hHitPatternLayer_rnd[3];

	for (int i = 0; i < 3; i++) {
		hHitPatternLayer_exc[i] = new TH2D(Form("hHitPatternLayer%d_exc", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer_rnd[i] = new TH2D(Form("hHitPatternLayer%d_rnd", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
        
        hHitPatternLayer_exc[i]->Sumw2();
        hHitPatternLayer_rnd[i]->Sumw2();
    }

	TString htitlex = "Normalized #frac{dN}{dx} vs x; x [mm]; #frac{dN}{dx} [mm^{-1}]";
	TString htitley = "Normalized #frac{dN}{dy} vs y; y [mm]; #frac{dN}{dy} [mm^{-1}]";
	TString htitleAx = "Normalized #frac{dN}{dA} vs x; x [mm]; #frac{dN}{dA} [mm^{-2}]";
	TString htitleAy = "Normalized #frac{dN}{dA} vs y; y [mm]; #frac{dN}{dA} [mm^{-2}]";


    TH1D *hBarHit_all = new TH1D("hBarHit_all", "", hitPatternBins, -1000, 1000);
    hBarHit_all->SetTitle(htitlex);
    hBarHit_all->Sumw2();
    TH1D *hBarHit_all_rb = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb");
    TH1D *hBarHit_all_rb2 = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb2");
    TH1D *hBarHit_all_rb3 = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb3");

    hBarHit_all_rb->Rebin(rbf);
    hBarHit_all_rb2->Rebin(rbf*2);
    hBarHit_all_rb3->Rebin(rbf*4);

    TH1D *hBarHit_all_rnd = (TH1D*)hBarHit_all->Clone("hBarHit_all_rnd");
    TH1D *hBarHit_all_rnd_rb = (TH1D*)hBarHit_all_rb->Clone("hBarHit_all_rnd_rb");
    TH1D *hBarHit_all_rnd_rb2 = (TH1D*)hBarHit_all_rb2->Clone("hBarHit_all_rnd_rb2");
    TH1D *hBarHit_all_rnd_rb3 = (TH1D*)hBarHit_all_rb3->Clone("hBarHit_all_rnd_rb3");


    TH1D *hBarHit_all_area = (TH1D*)hBarHit_all->Clone("hBarHit_all_area");
    TH1D *hBarHit_all_area_rb = (TH1D*)hBarHit_all_rb->Clone("hBarHit_all_area_rb");
    TH1D *hBarHit_all_area_rb2 = (TH1D*)hBarHit_all_rb2->Clone("hBarHit_all_area_rb2");
    TH1D *hBarHit_all_area_rb3 = (TH1D*)hBarHit_all_rb3->Clone("hBarHit_all_area_rb3");

    TH1D *hBarHit_all_rnd_area = (TH1D*)hBarHit_all_rnd->Clone("hBarHit_all_rnd_area");
    TH1D *hBarHit_all_rnd_area_rb = (TH1D*)hBarHit_all_rnd_rb->Clone("hBarHit_all_rnd_area_rb");
    TH1D *hBarHit_all_rnd_area_rb2 = (TH1D*)hBarHit_all_rnd_rb2->Clone("hBarHit_all_rnd_area_rb2");
    TH1D *hBarHit_all_rnd_area_rb3 = (TH1D*)hBarHit_all_rnd_rb3->Clone("hBarHit_all_rnd_area_rb3");



	TH1D *hBarHit[72];
	TH1D *hBarHit_exc[72];
    TH1D *hBarHit_rb[72];
    TH1D *hBarHit_rb2[72];
    TH1D *hBarHit_rb3[72];
    TH1D *hBarHit_area[72];
    TH1D *hBarHit_area_rb[72];
    TH1D *hBarHit_area_rb2[72];
    TH1D *hBarHit_area_rb3[72];
    TH1D *hBarHit_rnd[72];
    TH1D *hBarHit_rnd_rb[72];
    TH1D *hBarHit_rnd_rb2[72];
    TH1D *hBarHit_rnd_rb3[72];
    TH1D *hBarHit_rnd_area[72];
    TH1D *hBarHit_rnd_area_rb[72];
    TH1D *hBarHit_rnd_area_rb2[72];
    TH1D *hBarHit_rnd_area_rb3[72];

    for (int i=0;i<72;i++){
        hBarHit[i] = new TH1D(Form("hBarHit%d",i),"",200,-1000,1000);
        hBarHit[i]->SetTitle(htitlex);
        if (i>=24 && i<48){
            hBarHit[i]->SetTitle(htitley);
        }
        hBarHit[i]->Sumw2();
        hBarHit_rb[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb",i));
        hBarHit_rb2[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb2",i));
        hBarHit_rb3[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb3",i));

        hBarHit_rb[i]->Rebin(rbf);
        hBarHit_rb2[i]->Rebin(rbf*2);
        hBarHit_rb3[i]->Rebin(rbf*4);

        hBarHit_area[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_area",i));
        hBarHit_area_rb[i] = (TH1D*)hBarHit_rb[i]->Clone(Form("hBarHit%d_area_rb",i));
        hBarHit_area_rb2[i] = (TH1D*)hBarHit_rb2[i]->Clone(Form("hBarHit%d_area_rb2",i));
        hBarHit_area_rb3[i] = (TH1D*)hBarHit_rb3[i]->Clone(Form("hBarHit%d_area_rb3",i));


        hBarHit_exc[i] = new TH1D(Form("hBarHit_exc%d",i),"",200,-1000,1000);

        hBarHit_rnd[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rnd",i));
        hBarHit_rnd_rb[i] = (TH1D*)hBarHit_rb[i]->Clone(Form("hBarHit%d_rnd_rb",i));
        hBarHit_rnd_rb2[i] = (TH1D*)hBarHit_rb2[i]->Clone(Form("hBarHit%d_rnd_rb2",i));
        hBarHit_rnd_rb3[i] = (TH1D*)hBarHit_rb3[i]->Clone(Form("hBarHit%d_rnd_rb3",i));

        hBarHit_rnd_area[i] = (TH1D*)hBarHit_rnd[i]->Clone(Form("hBarHit%d_rnd_area",i));
        hBarHit_rnd_area_rb[i] = (TH1D*)hBarHit_rnd_rb[i]->Clone(Form("hBarHit%d_rnd_area_rb",i));
        hBarHit_rnd_area_rb2[i] = (TH1D*)hBarHit_rnd_rb2[i]->Clone(Form("hBarHit%d_rnd_area_rb2",i));
        hBarHit_rnd_area_rb3[i] = (TH1D*)hBarHit_rnd_rb3[i]->Clone(Form("hBarHit%d_rnd_area_rb3",i));

    }

    //convert histograms to counts and add them together
    for (size_t i=0;i<inputs.size();i++){        
        double scl_change = 1.0;

        if (inputs[i].find(nov_run) != std::string::npos){
            scl_change = scale_factor;
        }

        TFile *f = TFile::Open(inputs[i].c_str(),"READ");
        //add the all layers histograms
        std::string hhp = "hHitPattern";
        std::string hbhp = "hBarHitPattern";
        std::string hhpl = "hHitPatternLayer";

        h2_cmb(hHitPattern, f, hhp, scl_change, num_ent[i]);
        h2_cmb(hBarHitPattern, f, hbhp, scl_change, num_ent[i]);

        for (int j=0;j<3;j++){
            //add layer-wise histograms
            h2_cmb(hHitPatternLayer[j], f, hhpl+to_string(j), scl_change, num_ent[i]);
            h2_cmb(hHitPatternLayer_exc[j], f, hhpl+to_string(j) + "_exc", scl_change, num_ent[i]);
            h2_cmb(hHitPatternLayer_rnd[j], f, hhpl+to_string(j) + "_rnd", scl_change, num_ent[i]);
        }

        h2_cmb(hHitPattern_exc, f, hhp + "_exc", scl_change, num_ent[i]);
        h2_cmb(hBarHitPattern_exc, f, hbhp + "_exc", scl_change, num_ent[i]);

        h2_cmb(hHitPattern_rnd, f, hhp + "_rnd", scl_change, num_ent[i]);
        h2_cmb(hBarHitPattern_rnd, f, hbhp + "_rnd", scl_change, num_ent[i]);

        std::string hbha = "hBarHit_all";
        std::string hbhar = "hBarHit_all_rnd";
        std::string hbhaa = "hBarHit_all_area";
        std::string hbhara = "hBarHit_all_rnd_area";
        
        h_cmb(hBarHit_all, f, hbha, scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rb, f, hbha + "_rb", scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rb2, f, hbha + "_rb2", scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rb3, f, hbha + "_rb3", scl_change, num_ent[i]);

        ha_cmb(hBarHit_all_area, f, hbhaa, scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_area_rb, f, hbhaa + "_rb", scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_area_rb2, f, hbhaa + "_rb2", scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_area_rb3, f, hbhaa + "_rb3", scl_change, num_ent[i]);

        h_cmb(hBarHit_all_rnd, f, hbhar, scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rnd_rb, f, hbhar + "_rb", scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rnd_rb2, f, hbhar + "_rb2", scl_change, num_ent[i]);
        h_cmb(hBarHit_all_rnd_rb3, f, hbhar + "_rb3", scl_change, num_ent[i]);

        ha_cmb(hBarHit_all_rnd_area, f, hbhara, scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_rnd_area_rb, f, hbhara + "_rb", scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_rnd_area_rb2, f, hbhara + "_rb2", scl_change, num_ent[i]);
        ha_cmb(hBarHit_all_rnd_area_rb3, f, hbhara + "_rb3", scl_change, num_ent[i]);

        std::string hbh = "hBarHit";

        for (int j=0;j<72;j++){
            h_cmb(hBarHit[j], f, hbh + to_string(j), scl_change, num_ent[i]);
            h_cmb(hBarHit_rb[j], f, hbh + to_string(j) + "_rb", scl_change, num_ent[i]);
            h_cmb(hBarHit_rb2[j], f, hbh + to_string(j) + "_rb2", scl_change, num_ent[i]);
            h_cmb(hBarHit_rb3[j], f, hbh + to_string(j) + "_rb3", scl_change, num_ent[i]);

            ha_cmb(hBarHit_area[j], f, hbh + to_string(j) + "_area" , scl_change, num_ent[i]);
            ha_cmb(hBarHit_area_rb[j], f, hbh + to_string(j) + "_area_rb", scl_change, num_ent[i]);
            ha_cmb(hBarHit_area_rb2[j], f, hbh + to_string(j) + "_area_rb2", scl_change, num_ent[i]);
            ha_cmb(hBarHit_area_rb3[j], f, hbh + to_string(j) + "_area_rb3", scl_change, num_ent[i]);

            h_cmb(hBarHit_exc[j], f, hbh +"_exc" + to_string(j), scl_change, num_ent[i]);

            h_cmb(hBarHit_rnd[j], f, hbh + to_string(j) + "_rnd", scl_change, num_ent[i]);
            h_cmb(hBarHit_rnd_rb[j], f, hbh + to_string(j) + "_rnd_rb", scl_change, num_ent[i]);
            h_cmb(hBarHit_rnd_rb2[j], f, hbh + to_string(j) + "_rnd_rb2", scl_change, num_ent[i]);
            h_cmb(hBarHit_rnd_rb3[j], f, hbh + to_string(j) + "_rnd_rb3", scl_change, num_ent[i]);

            ha_cmb(hBarHit_rnd_area[j], f, hbh + to_string(j) + "_rnd_area", scl_change, num_ent[i]);
            ha_cmb(hBarHit_rnd_area_rb[j], f, hbh + to_string(j) + "_rnd_area_rb", scl_change, num_ent[i]);
            ha_cmb(hBarHit_rnd_area_rb2[j], f, hbh + to_string(j) + "_rnd_area_rb2", scl_change, num_ent[i]);
            ha_cmb(hBarHit_rnd_area_rb3[j], f, hbh + to_string(j) + "_rnd_area_rb3", scl_change, num_ent[i]);
        }
        f->Close();
    }

    //scale to density and write to file
    fo1->cd();
    h2_den_write(hHitPattern, fo1, tot_ent);
    h2_den_write(hBarHitPattern, fo1, tot_ent);
    for (int i=0;i<3;i++){
        h2_den_write(hHitPatternLayer[i], fo1, tot_ent);
        h2_den_write(hHitPatternLayer_exc[i], fo1, tot_ent);
        h2_den_write(hHitPatternLayer_rnd[i], fo1, tot_ent);
    }

    h2_den_write(hHitPattern_exc, fo1, tot_ent);
    h2_den_write(hBarHitPattern_exc, fo1, tot_ent);


    h2_den_write(hHitPattern_rnd, fo1, tot_ent);
    h2_den_write(hBarHitPattern_rnd, fo1, tot_ent);

    h_den_write(hBarHit_all, fo1, tot_ent);
    h_den_write(hBarHit_all_rb, fo1, tot_ent);
    h_den_write(hBarHit_all_rb2, fo1, tot_ent);
    h_den_write(hBarHit_all_rb3, fo1, tot_ent);

    ha_den_write(hBarHit_all_area, fo1, tot_ent);
    ha_den_write(hBarHit_all_area_rb, fo1, tot_ent);
    ha_den_write(hBarHit_all_area_rb2, fo1, tot_ent);
    ha_den_write(hBarHit_all_area_rb3, fo1, tot_ent);

    h_den_write(hBarHit_all_rnd, fo1, tot_ent);
    h_den_write(hBarHit_all_rnd_rb, fo1, tot_ent);
    h_den_write(hBarHit_all_rnd_rb2, fo1, tot_ent);
    h_den_write(hBarHit_all_rnd_rb3, fo1, tot_ent);

    ha_den_write(hBarHit_all_rnd_area, fo1, tot_ent);
    ha_den_write(hBarHit_all_rnd_area_rb, fo1, tot_ent);
    ha_den_write(hBarHit_all_rnd_area_rb2, fo1, tot_ent);
    ha_den_write(hBarHit_all_rnd_area_rb3, fo1, tot_ent);


    for (int j=0;j<72;j++){
        h_den_write(hBarHit[j], fo1, tot_ent);
        h_den_write(hBarHit_rb[j], fo1, tot_ent);
        h_den_write(hBarHit_rb2[j], fo1, tot_ent);
        h_den_write(hBarHit_rb3[j], fo1, tot_ent);
        ha_den_write(hBarHit_area[j], fo1, tot_ent);
        ha_den_write(hBarHit_area_rb[j], fo1, tot_ent);
        ha_den_write(hBarHit_area_rb2[j], fo1, tot_ent);
        ha_den_write(hBarHit_area_rb3[j], fo1, tot_ent);

        h_den_write(hBarHit_exc[j], fo1, tot_ent);

        h_den_write(hBarHit_rnd[j], fo1, tot_ent);
        h_den_write(hBarHit_rnd_rb[j], fo1, tot_ent);
        h_den_write(hBarHit_rnd_rb2[j], fo1, tot_ent);
        h_den_write(hBarHit_rnd_rb3[j], fo1, tot_ent);

        ha_den_write(hBarHit_rnd_area[j], fo1, tot_ent);
        ha_den_write(hBarHit_rnd_area_rb[j], fo1, tot_ent);
        ha_den_write(hBarHit_rnd_area_rb2[j], fo1, tot_ent);
        ha_den_write(hBarHit_rnd_area_rb3[j], fo1, tot_ent);

    }


    fo1->Close();
    return;
}