#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>


void canvas_write(TCanvas *c,TH1D *h1,TH1D *h2,TFile *f_out,std::string l_name1,std::string l_name2){
    c->cd();
    c->SetSupportGL(true);

    if (h1->GetMaximum()>h2->GetMaximum()){
        h1->SetMaximum(h1->GetMaximum()*1.1);
    } else {
        h1->SetMaximum(h2->GetMaximum()*1.1);
    }
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->Draw("E1 HIST");
    h2->SetLineColor(kGray);
    h2->SetLineWidth(2);
    h2->Draw("E1 HIST SAME");
    TLegend *legend = new TLegend(0.78, 0.6, 0.98, 0.75);
    legend->AddEntry(h1,l_name1.c_str(),"l");
    legend->AddEntry(h2,l_name2.c_str(),"l");
    legend->Draw();
    
    f_out->cd();
    c->Write();
    return;
}

double find_ratio(TH1D *h1, TH1D *h2){
    double ratio = h1->Integral()/h2->Integral();
    return ratio;
}

double find_ratio_lim(TH1D *h1, TH1D *h2, int i_start, int i_stop){
    double ratio = h1->Integral(i_start,i_stop)/h2->Integral(i_start,i_stop);
    return ratio;
}

void scale_check_v2(
    const std::vector<std::string> &inputs = {"june.root,nov.root"},
    const std::string &ofile = "output"
){
    TFile *fo1 = new TFile((ofile+".root").c_str(),"RECREATE");
    std::ofstream output;
    output.open((ofile+".txt").c_str(),std::ios_base::trunc);


    TFile *f1 = TFile::Open(inputs[0].c_str(),"READ");
    TFile *f2 = TFile::Open(inputs[1].c_str(),"READ");
    
    TH1D *hBarHit_1 = (TH1D*)f1->Get("hBarHit_all")->Clone("hBarHit_Jun");
    TH1D *hBarHit_1_rb = (TH1D*)f1->Get("hBarHit_all_rb")->Clone("hBarHit_Jun_rb");
    TH1D *hBarHit_1_rb2 = (TH1D*)f1->Get("hBarHit_all_rb2")->Clone("hBarHit_Jun_rb2");
    TH1D *hBarHit_1_rb3 = (TH1D*)f1->Get("hBarHit_all_rb3")->Clone("hBarHit_Jun_rb3");

    TH1D *hBarHit_1_rnd = (TH1D*)f1->Get("hBarHit_all_rnd")->Clone("hBarHit_Jun_rnd");
    TH1D *hBarHit_1_rnd_rb = (TH1D*)f1->Get("hBarHit_all_rnd_rb")->Clone("hBarHit_Jun_rnd_rb");
    TH1D *hBarHit_1_rnd_rb2 = (TH1D*)f1->Get("hBarHit_all_rnd_rb2")->Clone("hBarHit_Jun_rnd_rb2");
    TH1D *hBarHit_1_rnd_rb3 = (TH1D*)f1->Get("hBarHit_all_rnd_rb3")->Clone("hBarHit_Jun_rnd_rb3");


    TH1D *hBarHit_2 = (TH1D*)f2->Get("hBarHit_all")->Clone("hBarHit_Nov");
    TH1D *hBarHit_2_rb = (TH1D*)f2->Get("hBarHit_all_rb")->Clone("hBarHit_Nov_rb");
    TH1D *hBarHit_2_rb2 = (TH1D*)f2->Get("hBarHit_all_rb2")->Clone("hBarHit_Nov_rb2");
    TH1D *hBarHit_2_rb3 = (TH1D*)f2->Get("hBarHit_all_rb3")->Clone("hBarHit_Nov_rb3");

    TH1D *hBarHit_2_rnd = (TH1D*)f2->Get("hBarHit_all_rnd")->Clone("hBarHit_Nov_rnd");
    TH1D *hBarHit_2_rnd_rb = (TH1D*)f2->Get("hBarHit_all_rnd_rb")->Clone("hBarHit_Nov_rnd_rb");
    TH1D *hBarHit_2_rnd_rb2 = (TH1D*)f2->Get("hBarHit_all_rnd_rb2")->Clone("hBarHit_Nov_rnd_rb2");
    TH1D *hBarHit_2_rnd_rb3 = (TH1D*)f2->Get("hBarHit_all_rnd_rb3")->Clone("hBarHit_Nov_rnd_rb3");


    int n_bins = hBarHit_1->GetNbinsX();

    double avg_r=find_ratio(hBarHit_1,hBarHit_2);
    double avg_ri=find_ratio(hBarHit_2,hBarHit_1);


    std::cout<<avg_r<<"\n";
    std::cout<<avg_ri<<"\n";
    std::cout<<"\n";

    std::cout<<find_ratio(hBarHit_1_rb,hBarHit_2_rb)<<"\n";
    std::cout<<find_ratio(hBarHit_1_rb2,hBarHit_2_rb2)<<"\n";
    std::cout<<find_ratio(hBarHit_1_rb3,hBarHit_2_rb3)<<"\n";
    std::cout<<"\n";
    std::cout<<find_ratio(hBarHit_1_rnd,hBarHit_2_rnd)<<"\n";
    std::cout<<find_ratio(hBarHit_1_rnd_rb,hBarHit_2_rnd_rb)<<"\n";
    std::cout<<find_ratio(hBarHit_1_rnd_rb2,hBarHit_2_rnd_rb2)<<"\n";
    std::cout<<find_ratio(hBarHit_1_rnd_rb3,hBarHit_2_rnd_rb3)<<"\n";
    std::cout<<"\n";

    int i_start = hBarHit_1->FindFirstBinAbove(0.0);
    int i_stop = hBarHit_1->FindLastBinAbove(0.0);


    double avg_rl=find_ratio_lim(hBarHit_1,hBarHit_2,i_start,i_stop);
    double avg_rli=find_ratio_lim(hBarHit_2,hBarHit_1,i_start,i_stop);

    std::cout<<avg_rl<<"\n";
    std::cout<<avg_rli<<"\n";
    std::cout<<"\n";

    std::cout<<find_ratio_lim(hBarHit_1_rb,hBarHit_2_rb,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rb2,hBarHit_2_rb2,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rb3,hBarHit_2_rb3,i_start,i_stop)<<"\n";
    std::cout<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rnd,hBarHit_2_rnd,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rnd_rb,hBarHit_2_rnd_rb,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rnd_rb2,hBarHit_2_rnd_rb2,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_1_rnd_rb3,hBarHit_2_rnd_rb3,i_start,i_stop)<<"\n";
    std::cout<<"\n";

    std::cout<<find_ratio_lim(hBarHit_2_rnd,hBarHit_1_rnd,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_2_rnd_rb,hBarHit_1_rnd_rb,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_2_rnd_rb2,hBarHit_1_rnd_rb2,i_start,i_stop)<<"\n";
    std::cout<<find_ratio_lim(hBarHit_2_rnd_rb3,hBarHit_1_rnd_rb3,i_start,i_stop)<<"\n";
    std::cout<<"\n";


    output<<avg_ri<<"\n";
    output<<find_ratio_lim(hBarHit_2_rnd,hBarHit_1_rnd,i_start,i_stop)<<"\n";

    fo1->cd();
    hBarHit_1->Write();
    hBarHit_2->Write();

    hBarHit_1_rb->Write();
    hBarHit_2_rb->Write();
    hBarHit_1_rb2->Write();
    hBarHit_2_rb2->Write();
    hBarHit_1_rb3->Write();
    hBarHit_2_rb3->Write();

    hBarHit_1_rnd->Write();
    hBarHit_2_rnd->Write();
    hBarHit_1_rnd_rb->Write();
    hBarHit_2_rnd_rb->Write();
    hBarHit_1_rnd_rb2->Write();
    hBarHit_2_rnd_rb2->Write();
    hBarHit_1_rnd_rb3->Write();
    hBarHit_2_rnd_rb3->Write();


    std::string j1 = "June";
    std::string n1 = "November";
    std::string j2 = "June_scaled";
    std::string n2 = "November_scaled";

    std::vector<TH1D*> h_jn {hBarHit_1,hBarHit_1_rb,hBarHit_1_rb2,hBarHit_1_rb3,hBarHit_1_rnd,hBarHit_1_rnd_rb,hBarHit_1_rnd_rb2,hBarHit_1_rnd_rb3};
    std::vector<TH1D*> h_nv {hBarHit_2,hBarHit_2_rb,hBarHit_2_rb2,hBarHit_2_rb3,hBarHit_2_rnd,hBarHit_2_rnd_rb,hBarHit_2_rnd_rb2,hBarHit_2_rnd_rb3};

    for (size_t i=0;i<h_jn.size();i++){
        TCanvas *c = new TCanvas(Form("c%zu",i+1),Form("c%zu",i+1));
        canvas_write(c,h_jn[i],h_nv[i],fo1,j1,n1);

        TCanvas *c_scl = new TCanvas(Form("c_scl_%zu",i+1),Form("c_scl_%zu",i+1));
        TH1D *hcopy_j = (TH1D*)h_jn[i]->Clone(Form("%s_scaled",h_jn[i]->GetName()));
        hcopy_j->Scale(avg_ri);
        canvas_write(c_scl,hcopy_j,h_nv[i],fo1,j2,n1);

        TCanvas *c_scl_rev = new TCanvas(Form("c_scl_rev_%zu",i+1),Form("c_scl_rev_%zu",i+1));
        TH1D *hcopy_n = (TH1D*)h_nv[i]->Clone(Form("%s_scaled",h_nv[i]->GetName()));
        hcopy_n->Scale(avg_r);
        canvas_write(c_scl_rev,h_jn[i],hcopy_n,fo1,j1,n2);

    }



    f1->Close();
    f2->Close();

    fo1->Close();

    return;
}


