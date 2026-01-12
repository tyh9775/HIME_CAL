#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>


void cnv_write(TFile *f,TFile *fo1,std::string hname,std::string hname2,double x1,double x2){
    TCanvas *c1 = new TCanvas(hname.c_str(),hname.c_str());
    c1->cd();
    c1->SetSupportGL(true);

    TH1D *h1 = (TH1D*)f->Get(hname.c_str())->Clone();
    h1->SetLineColor(kRed);
    h1->SetLineWidth(2);
    h1->SetMaximum(h1->GetMaximum()*1.2);
    h1->Draw("E1 HIST");

    TH1D *h2 = (TH1D*)f->Get(hname2.c_str())->Clone();
    h2->SetLineColor(kGray);
    h2->SetLineWidth(2);
    h2->SetLineStyle(7);
    h2->Draw("E1 HIST SAME");
    
    gPad->Update();
    double max_y = gPad->GetFrame()->GetY2();
    double min_y = gPad->GetFrame()->GetY1();

    TLine *vline1 = new TLine(x1,min_y,x1,max_y);
    TLine *vline2 = new TLine(x2,min_y,x2,max_y);
    vline1->Draw("SAME");
    vline2->Draw("SAME");

    TLegend *legend1 = new TLegend(0.78, 0.6, 0.98, 0.75);
    legend1->AddEntry(h1,"shadow","l");
    legend1->AddEntry(h2,"no shadow","l");
    legend1->Draw();
    fo1->cd();
    c1->SetName(hname.c_str());
    c1->Write();

    return;
}

void g_write(TFile *f,TFile *fo1,std::vector<std::string> gname, TCanvas *c1,int sb_num){

    TLegend *legend1 = new TLegend(0.12, 0.72, 0.35, 0.88);

    for (size_t i=0;i<gname.size();i++){
        TGraph *g1 = (TGraph*)f->Get(gname[i].c_str())->Clone();
        g1->GetXaxis()->SetTitle("Momentum (MeV/c)");
        g1->GetYaxis()->SetTitle("Yield Fraction");
        g1->SetTitle(Form("Shadow Bar %d Yield Fraction vs Momentum", sb_num+2));
        g1->GetYaxis()->SetRangeUser(0.0,1.0);
        g1->GetXaxis()->SetLimits(0.0,500.0);
        g1->SetMarkerStyle(20);
        g1->SetMarkerSize(0.8);
        g1->SetMarkerColor(kRed+i*2);
        c1->cd();
        if (i==0){
            g1->Draw("AP");
        } else {
            g1->Draw("P SAME");
        }

        legend1->AddEntry(g1,Form("Layer %zu", i+1),"p");
    }
    c1->cd();
    legend1->Draw();

    fo1->cd();
    c1->Write();
    return;
}


void visualize_comp_v4(
    const std::string &file = "input1.root",
    const std::string &ofile = "output.root"
){
    gStyle->SetCanvasPreferGL(true);
    std::string hb = "hBarHit_mnt";
    //std::string hbe = "hBarHit_exc";
    TFile *f = TFile::Open(file.c_str(),"READ");

    std::vector<double> sb2_x {};
    std::vector<double> sb2_y {};
    std::vector<double> sb3_x {};
    std::vector<double> sb3_y {};
    std::vector<double> sb4_x {};
    std::vector<double> sb4_y {};

    TFile *fo1 = TFile::Open(ofile.c_str(),"RECREATE");
    std::vector<std::vector<int>> sbs {sb2,sb3,sb4};
    std::vector<double> sb_xmins {sb2_xmin,sb3_xmin,sb4_xmin};
    std::vector<double> sb_xmaxes {sb2_xmax,sb3_xmax,sb4_xmax};
    std::vector<double> sb_ymins {sb2_ymin,sb3_ymin,sb4_ymin};
    std::vector<double> sb_ymaxes {sb2_ymax,sb3_ymax,sb4_ymax};

    for (int i=0;i<3;i++){
        std::vector<std::string> gnames {};
        std::vector<std::string> gnames_rb {};        
        for (int j=0;j<3;j++){
            double xmin = 0.0;
            double xmax = 0.0;
            if (j==0 || j==2){
                xmin = sb_xmins[i];
                xmax = sb_xmaxes[i];
            } else {
                xmin = sb_ymins[i];
                xmax = sb_ymaxes[i];
            }


            for (int k=0;k<5;k++){
                std::string hname = hb+to_string(sbs[i][j])+"_sb_"+to_string(i+2)+"_sec_"+to_string(k);
                std::string hname2 = hb+to_string(sbs[i][j])+"_no_sb_"+to_string(i+2)+"_sec_"+to_string(k);

                std::string hnr = hb+to_string(sbs[i][j])+"_rb_sb_"+to_string(i+2)+"_sec_"+to_string(k);
                std::string hn2r = hb+to_string(sbs[i][j])+"_rb_no_sb_"+to_string(i+2)+"_sec_"+to_string(k);

                cnv_write(f,fo1,hname,hname2,xmin,xmax);
                cnv_write(f,fo1,hnr,hn2r,xmin,xmax);

            }

            gnames.push_back("hBarHit_mnt"+to_string(sbs[i][j])+"_yield_coarse_sb_"+to_string(i+2));
            gnames_rb.push_back("hBarHit_mnt"+to_string(sbs[i][j])+"_rb_yield_coarse_sb_"+to_string(i+2));

        }
        TCanvas *c1 = new TCanvas(("c_yield_sb_"+to_string(i+2)).c_str(),("c_yield_sb_"+to_string(i+2)).c_str());
        TCanvas *c1_rb = new TCanvas(("c_yield_rb_sb_"+to_string(i+2)).c_str(),("c_yield_rb_sb_"+to_string(i+2)).c_str());
        c1->SetTitle(Form("SB%d Yield Fraction vs Momentum", i+2));
        c1_rb->SetTitle(Form("SB%d RB Yield Fraction vs Momentum", i+2));

        g_write(f,fo1,gnames,c1,i);
        g_write(f,fo1,gnames_rb,c1_rb,i);
    }
    f->Close();
    fo1->Close();
    return;
}
