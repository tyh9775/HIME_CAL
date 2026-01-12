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


void visualize_comp_v3(
    const std::string &file = "input1.root",
    const std::string &ofile = "output.root"
){
    gStyle->SetCanvasPreferGL(true);
    std::string hb = "hBarHit";
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
            std::string hname = hb+to_string(sbs[i][j])+"_sb_"+to_string(i+2);
            std::string hname2 = hb+to_string(sbs[i][j])+"_no_sb_"+to_string(i+2);

            std::string hnr = hb+to_string(sbs[i][j])+"_rb_sb_"+to_string(i+2);
            std::string hn2r = hb+to_string(sbs[i][j])+"_rb_no_sb_"+to_string(i+2);

            std::string hnr2 = hb+to_string(sbs[i][j])+"_rb2_sb_"+to_string(i+2);
            std::string hn2r2 = hb+to_string(sbs[i][j])+"_rb2_no_sb_"+to_string(i+2);

            std::string hnr3 = hb+to_string(sbs[i][j])+"_rb3_sb_"+to_string(i+2);
            std::string hn2r3 = hb+to_string(sbs[i][j])+"_rb3_no_sb_"+to_string(i+2);

            cnv_write(f,fo1,hname,hname2,xmin,xmax);
            cnv_write(f,fo1,hnr,hn2r,xmin,xmax);
            cnv_write(f,fo1,hnr2,hn2r2,xmin,xmax);
            cnv_write(f,fo1,hnr3,hn2r3,xmin,xmax);


            std::string hna = hb+to_string(sbs[i][j])+"_area_sb_"+to_string(i+2);
            std::string hna2 = hb+to_string(sbs[i][j])+"_area_no_sb_"+to_string(i+2);

            std::string hna_r = hb+to_string(sbs[i][j])+"_area_rb_sb_"+to_string(i+2);
            std::string hna2_r = hb+to_string(sbs[i][j])+"_area_rb_no_sb_"+to_string(i+2);

            std::string hna_r2 = hb+to_string(sbs[i][j])+"_area_rb2_sb_"+to_string(i+2);
            std::string hna2_r2 = hb+to_string(sbs[i][j])+"_area_rb2_no_sb_"+to_string(i+2);

            std::string hna_r3 = hb+to_string(sbs[i][j])+"_area_rb3_sb_"+to_string(i+2);
            std::string hna2_r3 = hb+to_string(sbs[i][j])+"_area_rb3_no_sb_"+to_string(i+2);

            cnv_write(f,fo1,hna,hna2,xmin,xmax);
            cnv_write(f,fo1,hna_r,hna2_r,xmin,xmax);
            cnv_write(f,fo1,hna_r2,hna2_r2,xmin,xmax);
            cnv_write(f,fo1,hna_r3,hna2_r3,xmin,xmax);

            std::string hn_rnd = hb+to_string(sbs[i][j])+"_rnd_sb_"+to_string(i+2);
            std::string hn2_rnd = hb+to_string(sbs[i][j])+"_rnd_no_sb_"+to_string(i+2);

            std::string hn_rnd_r = hb+to_string(sbs[i][j])+"_rnd_rb_sb_"+to_string(i+2);
            std::string hn2_rnd_r = hb+to_string(sbs[i][j])+"_rnd_rb_no_sb_"+to_string(i+2);

            std::string hn_rnd2 = hb+to_string(sbs[i][j])+"_rnd_rb2_sb_"+to_string(i+2);
            std::string hn2_rnd2 = hb+to_string(sbs[i][j])+"_rnd_rb2_no_sb_"+to_string(i+2);

            std::string hn_rnd3 = hb+to_string(sbs[i][j])+"_rnd_rb3_sb_"+to_string(i+2);
            std::string hn2_rnd3 = hb+to_string(sbs[i][j])+"_rnd_rb3_no_sb_"+to_string(i+2);

            cnv_write(f,fo1,hn_rnd,hn2_rnd,xmin,xmax);
            cnv_write(f,fo1,hn_rnd_r,hn2_rnd_r,xmin,xmax);
            cnv_write(f,fo1,hn_rnd2,hn2_rnd2,xmin,xmax);
            cnv_write(f,fo1,hn_rnd3,hn2_rnd3,xmin,xmax);

        }
    }
    f->Close();
    fo1->Close();
    return;
}
