#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TLine.h>


void visualize_simple(
    const std::string &file = "input.root",
    const std::string &ofile = "output.root"
){
    auto *f = new TFile(file.c_str(),"READ");
    auto *fo1 = new TFile(ofile.c_str(),"RECREATE");

    std::vector<std::vector<double>> x = {{pro_0[0],deu_0[0],tri_0[0],hel_0[0],alp_0[0]},{pro_1[0],deu_1[0],tri_1[0],hel_1[0],alp_1[0]},{pro_2[0],deu_2[0],tri_2[0],hel_2[0],alp_2[0]}};
    std::vector<double> y = {pro_pt,deu_pt,tri_pt,hel_pt,alp_pt};
    for (int i=0;i<3;i++){
        TCanvas *c = new TCanvas(Form("c%d",i),Form("Layer %d",i));
        c->cd();
        TH2D *h = (TH2D*)f->Get(Form("hEVsTof_L%d_cal",i));
        TGraph *g = new TGraph(5,x[i].data(),y.data());
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kRed);
        h->Draw();
        g->Draw("P SAME");
        fo1->cd();
        c->Write();

        TCanvas *c2 = new TCanvas(Form("c%d_v2",i),Form("Layer %d",i));
        c2->cd();
        TH2D *h2 = (TH2D*)f->Get(Form("hEVsTof_L%d_cal_v2",i));
        h2->Draw();
        g->Draw("P SAME");
        fo1->cd();
        c2->Write();

        for (int j=1;j<=5;j++){
            TCanvas *c_proj = new TCanvas(Form("c_proj%d_L%d",j,i),Form("Layer %d Projections",i));
            c_proj->cd();
            int xbin = h->GetXaxis()->FindBin(x[i][j-1]);
            TH1D *py = h->ProjectionY(Form("proj_%d",j),xbin-1,xbin+1);
            py->Draw();
            double gmin = gPad->GetUymin();
            double gmax = gPad->GetUymax();
            TLine *line = new TLine(y[j-1],gmin,y[j-1],gmax);
            line->SetLineWidth(2);
            line->Draw("SAME");
            fo1->cd();
            c_proj->Write();
        }



        for (int j=1;j<=5;j++){
            TCanvas *c_proj = new TCanvas(Form("c_proj%d_L%d_v2",j,i),Form("Layer %d Projections",i));
            c_proj->cd();
            int xbin = h2->GetXaxis()->FindBin(x[i][j-1]);
            TH1D *py = h2->ProjectionY(Form("proj_%d",j),xbin-1,xbin+1);
            py->Draw();
            double gmin = gPad->GetUymin();
            double gmax = gPad->GetUymax();
            TLine *line = new TLine(y[j-1],gmin,y[j-1],gmax);
            line->SetLineWidth(2);
            line->Draw("SAME");
            fo1->cd();
            c_proj->Write();
        }
    
        /* TCanvas *c_proj2 = new TCanvas(Form("c_pro%d_v2",i),Form("Layer %d Projections",i));
        c_proj2->cd();
        c_proj2->Divide(3,1);
        for (int j=1;j<=3;j++){
            c_proj2->cd(j);
            int xbin = h2->GetXaxis()->FindBin(x[i][j-1]);
            TH1D *py = h2->ProjectionY(Form("proj_%d",j),xbin-1,xbin+1);
            py->Draw();
            double gmin = gPad->GetUymin();
            double gmax = gPad->GetUymax();
            TLine *line = new TLine(y[j-1],gmin,y[j-1],gmax);
            line->SetLineWidth(2);
            line->Draw("SAME");
        }
        fo1->cd();
        c_proj2->Write(); */
    
    }

    fo1->Close();

    return;
}