#include "myConst3.h"
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

std::vector<double> cut_range(const TCutG *c){
    double xmax = c->GetPointX(0);
    double xmin = c->GetPointX(0);
    double ymax = c->GetPointY(0);
    double ymin = c->GetPointY(0);
    for (int i=1;i<c->GetN();i++){
        double xval = c->GetPointX(i);
        if (xval<xmin){
            xmin=xval;
        } else if (xval>xmax){
            xmax=xval;
        }
        double yval = c->GetPointY(i);
        if (yval<ymin){
            ymin=yval;
        } else if (yval>ymax){
            ymax=yval;
        }        
    }
    return {xmin,xmax,ymin,ymax};
}

TH2D *htmp_make(const TH2D *h,const TCutG *cut){
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    TH2D *htmp = (TH2D*)h->Clone();
    for (int xi=1;xi<=nx;xi++){
        double xval = h->GetXaxis()->GetBinCenter(xi);
        for (int yi=1;yi<=ny;yi++){
            double yval = h->GetYaxis()->GetBinCenter(yi);
            if (!cut->IsInside(xval,yval)){
                htmp->SetBinContent(xi,yi,0);
                htmp->SetBinError(xi,yi,0);
            }
        }   
    }

    htmp->SetEntries(htmp->Integral(1,nx,1,ny));

    return htmp;
}


std::vector<std::vector<double>> data_load(const TH2D *h,const TCutG *cut){
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> err {};
    
    TH2D *htmp = htmp_make(h,cut);

    int xmin_i = htmp->GetXaxis()->FindFixBin(cut_range(cut)[0]);
    int xmax_i = htmp->GetXaxis()->FindFixBin(cut_range(cut)[1]);
    int ymin_i = htmp->GetYaxis()->FindFixBin(cut_range(cut)[2]);
    int ymax_i = htmp->GetYaxis()->FindFixBin(cut_range(cut)[3]);

    int nybins = htmp->GetNbinsY();
    for (int xj=xmin_i;xj<xmax_i;xj++){
        TH1D *py_tmp = htmp->ProjectionY(Form("py%d",xj),xj-1,xj+1);
        double xval=htmp->GetXaxis()->GetBinCenter(xj);
        for (int yj=1;yj<=nybins;yj++){
            double yval=htmp->GetYaxis()->GetBinCenter(yj);
            if (!cut->IsInside(xval,yval)){
                py_tmp->SetBinContent(yj,0);
                py_tmp->SetBinError(yj,0);
            }
        }
        
        py_tmp->SetEntries(py_tmp->Integral(ymin_i,ymax_i));
        
        if (py_tmp->GetEntries()==0){
            continue;
        }
        int max_ybin = py_tmp->GetMaximumBin();
        double max_y = py_tmp->GetBinCenter(max_ybin);
        double max_count = py_tmp->GetBinContent(max_ybin);

        int y_low = py_tmp->FindFirstBinAbove(0);
        int y_high = py_tmp->FindLastBinAbove(0);
        double y_l = py_tmp->GetBinCenter(y_low);
        double y_h = py_tmp->GetBinCenter(y_high);

        TF1 *f_1 = new TF1("f_1","gaus",y_l,y_h);
        f_1->SetParameters(max_count,max_y,1);
        f_1->SetParLimits(1,y_l,y_h);
        f_1->SetParLimits(2,0.1,50.0);
        py_tmp->Fit(f_1,"QRB+");
        double max_y_fit = f_1->GetParameter(1);
        double y_err = f_1->GetParError(1);

        if (/* abs(y_err)<=2.5 &&  */cut->IsInside(xval,max_y_fit)){
            x.push_back(xval);
            y.push_back(max_y_fit);
            err.push_back(y_err);
            
        } 
        py_tmp->Delete();
    }
    

    return {x,y,err};
}

std::vector<double> max_find(const TGraphErrors *graph){
    TGraphErrors *g = (TGraphErrors*)graph->Clone();
    int n = g->GetN();
    double x_max {0.0};
    double y_max {0.0};
    double err_max {0.0};
    double x0 = g->GetPointX(0);
    double xf = g->GetPointX(n-1);
    double y0 = g->GetPointY(0);
    double yf = g->GetPointY(n-1);
    for (int i=0;i<n;i++){
        double xi = g->GetPointX(i);
        double yi = g->GetPointY(i);
        double err_i = g->GetErrorY(i);
        /* if (err_i>2.5 || err_i<0.0){
            continue;
        } */
        if (yi>=y_max){
            y_max=yi;
            x_max=xi;
            err_max=err_i;
        }
    }
    return {x_max,y_max,err_max,x0,y0,xf,yf};
}


void reader(
    const std::string &file = "input.root",
    const std::string &ofile = "output.root",
    const std::vector<std::string> &cuts = {"cut1.root","cut2.root","cut3.root"},
    const std::vector<std::string> &cuts2 = {"cut12.root","cut22.root","cut32.root"}
){
    TFile *f = new TFile(file.c_str(),"READ");
    TFile *fo1 = new TFile(ofile.c_str(),"RECREATE");

    std::vector<std::vector<double>> x = {{pro_0[0],deu_0[0],tri_0[0],hel_0[0],alp_0[0]},{pro_1[0],deu_1[0],tri_1[0],hel_1[0],alp_1[0]},{pro_2[0],deu_2[0],tri_2[0],hel_2[0],alp_2[0]}};
    std::vector<double> y = {pro_pt,deu_pt,tri_pt,hel_pt,alp_pt};

    for (int i=0;i<3;i++){
        TH2D *h = (TH2D*)f->Get(Form("hTotVsTof_sh_t_L%d",i));
        TH2D *hc = (TH2D*)h->Clone(Form("hTotVsTof_L%d_cal",i));
        TFile *fcut = new TFile(cuts[i].c_str(),"READ");
        TCanvas *ci = new TCanvas(Form("hTotVsTof_L%d_cal",i));

        std::vector<double> tof {};
        std::vector<TGraphErrors*> graphs {};
        std::vector<TGraphErrors*> mxpts {};

        std::vector<double> tot_mx {};

        for (int j=1;j<=5;j++){
            TCutG *cut = (TCutG*)fcut->Get(Form("cut_s%d",j));
            std::vector<std::vector<double>> data = data_load(hc,cut);
            std::vector<double> x_j = data[0];
            std::vector<double> y_j = data[1];
            std::vector<double> err_j = data[2];

            TGraphErrors *g_j = new TGraphErrors(x_j.size(),x_j.data(),y_j.data(),nullptr,err_j.data());
            g_j->SetName(Form("g_s%d_L%d",j,i));
            std::vector<double> max_vals = max_find(g_j);
            double tof_j = max_vals[0];
            double tot_j = max_vals[1];
            double e_j = max_vals[2];
            TGraphErrors *mxpt = new TGraphErrors(1,&tof_j,&tot_j,nullptr,&e_j);
            mxpt->SetName(Form("mxpt_s%d_L%d",j,i));
            mxpt->SetMarkerStyle(20);
            mxpt->SetMarkerColor(kRed);

            graphs.push_back(g_j);
            //mxpts.push_back(mxpt);

            tof.push_back(tof_j);
            tot_mx.push_back(tot_j);
        }

        ci->cd();
        hc->Draw("COLZ");
        for (size_t gi=0;gi<graphs.size();gi++){
            graphs[gi]->Draw("SAME");
            //mxpts[gi]->Draw("SAME");
        }
        ci->Update();
        TGraph *gmx = new TGraph(5,tof.data(),tot_mx.data());
        gmx->SetMarkerStyle(20);
        gmx->SetMarkerColor(kRed);
        gmx->Draw("P SAME");
        ci->Update();

        fo1->cd();
        ci->Write();

        TH2D *h2 = (TH2D*)f->Get(Form("hEVsTof_L%d_cal",i));
        TFile *fcut2 = new TFile(cuts2[i].c_str(),"READ");
        TCanvas *c2i = new TCanvas(Form("hEVsTof_L%d_cal",i));

        std::vector<double> tof2 {};
        std::vector<TGraphErrors*> graphs2 {};
        std::vector<TGraphErrors*> mxpts2 {};

        std::vector<double> E_mx {};

        for (int j=1;j<=5;j++){
            TCutG *cut = (TCutG*)fcut2->Get(Form("cut_s%d",j));
            std::vector<std::vector<double>> data = data_load(h2,cut);
            std::vector<double> x_j = data[0];
            std::vector<double> y_j = data[1];
            std::vector<double> err_j = data[2];

            TGraphErrors *g_j = new TGraphErrors(x_j.size(),x_j.data(),y_j.data(),nullptr,err_j.data());
            g_j->SetName(Form("g2_s%d_L%d",j,i));
            std::vector<double> max_vals = max_find(g_j);
            double tof_j = max_vals[0];
            double E_j = max_vals[1];
            double e_j = max_vals[2];
            TGraphErrors *mxpt = new TGraphErrors(1,&tof_j,&E_j,nullptr,&e_j);
            mxpt->SetName(Form("mxpt2_s%d_L%d",j,i));
            mxpt->SetMarkerStyle(20);
            mxpt->SetMarkerColor(kRed);

            graphs2.push_back(g_j);
            //mxpts2.push_back(mxpt);

            tof2.push_back(tof_j);
            E_mx.push_back(E_j);
        }

        c2i->cd();
        h2->Draw("COLZ");
        for (size_t gi=0;gi<graphs2.size();gi++){
            graphs2[gi]->Draw("SAME");
            //mxpts2[gi]->Draw("SAME");
        }
        c2i->Update();

        TGraph *gE = new TGraph(5,tof2.data(),E_mx.data());
        gE->SetMarkerStyle(20);
        gE->SetMarkerColor(kRed);
        gE->Draw("P SAME");
        c2i->Update();

        TGraph *g = new TGraph(5,x[i].data(),y.data());
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlue);
        g->Draw("P SAME");
        c2i->Update();

        fo1->cd();
        c2i->Write();

    }
    f->Close();
    fo1->Close();

    return;
}