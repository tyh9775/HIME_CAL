#include "../myConst.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCutG.h>


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

    int xmin_i = h->GetXaxis()->FindFixBin(cut_range(cut)[0]);
    int xmax_i = h->GetXaxis()->FindFixBin(cut_range(cut)[1]);
    int ymin_i = h->GetYaxis()->FindFixBin(cut_range(cut)[2]);
    int ymax_i = h->GetYaxis()->FindFixBin(cut_range(cut)[3]);

    int nybins = h->GetNbinsY();
    for (int xj=xmin_i;xj<xmax_i;xj++){
        TH1D *py_tmp = htmp->ProjectionY(Form("py%d",xj),xj,xj);
        double xval=h->GetXaxis()->GetBinCenter(xj);
        for (int yj=1;yj<=nybins;yj++){
            double yval=h->GetYaxis()->GetBinCenter(yj);
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

        TF1 *f_1 = new TF1("f_1","gaus",py_tmp->GetBinCenter(y_low),py_tmp->GetBinCenter(y_high));
        f_1->SetParameters(max_count,max_y,1);
        py_tmp->Fit(f_1,"QRB+");
        double max_y_fit = f_1->GetParameter(1);
        double y_err = f_1->GetParError(1);

        if (abs(y_err)<=2.5 && abs(y_err)>0.0 && cut->IsInside(xval,max_y_fit)){
            x.push_back(xval);
            y.push_back(max_y_fit);
            err.push_back(y_err);
            
        } 

    }

    return {x,y,err};
}

void muon_load(
    const std::string &input = "input.root",
    const std::string &output = "output"
){
    TFile *f = new TFile(input.c_str(),"READ");
    TFile *fo1 = new TFile((output+".root").c_str(),"RECREATE");

    std::ofstream muon_output;
    muon_output.open((output+".txt").c_str(),std::ios_base::trunc);
    muon_output<<std::setw(8) << "Det Num " << std::setw(12) << "Avg. ToT" <<std::setw(12)<<"Error"<<"\n";

    double xmin {-20.0};
    double xmax {20.0};
    double ymin {16.0};
    double ymax {28.0};
    TCutG *cut = new TCutG("cut",5);
    cut->SetPoint(0,xmin,ymin);
    cut->SetPoint(1,xmax,ymin);
    cut->SetPoint(2,xmax,ymax);
    cut->SetPoint(3,xmin,ymax);
    cut->SetPoint(4,xmin,ymin);

    for (int i=0;i<72;i++){
        TCanvas *ci = new TCanvas(Form("c%d",i),Form("Det %d",i));
        TH2D *h = (TH2D*)f->Get(Form("hTotVsTDiff_mod_%d",i))->Clone();
        TH2D *hc = (TH2D*)h->Clone(Form("hTotVsTDiff_mod_%d_cut",i));
        for (int x=1;x<=hc->GetNbinsX();x++){
            for (int y=1;y<=hc->GetNbinsY();y++){
                double xx = hc->GetXaxis()->GetBinCenter(x);
                double yy = hc->GetYaxis()->GetBinCenter(y);
                if (!cut->IsInside(xx,yy)){
                    hc->SetBinContent(x,y,0);
                    hc->SetBinError(x,y,0);
                }
            }
        }
        hc->SetEntries(hc->Integral(1,hc->GetNbinsX(),1,hc->GetNbinsY()));

        std::vector<std::vector<double>> data = data_load(hc,cut);
        int n = data[0].size();
        TGraphErrors *g = new TGraphErrors(n,data[0].data(),data[1].data(),nullptr,data[2].data());
        g->SetName(Form("Det %d mu points",i));
        
        ci->cd();
        hc->Draw();
        g->Draw("SAME");
    
        double avg {0.0};
        double avg_err {0.0};
        double y_wgt {0.0};
        double w_sum {0.0};

        for (int xi=0;xi<n;xi++){
            if (data[2][xi] > 0){
                y_wgt += data[1][xi]/pow(data[2][xi],2);
                w_sum += 1.0/pow(data[2][xi],2);
            }
        }

        if (w_sum>0.0){
            avg = y_wgt/w_sum;
            avg_err = sqrt(1.0/w_sum);
        }

        muon_output<<std::setw(8) << i << std::setw(12) << avg << std::setw(8) << "+/-" << std::setw(12) << avg_err << "\n";

        fo1->cd();
        ci->Write();
    }

    fo1->Close();
    muon_output.close();
    
    return;
}