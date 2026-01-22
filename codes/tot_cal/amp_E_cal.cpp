#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>

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
        f_1->SetParLimits(2,0.1,100.0);
        py_tmp->Fit(f_1,"QRB+");
        double max_y_fit = f_1->GetParameter(1);
        double y_err = f_1->GetParError(1);

        if (cut->IsInside(xval,max_y_fit)){
            x.push_back(xval);
            y.push_back(max_y_fit);
            if (abs(y_err)<=500.0){   
                err.push_back(y_err);
            } else {
                err.push_back(500.0);
            }
        } 

    }

    return {x,y,err};
}


std::vector<double> max_find(const TGraphErrors *g){
    int n = g->GetN();
    double x_max {0};
    double y_max {0};
    double err_max {0};
    double x0 = g->GetPointX(0);
    double xf = g->GetPointX(n-1);
    double y0 = g->GetPointY(0);
    double yf = g->GetPointY(n-1);
    for (int i=0;i<n;i++){
        double xi = g->GetPointX(i);
        double yi = g->GetPointY(i);
        double err_i = g->GetErrorY(i);
        if (err_i>100.0 || err_i<0.0){
            continue;
        }
        if (yi>=y_max){
            y_max=yi;
            x_max=xi;
            err_max=err_i;
        }
    }
    
    return {x_max,y_max,err_max,x0,y0,xf,yf};
}

double avg_find(std::vector<double> v){
    double sum {0};
    for (size_t i=0;i<v.size();i++){
        sum+=v[i];
    }
    return sum/v.size();
}



double f_exp(double ToT,const std::vector<double> &parA){
    return parA[0]*exp(parA[1]*ToT)+parA[2];
}

double f_der(double ToT,const std::vector<double> &parA){
    return parA[0]*parA[1]*exp(parA[1]*ToT);
}

std::vector<double> mu_read(const std::string &tfile,int det){
    std::ifstream infile(tfile);
    std::string line;
    int di {0};
    double tot {0.0};
    double err {0.0};

    std::getline(infile,line);

    while (std::getline(infile,line)){
        if (line.empty()) {
            continue;
        }
        
        std::istringstream ss(line);
        std::string pm;
        ss>>di>>tot>>pm>>err;

        if (di==det){
            return {tot,err};
        }
    }
    return {};
}



void amp_E_cal(
    const std::string &file = "input",
    const std::string &ofile = "output",
    const std::string &cut_pattern = "TCuts/cal_l",
    bool He = false
){
    TFile *f = new TFile((file+".root").c_str(),"READ");
    std::string oroot;
    oroot = ofile+".root";
    TFile *fo1 = new TFile(oroot.c_str(),"RECREATE");

    ofstream output;
    std::string output_file = ofile+".txt";
    output.open(output_file.c_str(),std::ios_base::trunc);

    std::vector<double> mu_amp {};
    std::vector<double> mu_err {};
    std::vector<int> hbases {hbase1, hbase2, hbase3};    
    std::vector<std::vector<double>> tof_ls {L0_tof,L1_tof,L2_tof};

    for (int i=0;i<3;i++){
        std::vector<double> mu_i=mu_read(tfile,hbases[i]);
        //std::cout<<mu_i<<"\n";
        mu_amp.push_back(f_exp(mu_i[0],parA_mat[i]));
        mu_err.push_back(f_der(mu_i[0],parA_mat[i])*mu_i[1]);
    }

    std::vector<double> pnch_pts {cos_mu,pro_pt, deu_pt, tri_pt};
    if (He){
        pnch_pts.push_back(hel_pt);
        pnch_pts.push_back(alp_pt);
    }
    int num_sec = pnch_pts.size()-1;

    if (E_min>=0.0){
        pnch_pts.insert(pnch_pts.begin(),E_min);
    }

    for (int i=0;i<3;i++){
        TH2D *h = (TH2D*)f->Get(Form("hAmpVsTof_L%d_cal",i));

        TFile *fcut = new TFile(cuts[i].c_str(),"READ");
        std::vector<double> amps {};
        std::vector<double> amps2 {};
        std::vector<double> amps_err {};
        std::vector<double> amps2_err {};

        std::vector<TGraphErrors*> graphs {};
        std::vector<TGraph*> mxpts {};
        std::vector<TGraphErrors*> mxpts2 {};
        for (int j=1;j<=num_sec;j++){
            TCutG *cut = (TCutG*)fcut->Get(Form("cut_s%d",j));
            std::vector<std::vector<double>> data = data_load(h,cut);
            std::vector<double> x_j {data[0]};
            std::vector<double> y_j {data[1]};
            std::vector<double> err_j {data[2]};

            TGraphErrors *g_j = new TGraphErrors(x_j.size(),x_j.data(),y_j.data(),nullptr,err_j.data());
            g_j->SetName(Form("g_s%d_L%d",j,i));
            double tof_j = max_find(g_j)[0];
            double amp_j = max_find(g_j)[1];
            //double amp_err_j = max_find(g_j)[2];
            double amp_err_j = avg_find(err_j);
            amps.push_back(amp_j);
            amps_err.push_back(amp_err_j);
            TGraphErrors *gc = (TGraphErrors*)g_j->Clone();
            graphs.push_back(gc);
            TGraph *mxpt = new TGraph(1,&tof_j,&amp_j);
            mxpt->SetName(Form("mxpt_s%d_L%d",j+1,i));
            mxpt->SetMarkerStyle(20);
            mxpt->SetMarkerColor(kRed);
            mxpts.push_back(mxpt);

            TH2D *h2 = htmp_make(h,cut);
            double xi = tof_ls[i][j-1];
            int xbi = h2->GetXaxis()->FindBin(xi);
            TH1D *py_tmp = h2->ProjectionY(Form("py2_s%d_L%d",j,i),xbi-1,xbi+1);
            int max_ybin = py_tmp->GetMaximumBin();
            double max_y = py_tmp->GetBinCenter(max_ybin);
            double max_count = py_tmp->GetBinContent(max_ybin);

            int y_low = py_tmp->FindFirstBinAbove(0);
            int y_high = py_tmp->FindLastBinAbove(0);
            /* double y_l = cut_range(cut)[2];
            double y_h = cut_range(cut)[3]; */
            double y_l = py_tmp->GetBinCenter(y_low);
            double y_h = py_tmp->GetBinCenter(y_high);
            //TF1 *f_1 = new TF1("f_1","gaus",py_tmp->GetBinCenter(y_low),py_tmp->GetBinCenter(y_high));
            TF1 *f_1 = new TF1("f_1","gaus",y_l,y_h);
            f_1->SetParameters(max_count,max_y,1);
            f_1->SetParLimits(1,y_l,y_h);
            f_1->SetParLimits(2,0.1,100.0);
            py_tmp->Fit(f_1,"QRB");
            double max_y_fit = f_1->GetParameter(1);
            double y_err = f_1->GetParError(1);
            amps2.push_back(max_y_fit);
            amps2_err.push_back(y_err);
            TGraphErrors *mxpt2 = new TGraphErrors(1,&xi,&max_y_fit,nullptr,&y_err);
            mxpt2->SetName(Form("mxpt2_s%d_L%d",j,i));
            mxpt2->SetMarkerStyle(20);
            mxpt2->SetMarkerColor(kBlue);
            mxpts2.push_back(mxpt2);    

        }

        
        TH1D *py = (TH1D*)h->ProjectionY("py");
        int y_min_bin = py->FindFirstBinAbove(0);
        double amp_min = py->GetXaxis()->GetBinCenter(y_min_bin);
        double amp_min_err = py->GetXaxis()->GetBinWidth(y_min_bin);


        amps.insert(amps.begin(),mu_amp[i]);
        amps2.insert(amps2.begin(),mu_amp[i]);
        amps_err.insert(amps_err.begin(),mu_err[i]);
        amps2_err.insert(amps2_err.begin(),mu_err[i]);

        if (E_min>=0.0){
            amps.insert(amps.begin(),amp_min);
            amps2.insert(amps2.begin(),amp_min);
            amps_err.insert(amps_err.begin(),amp_min_err);
            amps2_err.insert(amps2_err.begin(),amp_min_err);
        }

        TCanvas *c0 = new TCanvas(Form("hAmpVsTof_cal_L%d",i));
        c0->cd();
        h->Draw(); 

        for (size_t g=0; g< graphs.size(); g++) {
            graphs[g]->Draw("SAME");
            mxpts[g]->Draw("P SAME");
            mxpts2[g]->Draw("P SAME");
        }

        fo1->cd();
        c0->Write();

        std::cout<<"Layer "<<i<<"\n";
        for (size_t j=0;j<pnch_pts.size();j++){
            std::cout << "Amp: " << amps[j] << "\n";
            std::cout << "E: " << pnch_pts[j] << "\n";
        }

        TGraphErrors *g_i = new TGraphErrors(amps.size(),amps.data(),pnch_pts.data(),amps_err.data(),nullptr);
        g_i->SetName(Form("g_L%d",i));
        g_i->SetMarkerStyle(20);
        TF1 *f_i = new TF1(Form("f_L%d",i),"pol1", g_i->GetXaxis()->GetXmin(), g_i->GetXaxis()->GetXmax());
        f_i->SetParameters(0.0,1.0);
        if (limit){
            f_i->FixParameter(0,0.0);
        }
        g_i->Fit(f_i,"QRB+");


        TGraphErrors *g2_i = new TGraphErrors(amps2.size(),amps2.data(),pnch_pts.data(),amps2_err.data(),nullptr);
        g2_i->SetName(Form("g2_L%d",i));
        g2_i->SetMarkerStyle(20);
        TF1 *f2_i = new TF1(Form("f2_L%d",i),"pol1", g2_i->GetXaxis()->GetXmin(), g2_i->GetXaxis()->GetXmax());
        f2_i->SetParameters(0.0,1.0);
        if (limit){
            f2_i->FixParameter(0,0.0);
        }
        g2_i->Fit(f2_i,"QRB+");

        fo1->cd();
        g_i->Write();
        g2_i->Write();

        output<<"Layer "<<i<<"\n";
        output<<"y-int: "<<f_i->GetParameter(0)<<"+/-"<<f_i->GetParError(0)<<", slope: "<<f_i->GetParameter(1)<<"+/-"<<f_i->GetParError(1)<<"\n";
        output<<"\n";
        output<<"v2 "<<"\n";
        output<<"y-int: "<<f2_i->GetParameter(0)<<"+/-"<<f2_i->GetParError(0)<<", slope: "<<f2_i->GetParameter(1)<<"+/-"<<f2_i->GetParError(1)<<"\n";
        output<<"\n";
        fcut->Close();
    }

    f->Close();
    output.close();
    fo1->Close();

    return;
}