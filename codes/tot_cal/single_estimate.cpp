#include "../myConst.h"
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSpectrum.h>

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
    std::vector<double> counts {};
    
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
        f_1->SetParLimits(2,0.1,5.0);
        py_tmp->Fit(f_1,"QRB+");
        double max_val = f_1->GetParameter(0);
        double max_y_fit = f_1->GetParameter(1);
        double y_err = f_1->GetParError(1);
        //y_err=y_err/(max_val*num_ent);

        if (/* abs(y_err)<=25 && */ cut->IsInside(xval,max_y_fit)){
            x.push_back(xval);
            y.push_back(max_y_fit);
            err.push_back(y_err);
            counts.push_back(max_val);
        } 

    }

    return {x,y,err,counts};
}


std::vector<double> max_find(const TGraphErrors *g){
    int n = g->GetN();
    double x_max {0.0};
    double y_max {0.0};
    double err_max {0.0};

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
    
    return {x_max,y_max,err_max};
}


TH2D *cust_axis(TH2D *h, const TF1 *f){
    int nx = h->GetNbinsX();
    //int ny = h->GetNbinsY();
    int ny=310;
    std::vector<double> new_edges(ny+1);
    for (int j=0;j<=ny;j++){
        double yE = h->GetYaxis()->GetBinLowEdge(j+1);
        new_edges[j]=f->Eval(yE);
    }

    std::sort(new_edges.begin(),new_edges.end());

    TH2D *hnew = new TH2D("hnew",h->GetTitle(),nx,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),ny,new_edges.data());
    
    for (int i=1;i<=nx;i++) {
        for (int j=1;j<=ny;j++) {
            double x = h->GetXaxis()->GetBinCenter(i);
            double y = h->GetYaxis()->GetBinCenter(j);
            double z = h->GetBinContent(i, j);

            hnew->Fill(x, f->Eval(y), z);
        }
    }

    return hnew;
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

TH1D *ht_1D(const TH1D *py,const TCutG *cut,const double xval){
    TH1D *py_tmp = (TH1D*)py->Clone();
    int n = py->GetNbinsX();
    for (int i=1;i<=n;i++){
        double yval = py->GetBinCenter(i);
        if (!cut->IsInside(xval,yval)){
            py_tmp->SetBinContent(i,0);
            py_tmp->SetBinError(i,0);
        }
    }
    py_tmp->SetEntries(py_tmp->Integral(1,n));
    return py_tmp;
}

double avg_find(std::vector<double> v){
    double sum {0};
    for (size_t i=0;i<v.size();i++){
        sum+=v[i];
    }
    return sum/v.size();
}

void single_estimate(
    const std::string &file = "input.root", //file containing hTotVsTof_det histograms
    const std::string &tfile = "input.txt", // file containing mu values
    const std::string &ofile = "output",
    const std::string &cut_pattern = "./TCuts/l",
    bool He = false
){
    TFile *f = new TFile(file.c_str(),"READ");
    std::string oroot;
    oroot = ofile+".root";
    TFile *fo1 = new TFile(oroot.c_str(),"RECREATE");

    ofstream output;
    std::string output_file = ofile+".txt";
    output.open(output_file.c_str(),std::ios_base::trunc);
    std::vector<double> mu {};
    std::vector<double> mu_err {};
    std::vector<int> hbases {hbase1, hbase2, hbase3};    

    for (int i=0;i<3;i++){
        std::vector<double> mu_i=mu_read(tfile,hbases[i]);
        mu.push_back(mu_i[0]);
        mu_err.push_back(mu_i[1]);
    }

    std::vector<std::vector<double>> tofs {L0_tof,L1_tof,L2_tof};
    
    std::vector<double> pnch_pts {cos_mu, pro_pt, deu_pt, tri_pt}; //energies to be used for parameter estimation
    
    int num_sec=3;
    if (He){
        pnch_pts.push_back(hel_pt);
        pnch_pts.push_back(alp_pt);
        num_sec=5;
    }

    std::vector<std::vector<double>> tots {};
    std::vector<std::vector<double>> tots_errs {};
    for (int i=0;i<3;i++){
        std::vector<double> tot_i {};
        std::vector<double> tot_err_i {};

        tot_i.push_back(mu[i]); //add mu value as first point
        tot_err_i.push_back(mu_err[i]);

        TH2D *h = (TH2D*)f->Get(Form("hTotVsTof_det%d",hbases[i]));
        TFile *fcut = new TFile((cut_pattern+std::to_string(i)+".root").c_str(),"READ");

        for (int j=1;j<=num_sec;j++){
            double tof_j = tofs[i][j-1];
            int tof_j_bin = h->GetXaxis()->FindBin(tof_j);
            TH1D *py = (TH1D*)h->ProjectionY(Form("py_L%d_s%d",i,j),tof_j_bin-1,tof_j_bin+1);
            py->SetTitle(Form("Y Proj at ToF %.1f",tof_j));

            TCutG *cut = (TCutG*)fcut->Get(Form("cut_s%d",j));

            for (int pyi=1;pyi<=py->GetNbinsX();pyi++){
                double yval = py->GetBinCenter(pyi);
                if (!cut->IsInside(tof_j,yval)){
                    py->SetBinContent(pyi,0);
                    py->SetBinError(pyi,0);
                }
            }

            TSpectrum *s = new TSpectrum(2);
            Int_t nfound = s->Search(py,2,"",0.1);
            Double_t *xpeaks = s->GetPositionX();
            Double_t *ypeaks = s->GetPositionY();
            for (Int_t k=0;k<nfound;k++){
                double peak_x = xpeaks[k];
                double peak_y = ypeaks[k];
                TF1 *f_gaus = new TF1(Form("f_gaus_L%d_s%d_p%d",i,j,k),"gaus",peak_x-1,peak_x+1);
                f_gaus->SetParameters(peak_y, peak_x, 1.0);
                //f_gaus->FixParameter(1,peak_x);
                py->Fit(f_gaus,"BRQ+");
                if (nfound==1){
                    tot_i.push_back(f_gaus->GetParameter(1));
                    tot_err_i.push_back(f_gaus->GetParError(1));
                    /* tot_i.push_back(peak_x);
                    tot_err_i.push_back(f_gaus->GetParameter(2)); */
                } else {
                    if (k==1){
                        tot_i.push_back(f_gaus->GetParameter(1));
                        tot_err_i.push_back(f_gaus->GetParError(1));
                        /* tot_i.push_back(peak_x);
                        tot_err_i.push_back(f_gaus->GetParameter(2)); */
                    }
                }
                delete f_gaus;
            }


            delete s;

            fo1->cd();
            py->Write();

        }

        tots.push_back(tot_i);
        tots_errs.push_back(tot_err_i);

    }

    for (size_t i=0;i<tots.size();i++){
        std::cout<<"Detector "<<hbases[i]<<":"<<"\n";
        for (size_t j=0;j<tots[i].size();j++){
            std::cout<<"\t"<<tots[i][j]<<"\t"<<tots_errs[i][j];
        }
        std::cout<<"\n";
    }



    

    fo1->Close();

    return;
}