#include "../myConst.h"
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TSpectrum.h>



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
    std::vector<double> pnch_pts_err {0.1, pro_pt_err, deu_pt_err, tri_pt_err};
    int num_sec=3;
    if (He){
        pnch_pts.push_back(hel_pt);
        pnch_pts.push_back(alp_pt);
        pnch_pts_err.push_back(hel_pt_err);
        pnch_pts_err.push_back(alp_pt_err);
        num_sec=5;
    }

    std::vector<std::vector<double>> tots {};
    std::vector<std::vector<double>> tots_errs {};

    //find the tot peaks from the y projections at certain tofs
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

            TSpectrum *s = new TSpectrum(3);
            Int_t nfound = s->Search(py,2,"",0.3);
            Double_t *xpeaks = s->GetPositionX();
            
            std::vector<double> peak_vals {};
            for (int k=0;k<nfound;k++){
                peak_vals.push_back(xpeaks[k]);
            }

            std::sort(peak_vals.begin(),peak_vals.end());

            for (size_t k=0;k<peak_vals.size();k++){
                double xpeak = peak_vals[k];
                int xbin = py->GetXaxis()->FindBin(xpeak);
                double xcenter = py->GetXaxis()->GetBinCenter(xbin);
                double A0 = py->GetBinContent(xbin);
                TF1 *f1 = new TF1("f1","gaus",xcenter-0.6,xcenter+0.6);
                f1->SetParameters(A0,xcenter,0.1);
                f1->SetParLimits(1,xcenter-0.6,xcenter+0.6);
                py->Fit(f1,"RQ+");
                if (nfound==1 && k==0){
                    double tot_val = f1->GetParameter(1);
                    double tot_err = f1->GetParError(1);
                    tot_i.push_back(tot_val);
                    tot_err_i.push_back(tot_err);

                } else if (nfound>1 && k==nfound-1){
                    double tot_val = f1->GetParameter(1);
                    double tot_err = f1->GetParError(1);
                    tot_i.push_back(tot_val);
                    tot_err_i.push_back(tot_err);

                }
                delete f1;
            }

            delete s;

            fo1->cd();
            py->Write();

        }

        tots.push_back(tot_i);
        tots_errs.push_back(tot_err_i);

    }

    /* for (size_t i=0;i<tots.size();i++){
        std::cout<<"Detector "<<hbases[i]<<":"<<"\n";
        for (size_t j=0;j<tots[i].size();j++){
            std::cout<<"\t"<<tots[i][j]<<"\t"<<tots_errs[i][j];
        }
        std::cout<<"\n";
    }
    */

    for (int i=0;i<3;i++){
        TCanvas *c = new TCanvas(Form("c_E_vs_ToT_det%d",hbases[i]),Form("Det %d",hbases[i]));
        TGraphErrors *g = new TGraphErrors(tots[i].size(),&tots[i][0],&pnch_pts[0],&tots_errs[i][0],&pnch_pts_err[0]);
        g->SetTitle(Form("Energy vs ToT Det %d",hbases[i]));
        g->GetXaxis()->SetTitle("ToT (ns)");
        g->GetXaxis()->SetLimits(0,40);
        g->GetYaxis()->SetTitle("Energy (MeV)");
        g->GetYaxis()->SetRangeUser(0,200);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.5);
        c->cd();
        g->Draw("APE");

        //recall amp=A*exp(k*tot)+B
        //E=C*amp+D
        //assuming C=1 and D=0 for simplicity

        TF1 *f1 = new TF1("f1", "[0]*exp([1]*x)+[2]",0,40);
        f1->SetParameters(1,0.1,0);
        f1->SetParLimits(0,0,100);
        f1->SetParLimits(1,0,1);
        f1->SetParLimits(2,0,100);
        g->Fit(f1,"RQ");

        fo1->cd();
        c->Write();

        output<<std::setw(10)<<f1->GetParameter(0)<<std::setw(5)<<" +/- "<<std::setw(10)<<f1->GetParError(0)<<"\t"
              <<std::setw(10)<<f1->GetParameter(1)<<std::setw(5)<<" +/- "<<std::setw(10)<<f1->GetParError(1)<<"\t"
              <<std::setw(10)<<f1->GetParameter(2)<<std::setw(5)<<" +/- "<<std::setw(10)<<f1->GetParError(2)<<"\n";
    }



    output.close();
    f->Close();
    fo1->Close();

    return;
}