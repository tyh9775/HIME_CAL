#include "../myConst.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <TList.h>
#include <TCutG.h>
#include <TKey.h>
#include <Math/Minimizer.h>



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

double y_pt_find(const TH2D *h,const TCutG *cut){
    TH2D *htmp = htmp_make(h,cut);
    int xfix = htmp->GetXaxis()->FindBin(15.0);
    TH1D *pytmp0 = htmp->ProjectionY("pytmp0",xfix,xfix);
    
    int max_ybin = pytmp0->GetMaximumBin();
    double max_y = pytmp0->GetBinCenter(max_ybin);
    double max_count = pytmp0->GetBinContent(max_ybin);

    int y_low = pytmp0->FindFirstBinAbove(0);
    int y_high = pytmp0->FindLastBinAbove(0);

    TF1 *f_1 = new TF1("f_1","gaus",pytmp0->GetBinCenter(y_low),pytmp0->GetBinCenter(y_high));
    f_1->SetParameters(max_count,max_y,1);
    pytmp0->Fit(f_1,"QRB");
    double y0 = f_1->GetParameter(1);

    return y0;
}


TCutG *cut_adjust(const TCutG *c, double ydiff){
    TCutG *c1 = (TCutG*)c->Clone();
    if (ydiff==0.0){
        return c1;
    }
    int n = c->GetN();
    c1->Set(n);
    for (int i=0;i<n;i++){
        double x = c->GetPointX(i);
        double y = c->GetPointY(i)-ydiff;

        c1->SetPoint(i,x,y);
    }

    return c1;    
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

//need to adjust the x ranges (use the walk corrected ones to be more consistent)
void mult_loader(
    const std::string &file = "input.root",
    const std::string &ofile = "output.root",
    const std::string &cut_pattern = "./TCuts/cal_l",
    bool cal=false,
    bool dyn_cut=true
){
    //open file
    auto *f = new TFile(file.c_str(),"READ");
    //create output file
    auto *fo1 = new TFile(ofile.c_str(),"RECREATE");

        
    gROOT->SetBatch(kTRUE);

    std::vector<int> hstarts {hstart1,hstart2,hstart3};
    std::vector<int> hstops {hstop1,hstop2,hstop3};
    std::vector<int> hbases {hbase1,hbase2,hbase3};

    /* TList *list[72];

    for (int i=0;i<72;i++){
        list[i] = new TList(Form("det%d",i));
    } */

    TH2D *h_all = new TH2D("hTotVsTof_all","",2000, -100, 300, 150, 0, 75);
    for (int L=0;L<1;L++){
        //Graphical cuts
        auto *fcut = new TFile((cut_pattern+to_string(L)+".root").c_str(),"READ");
        
        int num_sec {0};
        
        TKey *key;
        TIter next_key(fcut->GetListOfKeys());
        while ((key = (TKey*)next_key())){
            num_sec+=1;
        }


        TH2D *hL = new TH2D(Form("hTotVsTof_L%d",L), "", 2000, -100, 300, 150, 0, 75);
        if (cal){
            hL->SetName(Form("hTotVsTof_L%d_cal",L));
        }

        TH2D *hLc = new TH2D(Form("hTotVsTof_L%d_cut",L), "", 2000, -100, 300, 150, 0, 75);
        if (cal){
            hLc->SetName(Form("hTotVsTof_L%d_cut_cal",L));
        }
        

        double y0 {0.0};
        
        TH2D *hcopy0 = nullptr;
        if (!cal){
            hcopy0 = (TH2D*)f->Get(Form("hTotVsTofBar%d",hbases[L]));
            //find the offset for the cuts of the detectors by using the first section
            TCutG *cut0 = (TCutG*)fcut->Get("cut_s0");
            y0 = y_pt_find(hcopy0,cut0);         
        } else {
            hcopy0 = (TH2D*)f->Get(Form("hTotVsTofBar_sh_t%d",hbases[L]));
            hcopy0->RebinY(5);
        }

        std::vector<int> degrees {};

        for (int j=0;j<num_sec;j++){
            /* if (std::find(sec_omit.begin(),sec_omit.end(),j+1) != sec_omit.end()){
                degrees.push_back(0);
                continue;
            } */
            TCutG *cutj = (TCutG*)fcut->Get(Form("cut_s%d",j+1));
            std::cout<<j<<"\n";
            std::vector<std::vector<double>> data0 = data_load(hcopy0,cutj);

            double bic_start {1e6};
            int n_j = data0[0].size();            
            TGraphErrors* gj = new TGraphErrors(n_j,data0[0].data(),data0[1].data(),nullptr,data0[2].data());
            gj->SetName(Form("g_s%d_dbase",j+1));
            
            int best_deg {0};
            TF1 *best_fit = nullptr;
            int dstrt {min_deg};
            int dstp {max_deg};
            if (j==1 || j==2){
                dstrt=2;
            } else if (j>=5){
                dstrt=1;
                dstp=1;
            } 
            for (int deg=dstrt;deg<=dstp;deg++){
                TF1 *fit_f = new TF1(Form("fit_deg_%d",deg),Form("pol%d",deg),data0[0][0],data0[0][n_j-1]);
                //fit_f->SetParLimits(0,-10,10);
                gj->Fit(fit_f,"QR");
                double chi2 {fit_f->GetChisquare()};
                int k {deg+1};
                double bic {k*chi2*log(n_j)+chi2};
                if (bic<bic_start){
                    bic_start=bic;
                    best_deg=deg;
                    if (best_fit){
                        delete best_fit;
                    }
                    fit_f->SetName(Form("f_s%d_dbase",j+1));
                    best_fit = (TF1*)fit_f->Clone();
                }

                if (fit_f != best_fit){
                    delete fit_f;
                }          
            }

            degrees.push_back(best_deg);
        }




        for (int i=hstarts[L];i<=hstops[L];i++){
            TH2D *h = nullptr;

            if (!cal){
                h = (TH2D*)f->Get(Form("hTotVsTofBar%d",i));
            } else {
                if (i==hbases[L]){
                    h = (TH2D*)f->Get(Form("hTotVsTofBar_sh_t%d", i))->Clone();
                } else {
                    h = (TH2D*)f->Get(Form("hTotVsTofBar_sh_t%d", i))->Clone();
                    h->RebinY(5);
                }

            }
        
            TH2D *hcopy = (TH2D*)h->Clone();
            if (!cal){
                hcopy->SetName(Form("hcopy%d",i));
                hcopy->SetTitle(Form("Tot vs Tof Det %d",i));
            } else {
                hcopy->SetName(Form("hcopy%d_cal",i));
                hcopy->SetTitle(Form("Tot vs Tof Det %d (cal)",i));
            }            

            hL->Add(hcopy);

            int nxbins = hcopy->GetNbinsX();
            int nybins = hcopy->GetNbinsY();
            

            TH2D *hcut = new TH2D(Form("hcut%d",i), "", 2000, -100, 300, 150, 0, 75);

            std::cout<<"\n";
            std::cout<<"Det "<<i<<"\n";

            std::vector<TGraphErrors*> graphs {};

            TCutG *cuti =(TCutG*)fcut->Get("cut_s1");

            double yi {0.0};
            double ydiff_i {0.0};
            
            if (!cal && dyn_cut){
                yi = y_pt_find(hcopy,cuti);
                ydiff_i = y0-yi;
            }

            for (int j=0;j<num_sec;j++){
                /* if (std::find(sec_omit.begin(),sec_omit.end(),j+1) != sec_omit.end()){
                    continue;
                }
 */
                TCutG *cut_ini = (TCutG*)fcut->Get(("cut_s"+to_string(j+1)).c_str());
                TCutG *cut = (TCutG*)cut_ini->Clone();

                cut = cut_adjust(cut_ini,ydiff_i);

                TH2D *htmp = htmp_make(hcopy,cut);

                hcut->Add(htmp);


                std::vector<std::vector<double>> data = data_load(hcopy,cut);
            

                /* if (std::find(sec_omit.begin(),sec_omit.end(),j+1) != sec_omit.end()){
                    continue;
                } */
                int n_j = data[0].size();  
                TGraphErrors* gj = new TGraphErrors(n_j,data[0].data(),data[1].data(),nullptr,data[2].data());
                gj->SetName(Form("g_s%d_d%d",j,i));

                TF1 *fit_f = new TF1(Form("fit_deg_%d",degrees[j]),Form("pol%d",degrees[j]),data[0][0],data[0][n_j-1]);

                std::cout<<"section "<<j<<" best fit"<<"\n";
                gj->Fit(fit_f,"R");

                std::cout<<"****************************************"<<"\n";
                graphs.push_back(gj);

            }

            hLc->Add(hcut);

            TCanvas *hc1 = new TCanvas(Form("hcopy%d",i),Form("hTotVsTof Det %d (log z)",i));
            hc1->cd();
            gPad->SetLogz();
            hcopy->Draw();

            TCanvas *hc2 = new TCanvas(Form("hcut%d",i),Form("hTotVsTof Det %d (cut, log z)",i));
            hc2->cd();
            gPad->SetLogz();
            hcut->Draw();
            for (size_t j=0;j<graphs.size();j++){
                hc1->cd();
                graphs[j]->Draw("SAME");
                hc2->cd();
                graphs[j]->Draw("SAME");
            }
            
            /* list[i]->Add(hcopy);
            list[i]->Add(hcut); */

            fo1->cd();
            hc1->Write();
            hc2->Write();

        }

        fo1->cd();
        hL->Write();
        hLc->Write();
        h_all->Add(hL);
        
        fcut->Close();
        
    
    }
    fo1->cd();
    h_all->Write();
    fo1->Close();
    f->Close();    
    gROOT->SetBatch(kFALSE);

    //abort();
    return;
}
