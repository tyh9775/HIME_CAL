#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TF2.h>
#include <TMatrixD.h>
#include <TStyle.h>
#include <TObjArray.h>
#include <TObjString.h>


double func(const double x, const double y, const double a, const double b, const double c){
    return a*x+b*y+c;
}


void h_chunk_mod(const TH2D *h0,TH2D *h,const double scl,const int x_l,const int x_u,const int y_l,const int y_u){
    for (int xi=x_l;xi<x_u;xi++){
        for (int yi=y_l;yi<y_u;yi++){
            double count=h0->GetBinContent(xi,yi);
            if (count==0){
                continue;
            }
            double err=h0->GetBinError(xi,yi);
            double new_cnt=count*scl;
            if (new_cnt<0){
                h->SetBinContent(xi,yi,0.0);
            } else {
                h->SetBinContent(xi,yi,new_cnt);
            }            
            h->SetBinError(xi,yi,err*scl);
        }
    }    
    return;
}

void grid_make_v2(
    const std::string &input = "input.root",    
    const std::string &output = "output",
    const std::string per_txt = "percentages.txt",
    const std::vector<int> &num_ent = {1,1},
    const std::string &hname = "hHitPattern"
){
    TFile *f = TFile::Open(input.c_str(),"READ");
    TFile *fo1 = new TFile((output+".root").c_str(),"RECREATE");
    
    std::ofstream otxt;
    otxt.open((output+".txt").c_str(),std::ios_base::trunc);
    
    std::ofstream otxt2;
    otxt2.open((output+"_mult_est.txt").c_str(),std::ios_base::trunc);

    int tot_ent=0;
    for (size_t i=0;i<num_ent.size();i++){
        tot_ent += num_ent[i];
    }

    std::cout<<"Total entries: "<<tot_ent<<"\n";

    std::vector<std::vector<double>> sb_posx_mat {{sb2_xmin,sb2_xmax},{sb3_xmin,sb3_xmax},{sb4_xmin,sb4_xmax}};
    std::vector<std::vector<double>> sb_posy_mat {{sb2_ymin,sb2_ymax},{sb3_ymin,sb3_ymax},{sb4_ymin,sb4_ymax}};

    std::ifstream per_f(per_txt.c_str());

    std::vector<std::vector<double>> per_mat {};
    std::vector<std::vector<double>> err_mat {};
    std::vector<double> per_vec_tmp {};
    std::vector<double> per_err_tmp {};

    std::string line;
    while (std::getline(per_f,line)){
        size_t pos = line.find("+/-");
        double val = std::stod(line.substr(0,pos));
        double err = std::stod(line.substr(pos+3));
        per_vec_tmp.push_back(val);
        per_err_tmp.push_back(err);
        if (per_vec_tmp.size()==3){
            per_mat.push_back(per_vec_tmp);
            per_vec_tmp.clear();
            err_mat.push_back(per_err_tmp);
            per_err_tmp.clear();
        }
    }

    TMatrixD P(3,3);
    TMatrixD P_err(3,3);
    for (int i=0;i<3;i++){
        for (int j=0;j<3;j++){
            P(j,i) = per_mat[i][j];
            P_err(j,i) = err_mat[i][j];
        }
    }

    double x_c_sb2 = (sb2_xmax+sb2_xmin)/2.0;
    double x_c_sb3 = (sb3_xmax+sb3_xmin)/2.0;
    double x_c_sb4 = (sb4_xmax+sb4_xmin)/2.0;

    double x_e_sb2 = (sb2_xmax-sb2_xmin)/2.0;
    double x_e_sb3 = (sb3_xmax-sb3_xmin)/2.0;
    double x_e_sb4 = (sb4_xmax-sb4_xmin)/2.0;

    double y_c_sb2 = (sb2_ymax+sb2_ymin)/2.0;
    double y_c_sb3 = (sb3_ymax+sb3_ymin)/2.0;
    double y_c_sb4 = (sb4_ymax+sb4_ymin)/2.0;

    double y_e_sb2 = (sb2_ymax-sb2_ymin)/2.0;
    double y_e_sb3 = (sb3_ymax-sb3_ymin)/2.0;
    double y_e_sb4 = (sb4_ymax-sb4_ymin)/2.0;

    std::vector<double> x {x_c_sb2,x_c_sb3,x_c_sb4};
    std::vector<double> y {y_c_sb2,y_c_sb3,y_c_sb4};
    std::vector<double> x_err {x_e_sb2,x_e_sb3,x_e_sb4};
    std::vector<double> y_err {y_e_sb2,y_e_sb3,y_e_sb4};

    TMatrixD A(3,3);
    A(0,0)=x[0]; A(0,1)=y[0]; A(0,2)=1.0;
    A(1,0)=x[1]; A(1,1)=y[1]; A(1,2)=1.0;
    A(2,0)=x[2]; A(2,1)=y[2]; A(2,2)=1.0;


    TH2D *h_all = (TH2D*)f->Get(hname.c_str())->Clone();

    fo1->cd();
    h_all->Write();

    TH2D *h_add = (TH2D*)h_all->Clone("hHitPattern_add");
    h_add->Reset();
    TH2D *h_add_fine = (TH2D*)h_add->Clone("hHitPattern_add_fine");
    TH2D *h_add_sb = (TH2D*)h_add->Clone("hHitPattern_add_sb");
    TH2D *h_add_coarse = (TH2D*)h_add->Clone("hHitPattern_add_coarse");

    otxt2<<std::setw(15)<<" "<<std::setw(12)<<" Original "<<std::setw(12)<<" Fine mod "<<std::setw(12)<<" Coarse mod "<<"\n";


    for (int i=0;i<3;i++){
        std::cout<<"Layer "<<i<<"\n";
        TMatrixD p = P.GetSub(i,i,0,2);
        TMatrixD p_err = P_err.GetSub(i,i,0,2);

        //iterate to see effect of x and y uncertainties in the final fit
        //initial guess, no effect
        double a=1.0;
        double b=1.0;
        double c=0.0;
        for (int j=0;j<5;j++){
            TMatrixD W(3,3);
            TMatrixD Z(3,1);
            for (int k=0;k<3;k++){
                double sig2 = p_err(0,k)*p_err(0,k) + a*a*x_err[k]*x_err[k] + b*b*y_err[k]*y_err[k];
                W(k,k) = 1.0/sig2;
                Z(k,0) = p(0,k);
            }
            TMatrixD AT(TMatrixD::kTransposed, A);
            TMatrixD ATWA = AT*W*A;
            TMatrixD ATWZ = AT*W*Z;
            TMatrixD params = ATWA.Invert()*ATWZ;
            a = params(0,0);
            b = params(1,0);
            c = params(2,0);
            std::cout<<"iter "<<j<<" a: "<<a<<" b: "<<b<<" c: "<<c<<"\n";
        }
        otxt<<a<<", "<<b<<", "<<c<<"\n";

        std::string hname_layer = hname+"Layer"+std::to_string(i);
        if (hname=="hHitPattern_rnd"){
            hname_layer = "hHitPatternLayer"+std::to_string(i)+"_rnd";
        }

        TH2D *h = (TH2D*)f->Get(hname_layer.c_str())->Clone();
        TH2D *hcopy = (TH2D*)h->Clone();
        hcopy->Reset();
        hcopy->SetName(Form("hHitPatternLayer%d_modified_fine",i));
        int nbinsx = hcopy->GetNbinsX();
        int nbinsy = hcopy->GetNbinsY();
        for (int xi=1;xi<=nbinsx;xi++){
            for (int yi=1;yi<=nbinsy;yi++){
                double count = h->GetBinContent(xi,yi);
                if (count==0.0){
                    continue;
                }
                double err = h->GetBinError(xi,yi);
                double xval=hcopy->GetXaxis()->GetBinCenter(xi);
                double yval=hcopy->GetYaxis()->GetBinCenter(yi);
                double sclr=func(xval,yval,a,b,c);
                double new_cnt=count*sclr;
                if (new_cnt<0.0){
                    hcopy->SetBinContent(xi,yi,0.0);
                    hcopy->SetBinError(xi,yi,err*sclr);
                } else {
                    hcopy->SetBinContent(xi,yi,count*sclr);
                    hcopy->SetBinError(xi,yi,err*sclr);
                }

            }
        }


        TH2D *hcopy2 = (TH2D*)h->Clone();
        hcopy2->Reset();
        hcopy2->SetName(Form("hHitPatternLayer%d_sbs",i));
        int xi_sb2_l = h->GetXaxis()->FindBin(sb2_xmin);
        int xi_sb2_u = h->GetXaxis()->FindBin(sb2_xmax);
        int yi_sb2_l = h->GetYaxis()->FindBin(sb2_ymin);
        int yi_sb2_u = h->GetYaxis()->FindBin(sb2_ymax);

        int xi_sb3_l = h->GetXaxis()->FindBin(sb3_xmin);
        int xi_sb3_u = h->GetXaxis()->FindBin(sb3_xmax);
        int yi_sb3_l = h->GetYaxis()->FindBin(sb3_ymin);
        int yi_sb3_u = h->GetYaxis()->FindBin(sb3_ymax);


        int xi_sb4_l = h->GetXaxis()->FindBin(sb4_xmin);
        int xi_sb4_u = h->GetXaxis()->FindBin(sb4_xmax);
        int yi_sb4_l = h->GetYaxis()->FindBin(sb4_ymin);
        int yi_sb4_u = h->GetYaxis()->FindBin(sb4_ymax);


        double sclr_c2 = func(x_c_sb2,y_c_sb2,a,b,c);
        h_chunk_mod(h,hcopy2,sclr_c2,xi_sb2_l,xi_sb2_u,yi_sb2_l,yi_sb2_u);

        double sclr_c3 = func(x_c_sb3,y_c_sb3,a,b,c);
        h_chunk_mod(h,hcopy2,sclr_c3,xi_sb3_l,xi_sb3_u,yi_sb3_l,yi_sb3_u);

        double sclr_c4 = func(x_c_sb4,y_c_sb4,a,b,c);
        h_chunk_mod(h,hcopy2,sclr_c4,xi_sb4_l,xi_sb4_u,yi_sb4_l,yi_sb4_u);

        int chunk_size_x = xi_sb2_u-xi_sb2_l;
        int chunk_size_y = yi_sb2_u-yi_sb2_l;
        int n_chunks_x = nbinsx/chunk_size_x;
        int n_chunks_y = nbinsy/chunk_size_y;

        int xi_c_start = (chunk_size_x+1)/2;
        int yi_c_start = (chunk_size_y+1)/2;
        

        TH2D *hcopy3 = (TH2D*)hcopy2->Clone();
        hcopy3->SetName(Form("hHitPatternLayer%d_modified_coarse",i));
        hcopy3->Reset();
        for (int xni=0;xni<n_chunks_x;xni++){
            int xi = xni*chunk_size_x+xi_c_start;
            double x_center = h->GetXaxis()->GetBinCenter(xi);
            int x_l = xi-chunk_size_x/2;
            if (x_center>0){
                x_l+=1;
            }
            int x_u = x_l+chunk_size_x;
            if (xni==10){
                x_u+=1;    
            }
            for (int yni=0;yni<n_chunks_y;yni++){
                int yi = yni*chunk_size_y+yi_c_start;
                double y_center = h->GetYaxis()->GetBinCenter(yi);
                double sclr = func(x_center,y_center,a,b,c);
                int y_l = yi-chunk_size_y/2;
                if (y_center>0){
                    y_l+=1;
                }
                int y_u = y_l+chunk_size_y;
                if (yni==10){
                    y_u+=1;
                }
                h_chunk_mod(h,hcopy3,sclr,x_l,x_u,y_l,y_u);                
            }
        }



        double h_int_w = h->Integral("width");
        double hcopy_int_w = hcopy->Integral("width");
        double hcopy2_int_w = hcopy2->Integral("width");
        double hcopy3_int_w = hcopy3->Integral("width");

        std::cout<<h_int_w<<" "<<hcopy_int_w<<" "<<hcopy3_int_w<<"\n";
        std::cout<<h_int_w*tot_ent<<" "<<hcopy_int_w*tot_ent<<" "<<hcopy3_int_w*tot_ent<<"\n";

        otxt2<<Form("Layer %d",i)<<"\n";
        otxt2<<std::setw(12)<<"Est avg mult: "<<std::setw(12)<<h_int_w<<std::setw(12)<<hcopy_int_w<<std::setw(12)<<hcopy3_int_w<<"\n";
        otxt2<<std::setw(12)<<"Est tot mult: "<<std::setw(12)<<h_int_w*tot_ent<<std::setw(12)<<hcopy_int_w*tot_ent<<std::setw(12)<<hcopy3_int_w*tot_ent<<"\n";

        h->SetEntries(h_int_w*tot_ent);
        hcopy->SetEntries(hcopy_int_w*tot_ent);
        hcopy2->SetEntries(hcopy2_int_w*tot_ent);
        hcopy3->SetEntries(hcopy3_int_w*tot_ent);
        
        h_add->Add(h);

        h_add_fine->Add(hcopy);
        h_add_sb->Add(hcopy2);
        h_add_coarse->Add(hcopy3);

        fo1->cd();
        h->Write();
        hcopy->Write();
        hcopy2->Write();
        hcopy3->Write();

    }

    double h_add_int_w = h_add->Integral("width");
    double h_add_fine_int_w = h_add_fine->Integral("width");
    double h_add_coarse_int_w = h_add_coarse->Integral("width");

    std::cout<<"add"<<"\n";
    std::cout<<h_add_int_w<<" "<<h_add_fine_int_w<<" "<<h_add_coarse_int_w<<"\n";
    std::cout<<h_add_int_w*tot_ent<<" "<<h_add_fine_int_w*tot_ent<<" "<<h_add_coarse_int_w*tot_ent<<"\n";

    otxt2<<"add"<<"\n";
    otxt2<<std::setw(12)<<"Est avg mult: "<<std::setw(12)<<h_add_int_w<<std::setw(12)<<h_add_fine_int_w<<std::setw(12)<<h_add_coarse_int_w<<"\n";
    otxt2<<std::setw(12)<<"Est tot mult: "<<std::setw(12)<<h_add_int_w*tot_ent<<std::setw(12)<<h_add_fine_int_w*tot_ent<<std::setw(12)<<h_add_coarse_int_w*tot_ent<<"\n";

    double h_all_int_w = h_all->Integral("width");

    std::cout<<"all"<<"\n";
    std::cout<<h_all_int_w<<"\n";
    std::cout<<h_all_int_w*tot_ent<<"\n";

    otxt2<<"all"<<"\n";
    otxt2<<std::setw(12)<<"Est avg mult: "<<std::setw(12)<<h_all_int_w<<"\n";
    otxt2<<std::setw(12)<<"Est tot mult: "<<std::setw(12)<<h_all_int_w*tot_ent<<"\n";

    fo1->cd();
    h_add->Write();
    h_add_fine->Write();
    h_add_sb->Write();
    h_add_coarse->Write();


    otxt.close();
    otxt2.close();
    f->Close();
    fo1->Close();
    return;
}