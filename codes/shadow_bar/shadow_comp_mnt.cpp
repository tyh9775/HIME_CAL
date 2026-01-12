#include "myConst3.h"
#include "spirit.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <filesystem>
#include <functional>


std::vector<TH2D*> h_rb(TFile *f, std::string hname){
    std::vector<TH2D*> h_list {};

    TH2D *hBarHit = (TH2D*)f->Get(hname.c_str())->Clone();
    h_list.push_back(hBarHit);
    TH2D *hBarHit_rb = (TH2D*)f->Get((hname+"_rb").c_str())->Clone();
    h_list.push_back(hBarHit_rb);

    return h_list;
}


std::vector<std::vector<std::vector<TH2D*>>> h_load(std::vector<int> sbl,TFile *fsb, TFile *fnsb,std::string hname){
    std::vector<std::vector<TH2D*>> h_set {}; //detectors inc sb
    std::vector<std::vector<TH2D*>> h2_set {}; //detectors without sb


    for (size_t i=0;i<sbl.size();i++){
        int det_i=sbl[i];
        std::string hname_det = hname + to_string(det_i);
        h_set.push_back(h_rb(fsb,hname_det));
        h2_set.push_back(h_rb(fnsb,hname_det));
    }

    return {h_set,h2_set};
}

double int_calc(TH2D *h,int p_min,int p_max, int bin_sb_min, int bin_sb_max, double &err){
    TH1D *hp = h->ProjectionX("hp",p_min,p_max);

    double integral = hp->IntegralAndError(bin_sb_min,bin_sb_max,err,"width");

    return integral;
}


double err_calc(double no_sb_val,double inc_sb_val,double no_sb_err,double inc_sb_err){
    return inc_sb_err/no_sb_val + inc_sb_val*no_sb_err/pow(no_sb_val,2);
}


void shadow_comp_v5(
    //const std::string &file, //file with all shadow bars
    const std::vector<std::string> &files_inc = {"inc1.root","inc2.root"}, //files including certain shadow bars
    const std::vector<std::string> &files = {"input1.root","input2.root"}, //files without certain shadow bars
    const std::string &ofile = "output"
){

    //TFile *f_all = TFile::Open(file.c_str(),"READ");
    std::vector<TFile*> f_sb;
    std::vector<TFile*> f_no_sb;
    for (size_t i=0;i<files_inc.size();i++){
        TFile *f = TFile::Open(files_inc[i].c_str(),"READ");
        f_sb.push_back(f);
    }

    for (size_t i=0;i<files.size();i++){
        TFile *f = TFile::Open(files[i].c_str(),"READ");
        f_no_sb.push_back(f);
    }
    
    auto *fo1 = new TFile((ofile+".root").c_str(),"RECREATE");

    std::string dir = ofile + "_res/";
    std::filesystem::create_directory(dir.c_str());
    std::ofstream sum_output;
    sum_output.open((dir+"all.txt").c_str(),std::ios_base::trunc);
    
    std::ofstream sum_output_og;
    sum_output_og.open((dir+"ogb.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb;
    sum_output_rb.open((dir+"rb.txt").c_str(),std::ios_base::trunc);


    std::vector<std::reference_wrapper<std::ofstream>> sum_outputs_og {sum_output_og,sum_output_rb};


    std::vector<std::vector<int>> sb_mat {sb2,sb3,sb4};
    std::vector<std::vector<double>> sb_posx_mat {{sb2_xmin,sb2_xmax},{sb3_xmin,sb3_xmax},{sb4_xmin,sb4_xmax}};
    std::vector<std::vector<double>> sb_posy_mat {{sb2_ymin,sb2_ymax},{sb3_ymin,sb3_ymax},{sb4_ymin,sb4_ymax}};

    for (size_t i=0;i<files.size();i++){ //i should be shadow bar number minus 2

        std::vector<std::vector<std::vector<TH2D*>>> h_sb {h_load(sb_mat[i],f_sb[i],f_no_sb[i], "hBarHit_mnt")};

        for (size_t hi=0;hi<h_sb.size();hi++){
            if(hi%2!=0){
                continue;
            }
            std::vector<std::vector<TH2D*>> h_sb_set = h_sb[hi];
            std::vector<std::vector<TH2D*>> h_no_sb_set = h_sb[hi+1];

            for (size_t hj=0;hj<h_sb_set.size();hj++){ //hj is detector iterator
                std::vector<TH2D*> h_list = h_sb_set[hj];
                std::vector<TH2D*> h2_list = h_no_sb_set[hj];
                for (size_t hk=0;hk<h_list.size();hk++){ //hk is rb iterator
                    int bin_sb_min = 0;
                    int bin_sb_max = 0;

                    if (hj==0 || hj==2){
                        bin_sb_min = h_list[hk]->GetXaxis()->FindBin(sb_posx_mat[i][0]);
                        bin_sb_max = h_list[hk]->GetXaxis()->FindBin(sb_posx_mat[i][1]);
                    } else if (hj==1){
                        bin_sb_min = h_list[hk]->GetXaxis()->FindBin(sb_posy_mat[i][0]);
                        bin_sb_max = h_list[hk]->GetXaxis()->FindBin(sb_posy_mat[i][1]);
                    }

                    const char *hn = h_list[hk]->GetName();
                    std::string hname = std::string(hn);
                    std::string hname1 = hname +"_sb_"+to_string(i+2);
                    std::string hname2 = hname +"_no_sb_"+to_string(i+2);

                    h_list[hk]->SetName(hname1.c_str());
                    h2_list[hk]->SetName(hname2.c_str());

                    fo1->cd();
                    h_list[hk]->Write();
                    h2_list[hk]->Write();

                    int n_bins_p = h_list[hk]->GetYaxis()->GetNbins();
                    
                    TGraphErrors *g_yield_fine = new TGraphErrors();
                    g_yield_fine->SetName((hname+"_yield_sb_"+to_string(i+2)).c_str());
                    
                    TGraphErrors *g_yield_coarse = new TGraphErrors();
                    g_yield_coarse->SetName((hname+"_yield_coarse_sb_"+to_string(i+2)).c_str());
                    int point_i =0;
                    for (int p = 1; p <= n_bins_p; p+=10){
                        double err_sb {0.0};
                        double err_no_sb {0.0};

                        double integral_sb = int_calc(h_list[hk],p,p+10,bin_sb_min,bin_sb_max,err_sb);
                        double integral_no_sb = int_calc(h2_list[hk],p,p+10,bin_sb_min,bin_sb_max,err_no_sb);

                        if (integral_no_sb == 0){
                            point_i++;
                            continue;
                        }

                        double yield_frac = (1.0 - integral_sb/integral_no_sb);
                        double yield_err = err_calc(integral_no_sb,integral_sb,err_no_sb,err_sb);
                        if (yield_frac < 0.0){
                            yield_frac = 0.0;
                            yield_err = 0.0;
                        }


                        g_yield_fine->SetPoint(point_i, h_list[hk]->GetYaxis()->GetBinCenter(p), yield_frac);
                        g_yield_fine->SetPointError(point_i, 0.0, yield_err);
                        point_i++;
                    }


                    TH1D *h_sec[5];
                    TH1D *h2_sec[5];

                    for (int pi = 0; pi <5; pi++){
                        int p_start = 1 + pi*50;;
                        int p_end = p_start + 50;
                        double err_sb {0.0};
                        double err_no_sb {0.0};

                        double integral_sb = int_calc(h_list[hk],p_start,p_end,bin_sb_min,bin_sb_max,err_sb);
                        double integral_no_sb = int_calc(h2_list[hk],p_start,p_end,bin_sb_min,bin_sb_max,err_no_sb);
                        
                        h_sec[pi] = h_list[hk]->ProjectionX(Form("h_sec_%d",pi),p_start,p_end);
                        h2_sec[pi] = h2_list[hk]->ProjectionX(Form("h2_sec_%d",pi),p_start,p_end);

                        h_sec[pi]->SetName((hname1+"_sec_"+to_string(pi)).c_str());
                        h2_sec[pi]->SetName((hname2+"_sec_"+to_string(pi)).c_str());

                        if (integral_no_sb == 0){
                            continue;
                        }
                        double yield_frac = (1.0 - integral_sb/integral_no_sb);
                        double yield_err = err_calc(integral_no_sb,integral_sb,err_no_sb,err_sb);
                        if (yield_frac < 0.0){
                            yield_frac = 0.0;
                            yield_err = 0.0;
                        }
                        g_yield_coarse->SetPoint(pi, (h_list[hk]->GetYaxis()->GetBinCenter(p_start) + h_list[hk]->GetYaxis()->GetBinCenter(p_end))/2.0, yield_frac);
                        g_yield_coarse->SetPointError(pi, 0.0, yield_err); 

                        fo1->cd();
                        h_sec[pi]->Write();
                        h2_sec[pi]->Write();
                    }

                    g_yield_fine->SetTitle("Shadow Bar Yield Fraction;Momentum (MeV/c);Yield Fraction");
                    g_yield_fine->Draw("AP");
                    g_yield_fine->SetMarkerStyle(20);
                    g_yield_fine->SetMarkerSize(0.8);
                    

                    fo1->cd();
                    g_yield_fine->Write();
                    g_yield_coarse->Write();

                    double err {0.0};
                    double err2 {0.0};

                    double h_integral = int_calc(h_list[hk],1,n_bins_p,bin_sb_min,bin_sb_max,err);
                    double h2_integral = int_calc(h2_list[hk],1,n_bins_p,bin_sb_min,bin_sb_max,err2);

                    sum_output<<hname<<" over all momentum \n";
                    sum_output<<"sb integral ("<<h_list[hk]->GetXaxis()->GetBinLowEdge(bin_sb_min)<<", "<<h_list[hk]->GetXaxis()->GetBinUpEdge(bin_sb_max)<<"): "<<h_integral<<" +/- "<<err<<"\n";
                    sum_output<<"no sb integral ("<<h2_list[hk]->GetXaxis()->GetBinLowEdge(bin_sb_min)<<", "<<h2_list[hk]->GetXaxis()->GetBinUpEdge(bin_sb_max)<<"): "<<h2_integral<<" +/- "<<err2<<"\n";
                    sum_output<<"yield percentage: "<<(1.0 - h_integral/h2_integral)*100.0<<"% +/- "<<err_calc(h2_integral,h_integral,err2,err)*100.0<<"% \n";
                    sum_output<<"\n";

                    if (hi==0){
                        sum_outputs_og[hk].get()<<1.0 - h_integral/h2_integral<<"+/-"<<err_calc(h2_integral,h_integral,err2,err)<<"\n";
                    } 
                }

            }

        }

    }


    sum_output.close();
    
    sum_output_og.close();
    sum_output_rb.close();

    fo1->Close();

    return;
}
