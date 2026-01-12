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


std::vector<TH1D*> h_rb(TFile *f, std::string hname){
    std::vector<TH1D*> h_list {};

    TH1D *hBarHit = (TH1D*)f->Get(hname.c_str())->Clone();
    h_list.push_back(hBarHit);
    TH1D *hBarHit_rb = (TH1D*)f->Get((hname+"_rb").c_str())->Clone();
    h_list.push_back(hBarHit_rb);
    TH1D *hBarHit_rb2 = (TH1D*)f->Get((hname+"_rb2").c_str())->Clone();
    h_list.push_back(hBarHit_rb2);
    TH1D *hBarHit_rb3 = (TH1D*)f->Get((hname+"_rb3").c_str())->Clone();
    h_list.push_back(hBarHit_rb3);
    
    return h_list;
}


std::vector<std::vector<std::vector<TH1D*>>> h_load(std::vector<int> sbl,TFile *fsb, TFile *fnsb,std::string hname){
    std::vector<std::vector<TH1D*>> h_set {}; //detectors inc sb
    std::vector<std::vector<TH1D*>> h2_set {}; //detectors without sb

    std::vector<std::vector<TH1D*>> ha_set {}; //inc sb
    std::vector<std::vector<TH1D*>> ha2_set {}; //without sb

    std::vector<std::vector<TH1D*>> h_set_rnd {}; //inc sb
    std::vector<std::vector<TH1D*>> h2_set_rnd {}; //without sb

    for (size_t i=0;i<sbl.size();i++){
        int det_i=sbl[i];
        std::string hname_det = hname + to_string(det_i);
        h_set.push_back(h_rb(fsb,hname_det));
        h2_set.push_back(h_rb(fnsb,hname_det));
        ha_set.push_back(h_rb(fsb,hname_det+"_area"));
        ha2_set.push_back(h_rb(fnsb,hname_det+"_area"));
        h_set_rnd.push_back(h_rb(fsb,hname_det+"_rnd"));
        h2_set_rnd.push_back(h_rb(fnsb,hname_det+"_rnd"));
    }

    return {h_set,h2_set,ha_set,ha2_set,h_set_rnd,h2_set_rnd};
}

double err_calc(double no_sb_val,double inc_sb_val,double no_sb_err,double inc_sb_err){
    return inc_sb_err/no_sb_val + inc_sb_val*no_sb_err/pow(no_sb_val,2);
}




void shadow_comp_v4(
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
    std::ofstream sum_output_rb2;
    sum_output_rb2.open((dir+"rb2.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb3;
    sum_output_rb3.open((dir+"rb3.txt").c_str(),std::ios_base::trunc);

    std::ofstream sum_output_og_a;
    sum_output_og_a.open((dir+"ogb_a.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb_a;
    sum_output_rb_a.open((dir+"rb_a.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb2_a;
    sum_output_rb2_a.open((dir+"rb2_a.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb3_a;
    sum_output_rb3_a.open((dir+"rb3_a.txt").c_str(),std::ios_base::trunc);


    std::ofstream sum_output_og_rnd;
    sum_output_og_rnd.open((dir+"ogb_rnd.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb_rnd;
    sum_output_rb_rnd.open((dir+"rb_rnd.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb2_rnd;
    sum_output_rb2_rnd.open((dir+"rb2_rnd.txt").c_str(),std::ios_base::trunc);
    std::ofstream sum_output_rb3_rnd;
    sum_output_rb3_rnd.open((dir+"rb3_rnd.txt").c_str(),std::ios_base::trunc);

    std::vector<std::reference_wrapper<std::ofstream>> sum_outputs_og {sum_output_og,sum_output_rb,sum_output_rb2,sum_output_rb3};
    std::vector<std::reference_wrapper<std::ofstream>> sum_outputs_a {sum_output_og_a,sum_output_rb_a,sum_output_rb2_a,sum_output_rb3_a};
    std::vector<std::reference_wrapper<std::ofstream>> sum_outputs_rnd {sum_output_og_rnd,sum_output_rb_rnd,sum_output_rb2_rnd,sum_output_rb3_rnd};

    std::vector<std::vector<int>> sb_mat {sb2,sb3,sb4};
    std::vector<std::vector<double>> sb_posx_mat {{sb2_xmin,sb2_xmax},{sb3_xmin,sb3_xmax},{sb4_xmin,sb4_xmax}};
    std::vector<std::vector<double>> sb_posy_mat {{sb2_ymin,sb2_ymax},{sb3_ymin,sb3_ymax},{sb4_ymin,sb4_ymax}};

    for (size_t i=0;i<files.size();i++){ //i should be shadow bar number minus 2

        std::vector<std::vector<std::vector<TH1D*>>> h_sb {h_load(sb_mat[i],f_sb[i],f_no_sb[i], "hBarHit")};

        for (size_t hi=0;hi<h_sb.size();hi++){
            if(hi%2!=0){
                continue;
            }
            std::vector<std::vector<TH1D*>> h_sb_set = h_sb[hi];
            std::vector<std::vector<TH1D*>> h_no_sb_set = h_sb[hi+1];

            for (size_t hj=0;hj<h_sb_set.size();hj++){ //hj is detector iterator
                std::vector<TH1D*> h_list = h_sb_set[hj];
                std::vector<TH1D*> h2_list = h_no_sb_set[hj];

                for (size_t hk=0;hk<h_list.size();hk++){ //hk is rb iterator
                    int bin_sb_min = 0;
                    int bin_sb_max = 0;

                    if (hj==0 || hj==2){
                        bin_sb_min = h_list[hk]->FindBin(sb_posx_mat[i][0]);
                        bin_sb_max = h_list[hk]->FindBin(sb_posx_mat[i][1]);
                    } else if (hj==1){
                        bin_sb_min = h_list[hk]->FindBin(sb_posy_mat[i][0]);
                        bin_sb_max = h_list[hk]->FindBin(sb_posy_mat[i][1]);
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

                    double err {0.0};
                    double err2 {0.0};

                    double h_integral = h_list[hk]->IntegralAndError(bin_sb_min,bin_sb_max,err,"width");
                    double h2_integral = h2_list[hk]->IntegralAndError(bin_sb_min,bin_sb_max,err2,"width");

                    sum_output<<hname<<"\n";
                    sum_output<<"sb integral ("<<h_list[hk]->GetXaxis()->GetBinLowEdge(bin_sb_min)<<", "<<h_list[hk]->GetXaxis()->GetBinUpEdge(bin_sb_max)<<"): "<<h_integral<<" +/- "<<err<<"\n";
                    sum_output<<"no sb integral ("<<h2_list[hk]->GetXaxis()->GetBinLowEdge(bin_sb_min)<<", "<<h2_list[hk]->GetXaxis()->GetBinUpEdge(bin_sb_max)<<"): "<<h2_integral<<" +/- "<<err2<<"\n";
                    sum_output<<"yield percentage: "<<(1.0 - h_integral/h2_integral)*100.0<<"% +/- "<<err_calc(h2_integral,h_integral,err2,err)*100.0<<"% \n";
                    sum_output<<"\n";

                    if (hi==0){
                        sum_outputs_og[hk].get()<<1.0 - h_integral/h2_integral<<"+/-"<<err_calc(h2_integral,h_integral,err2,err)<<"\n";
                    } else if (hi==2){
                        sum_outputs_a[hk].get()<<1.0 - h_integral/h2_integral<<"+/-"<<err_calc(h2_integral,h_integral,err2,err)<<"\n";
                    } else if (hi==4){
                        sum_outputs_rnd[hk].get()<<1.0 - h_integral/h2_integral<<"+/-"<<err_calc(h2_integral,h_integral,err2,err)<<"\n";
                    }
                }

            }

        }

    }


    sum_output.close();
    
    sum_output_og.close();
    sum_output_rb.close();
    sum_output_rb2.close();
    sum_output_rb3.close();

    sum_output_og_a.close();
    sum_output_rb_a.close();
    sum_output_rb2_a.close();
    sum_output_rb3_a.close();

    sum_output_og_rnd.close();
    sum_output_rb_rnd.close();
    sum_output_rb2_rnd.close();
    sum_output_rb3_rnd.close();

    fo1->Close();

    return;
}
