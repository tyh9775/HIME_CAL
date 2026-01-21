#include "../myConst.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TStyle.h>
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <TMath.h>
#include <TCutG.h>


std::vector<std::vector<std::vector<double>>> parA_load(const std::string &s_est){
    std::vector<std::vector<double>> parA {};
    std::vector<std::vector<double>> par_err {};

    std::ifstream f_in(s_est);
    std::string line;
    while (std::getline(f_in, line)) {
        std::vector<double> row;
        std::vector<double> row_err;
        std::istringstream iss(line);
        std::string token;
        
        while (iss >> token) {
            // Skip "+/-" tokens
            if (token == "+/-") {
                iss >> token;
                double value = std::stod(token);
                row_err.push_back(value);
            }
            
            else {
                double value = std::stod(token);
                row.push_back(value);
            } 
        }
        
        if (!row.empty()) {
            parA.push_back(row);
            par_err.push_back(row_err);
        }
    }
    
    f_in.close();
    return {parA, par_err};
}



double f_cal(const double *param, double ToT, const std::vector<double> &parA){
    double A_adj {param[0]};
    double k_adj {param[1]};
    double B_adj {param[2]};

    double A=parA[0];
    double k=parA[1];
    double B=parA[2];


    double En {A_adj*exp(k_adj*ToT)+B_adj};
    double ToT_cal {log((En-B)/A)/k};

    return ToT_cal;
}


TGraph *f_change(const double *param,const TGraphErrors *g,const std::vector<double> &parA){
    std::vector<double> ToT_cal {};
    std::vector<double> ToF {};
    int n = g->GetN();
    for (int i=0;i<n;i++){
        double tof = g->GetPointX(i);
        double tot = g->GetPointY(i);
        ToF.push_back(tof);
        ToT_cal.push_back(f_cal(param,tot,parA));
    }
    TGraph *gnew = new TGraph(n,ToF.data(),ToT_cal.data()); 

    return gnew;
}



struct DataPoint {
    double x;
    double y;
    double error;
};

struct MatchedPoint {
    double x;
    double y1;
    double err1;
    double y2;
    double err2;
};


std::vector<std::vector<MatchedPoint>> data_sort(const std::vector<std::vector<DataPoint>> &data_set_1, const std::vector<std::vector<DataPoint>> &data_set_2){
    std::vector<std::vector<MatchedPoint>> matched_points_set {};
    for (size_t i=0;i<data_set_1.size();i++){
        std::vector<DataPoint> data_1 {data_set_1[i]};
        std::vector<DataPoint> data_2 {data_set_2[i]};

        std::unordered_map<double, std::pair<double, double>> map2;
        for (const auto &point : data_2){
            map2[point.x] = std::make_pair(point.y, point.error);
        }

        std::vector<MatchedPoint> matched_points {};
        for (const auto &point1 :data_1){
            auto it = map2.find(point1.x);
            if (it != map2.end()) {
                matched_points.push_back({point1.x, point1.y, point1.error, it->second.first, it->second.second});
            }
        } 
        matched_points_set.push_back(matched_points);
    
    }
    
    return matched_points_set;
}




double f_opt(const double *param, const std::vector<double>& ToT1, const std::vector<double>& ToT2,const std::vector<double>& err1, const std::vector<double>& err2, const std::vector<double> &parA){
    double A_adj {param[0]};
    double k_adj {param[1]};
    double B_adj {param[2]};    
    double sum {0};

    double A=parA[0];
    double k=parA[1];
    double B=parA[2];

    double cd = (A+B)-(A_adj+B_adj);
    

    for (size_t i=0;i<ToT1.size();i++){
        sum+=cd*cd;
        double En1 {A*exp(k*ToT1[i])+B};

        double En2 {A_adj*exp(k_adj*ToT2[i])+B_adj};
        double diff {En2-En1};

        double dA1_dT1 {A*k*exp(k*ToT1[i])};

        double dA2_dT2 {A_adj*k_adj*exp(k_adj*ToT2[i])};

        double sig_A1 {};
        double sig_A2 {};

        if (err1.size()==0){
            sig_A1=1.0;
        } else {
            sig_A1=dA1_dT1*err1[i];
        }
        if (err2.size()==0){
            sig_A2=1.0;
        } else {
            sig_A2=dA2_dT2*err2[i];
        }

        double sig2 {sig_A1*sig_A1+sig_A2*sig_A2};
        
        if (sig2<1e-8){
            sig2=1e-8;
        }

        sum+=diff*diff/sig2; 
        
    }

    return pow(sum,0.5);
}

std::unique_ptr<ROOT::Math::Minimizer> min_func(const std::vector<double>& ToT1, const std::vector<double>& ToT2,const std::vector<double>& err1, const std::vector<double>& err2,const std::vector<double> &parA){
    std::unique_ptr<ROOT::Math::Minimizer> min_1(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    );
    min_1->SetMaxFunctionCalls(1000000);
    min_1->SetTolerance(1e-6);

    auto f_opt_1 = [&](const double *param){
        return f_opt(param,ToT1,ToT2,err1,err2,parA);
    };
    
    ROOT::Math::Functor f1(f_opt_1,3);

    double A=parA[0];
    double k=parA[1];
    double B=parA[2];


    min_1->SetFunction(f1);
    min_1->SetVariable(0,"A_adj",A,0.1); 
    min_1->SetVariable(1,"k_adj",k,0.1);
    min_1->SetVariable(2,"B_adj",B,0.01);


    min_1->SetVariableLimits(0,0.0,1.0);
    min_1->SetVariableLimits(1,0.0,1.0);
    min_1->SetVariableLimits(2,0.0,1.0);

    min_1->Minimize();

    const double *p1 = min_1->X();
    const double *p1_err = min_1->Errors();
    if (min_1->Status()==0){
        std::cout<<"successfully converged"<<"\n";
    } else if (min_1->Status()==1){
        std::cout<<"cov mat forced pos-def"<<'\n';
    } else if (min_1->Status()==2){
        std::cout<<"Hesse mat (2nd deriv) invalid"<<'\n';
    } else if (min_1->Status()==3){
        std::cout<<"Est Dist to Min too large"<<'\n';
    } else if (min_1->Status()==4){
        std::cout<<"max calls exceeded"<<'\n';
    } else if (min_1->Status()==5){
        std::cout<<"parameters at limit"<<'\n';
    } else {
        std::cout<<"unknown error"<<'\n';
    }
    
    std::cout<<"A_adj="<<p1[0]<<" +/- " <<p1_err[0]<<"\n";
    std::cout<<"k_adj="<<p1[1]<<" +/- " <<p1_err[1]<<"\n";
    std::cout<<"B_adj="<<p1[2]<<" +/- " <<p1_err[2]<<"\n"; 
    std::cout<<"\n";
    
    return min_1;
}


std::vector<std::vector<double>> data_load(const TGraphErrors *g){
    std::vector<double> tof {};
    std::vector<double> tot {};
    std::vector<double> err {};

    int n = g->GetN();
    for (int j=0;j<n;j++){
        double tof_j = g->GetPointX(j);
        double tot_j = g->GetPointY(j);
        double err_j = g->GetErrorY(j);
        tof.push_back(tof_j);
        tot.push_back(tot_j);
        err.push_back(err_j);
    }
    return {tof,tot,err};
}




void ToT_cal(
    const std::string &file = "input.root",
    const std::string &ofile = "output",
    const std::string &s_est = "single_est"
){
    //open loaded files
    auto *f = new TFile(file.c_str(),"READ");
    //create output file
    auto *fo1 = new TFile((ofile+".root").c_str(),"RECREATE");
    
    std::vector<int> hstarts {hstart1,hstart2,hstart3};
    std::vector<int> hstops {hstop1,hstop2,hstop3};
    std::vector<int> hbases {hbase1,hbase2,hbase3};

    //load estimated parameters
    std::vector<std::vector<double>> parA_mat = parA_load((s_est+".txt"))[0];
    std::vector<std::vector<double>> parA_err_mat = parA_load((s_est+".txt"))[1];

    std::cout<<"Loaded calibration parameters from "<<s_est+".txt"<<"\n";
    //for saving the parameters
    std::ofstream param_output;
    param_output.open((ofile+".txt").c_str(),std::ios_base::trunc);

    param_output << std::left << std::setw(8) << "Det" 
                 << std::setw(28) << "A" 
                 << std::setw(28) << "k" 
                 << std::setw(28) << "B" << "\n";

    for (int L=0;L<3;L++){
        std::vector<double> parA=parA_mat[L];
        std::vector<double> parA_err=parA_err_mat[L];
        TCanvas *c0 = (TCanvas*)f->Get(Form("hcopy%d",hbases[L]));
        int n0 = c0->GetListOfPrimitives()->GetSize();

        TH2D *h0 = nullptr;
        for (int i=0; i<n0; i++){
            TObject* obj = c0->GetListOfPrimitives()->At(i);
            if (obj->InheritsFrom("TH2")) {
                h0 = (TH2D*)obj;
                break;
            }
        }
        

        std::vector<TGraphErrors*> gbase {};

        std::vector<std::vector<DataPoint>> data_set_0 {};

        std::vector<double> starts {};
        std::vector<double> stops {};

        for (int i=0;i<n0;i++){
            TGraphErrors *g0 = nullptr;
            TObject* obj = c0->GetListOfPrimitives()->At(i);
            if (obj->InheritsFrom("TGraph")) {
                g0 = (TGraphErrors*)obj;
            } else {
                continue;
            }
            TGraphErrors *gc0 = (TGraphErrors*)g0->Clone();
            gbase.push_back(gc0);

            std::vector<std::vector<double>> data0 = data_load(gc0);
            std::vector<DataPoint> data_points0 {};
            for (size_t j=0; j<data0[0].size(); j++){
                data_points0.push_back({data0[0][j], data0[1][j], data0[2][j]});
            }
            data_set_0.push_back(data_points0);
            starts.push_back(data0[0][0]);
            stops.push_back(data0[0][data0[0].size()-1]);
        }
        

        for (int i=hstarts[L];i<=hstops[L];i++){
            if(i==hbases[L]){
                param_output << std::setw(8) << i
                << std::setw(11) << parA[0] << "+/- " << std::setw(12) << parA_err[0] << " "
                << std::setw(11) << parA[1] << "+/- " << std::setw(12) << parA_err[1] << " "
                << std::setw(11) << parA[2] << "+/- " << std::setw(12) << parA_err[2] << "\n";
                continue;
            }

            TCanvas *ci = (TCanvas*)f->Get(Form("hcopy%d",i));
            TH2D *hi = nullptr;
            int ni = ci->GetListOfPrimitives()->GetSize();
            for (int j=0; j<ni; j++) {
                TObject* obj = ci->GetListOfPrimitives()->At(j);
                if (obj->InheritsFrom("TH2")) {
                    hi = (TH2D*)obj;
                    break;
                }
            }
            
            std::vector<TGraphErrors*> g_i {};

            std::vector<std::vector<DataPoint>> data_set_i {};


            for (int j=0;j<ni;j++){
                TGraphErrors *g = nullptr;
                TObject* obj = ci->GetListOfPrimitives()->At(j);
                if (obj->InheritsFrom("TGraph")) {
                    g = (TGraphErrors*)obj;
                } else {
                    continue;
                }
                TGraphErrors *gc = (TGraphErrors*)g->Clone();
                g_i.push_back(gc);
                std::vector<std::vector<double>> data = data_load(gc);

                std::vector<DataPoint> data_points {};
                for (size_t j=0; j<data[0].size(); j++){
                    data_points.push_back({data[0][j], data[1][j], data[2][j]});
                }
                data_set_i.push_back(data_points);
            }


            std::vector<std::vector<MatchedPoint>> matched_points_i = data_sort(data_set_0, data_set_i);


            std::vector<double> tot0 {};
            std::vector<double> tot_i {};
            std::vector<double> err0 {};
            std::vector<double> err_i {};

            for (size_t j=0;j<matched_points_i.size();j++){
                const auto &matched_points = matched_points_i[j];
                for (int k=0;k<matched_points.size();k++){
                    const auto &point = matched_points[k];
                    tot0.push_back(point.y1);
                    tot_i.push_back(point.y2);
                    err0.push_back(point.err1);
                    err_i.push_back(point.err2);
                }
            }


      

            std::cout<<"*******************************"<<'\n';
            std::cout<<"Param for det"<<hbases[L]<<" and det"<<i<<"\n";
            auto min1 = min_func(tot0,tot_i,err0,err_i,parA);
            std::cout<<"*******************************"<<"\n";
            const double *p_i = min1->X();
            const double *p_i_err = min1->Errors();
            param_output << std::setw(8) << i
             << std::setw(11) << p_i[0] << "+/- " << std::setw(12) << p_i_err[0] << " "
             << std::setw(11) << p_i[1] << "+/- " << std::setw(12) << p_i_err[1] << " "
             << std::setw(11) << p_i[2] << "+/- " << std::setw(12) << p_i_err[2]
             << "\n";

            
            TCanvas *c_comp = new TCanvas(Form("c_comp%d",i),Form("comparison det %d",i));
            c_comp->Divide(2,2);
            c_comp->cd(1);
            gPad->SetLogz();
            h0->Draw();
            for (size_t j=0;j<gbase.size();j++){
                gbase[j]->Draw("SAME");
                gPad->Update();
            }
            
            c_comp->cd(2);
            gPad->SetLogz();
            hi->Draw();
            for (size_t j=0;j<g_i.size();j++){
                g_i[j]->Draw("SAME");
                gPad->Update();
            }
            
            c_comp->cd(3);
            gPad->SetLogz();
            h0->Draw();
            std::vector<TGraph*> gnew_list {};
            for (size_t j=0;j<g_i.size();j++){
                TGraphErrors *gj = (TGraphErrors*)g_i[j]->Clone();
                TGraph *gnew = f_change(p_i,gj,parA);
                TF1 *fbase = (TF1*)gbase[j]->GetListOfFunctions()->At(0);
                int npar = fbase->GetNpar();
                TF1 *ffit = new TF1(Form("fit_deg_%d",npar-1),Form("pol%d",npar-1),fbase->GetXmin(),fbase->GetXmax());
                ffit->SetLineColor(kBlue);
                gnew->Fit(ffit,"QRB");
                gnew->Draw("SAME");
                fbase->Draw("SAME");
                gnew_list.push_back(gnew);
                gPad->Update();
                if (j==0){
                    TLegend *legend1 = new TLegend(0.78, 0.6, 0.98, 0.7);
                    legend1->AddEntry(fbase,"base","l");
                    legend1->AddEntry(ffit,"calibrated","l");
                    legend1->Draw();
                    
                }
            }

            c_comp->cd(4);
            gPad->SetLogz();
            hi->Draw();
            for (size_t j=0;j<gnew_list.size();j++){
                gnew_list[j]->Draw("SAME");
                TF1 *ffit = (TF1*)gnew_list[j]->GetListOfFunctions()->At(0);
                TF1 *fdet=(TF1*)g_i[j]->GetListOfFunctions()->At(0);
                fdet->Draw("SAME");
                gPad->Update();
                if (j==0){
                    TLegend *legend1 = new TLegend(0.78, 0.6, 0.98, 0.7);
                    legend1->AddEntry(fdet,"detector","l");
                    legend1->AddEntry(ffit,"calibrated","l");
                    legend1->Draw();
                }
            }

            fo1->cd();
            c_comp->Write();

        }

        
    }
    
    param_output.close();    
    
    fo1->Close();
    f->Close();
    return;
}
        
