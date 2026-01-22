#include "../myConst.h"
#include <TCutG.h>
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>

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




void recut(
    const std::string &cut_pat = "pattern",
    const std::string &pfile = "parA"

){
    for (int l=0;l<3;l++){
        TFile *f = new TFile((cut_pat+std::to_string(l)+".root").c_str(),"READ");
        TFile *fo1 = new TFile((cut_pat+std::to_string(l)+"_new.root").c_str(),"RECREATE");
        std::vector<double> aparam = parA_load((pfile+".txt").c_str())[0][l];
        std::cout<<aparam[0]<<" "<<aparam[1]<<" "<<aparam[2]<<"\n";

        int num_sec=0;

        TIter next(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next())) {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (cl->InheritsFrom("TCutG")) {
                num_sec++;
            }
        }
        for (int i=0;i<num_sec;i++){
            TCutG *cut = (TCutG*)f->Get(Form("cut_s%d",i));
            std::cout<<"original cut"<<"\n";
            cut->Print();
            int n = cut->GetN();
            TCutG *cut_new = new TCutG();
            cut_new->SetName(Form("cut_s%d",i));
            for (int j=0;j<n;j++){
                double x = cut->GetPointX(j);
                double y = cut->GetPointY(j);
                double y_new = aparam[0]*exp(aparam[1]*y)+aparam[2];
                cut_new->SetPoint(j,x,y_new);
            }
            std::cout<<"remade cut"<<"\n";
            cut_new->Print();
            fo1->cd();
            cut_new->Write();
        }
        fo1->Close();
        f->Close();
    }

    return;
}