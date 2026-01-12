#include <TFile.h>
#include <TCutG.h>
#include <string>
#include <fstream>
#include <iostream>

//for loading the graphical cuts made
void cut_load(
    const std::vector<std::string> &cuts = {"cut1","cut2"},
    const std::string &ofile = "output.root"
){
    TFile *f_out = new TFile(ofile.c_str(),"RECREATE");
    
    for (size_t i=0;i<cuts.size();i++){
        std:: string cut_exc = ".x ";
        gROOT->ProcessLine((cut_exc+cuts[i]).c_str());
        std::string cut_name = "cut_s" + to_string(i+1);
        TCutG *cut_i = (TCutG*)gROOT->FindObject("CUTG");
        cut_i->SetName(cut_name.c_str());
        f_out->cd();
        cut_i->Write();
        cut_i->Delete();
    }
    
    
    f_out->cd();
    f_out->Close();

    return ;
}
