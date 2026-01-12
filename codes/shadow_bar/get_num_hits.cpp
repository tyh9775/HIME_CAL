#include "spirit.h"


void get_num_hits(
    const std::vector<std::string> &files = {"input1.root","input2.root"}, 
    const std::string &ofile = "num_entries.txt"
){
    std::ofstream output;
    output.open(ofile.c_str(),std::ios_base::trunc);  
    for (size_t i=0;i<files.size();i++){
        TFile *f = TFile::Open(files[i].c_str(),"READ");
        TH2D *hHitPattern = (TH2D*)f->Get("hHitPattern");
        double nhits = hHitPattern->GetEntries();
        std::cout<<"Number of Hits:"<<nhits<<"\n";
        output<<nhits<<"\n";
    }
    output.close();
}