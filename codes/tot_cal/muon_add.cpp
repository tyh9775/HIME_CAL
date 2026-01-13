#include "../myConst.h"
#include <TFile.h>



void muon_add(
    const std::string &input = "input",
    const std::string &ofile = "output.root"
){
    std::vector<std::string> files {};
    for (int i=0;i<15;i++){
        std::string f_i = input;
        if (i<10){
            f_i += "0"+to_string(i)+".root";
        } else {
            f_i += to_string(i)+".root";
        }
        files.push_back(f_i);
    }
    TFile *fo1 = new TFile(ofile.c_str(),"RECREATE");

    TH2F *hTotVsTDiff = new TH2F("hTotVsTDiff","",320,-40,40,100,0,50);
    hTotVsTDiff->Sumw2();

    std::vector<TH2F*> hTotVsTDiffMod {};
    for (int i = 0; i < 72; i++) {
		TH2F *hTotVsTDiff_i = new TH2F(Form("hTotVsTDiff_mod_%d", i), "", 320,-40,40,100,0,50);
		hTotVsTDiff_i->Sumw2();
		hTotVsTDiffMod.push_back(hTotVsTDiff_i);

	}

    for (size_t i=0;i<files.size();i++){
        TFile *f = new TFile(files[i].c_str(),"READ");
        for (int j=0;j<72;j++){
            TH2F *h = (TH2F*)f->Get(Form("hTotVsTDiff_module_0%02d",j))->Clone();
            hTotVsTDiffMod[j]->Add(h);
            hTotVsTDiff->Add(h);
        }
        
    }

    fo1->cd();
    hTotVsTDiff->Write();

    for (int i=0;i<72;i++){
        fo1->cd();
        hTotVsTDiffMod[i]->Write();
    }
    fo1->Close();
    return;
}


