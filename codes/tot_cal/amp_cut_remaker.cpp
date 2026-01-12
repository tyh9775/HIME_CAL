#include "myConst3.h"
#include <TCutG.h>
#include <TFile.h>
#include <TKey.h>
#include <TClass.h>

void amp_cut_remaker(
    const std::vector<std::string> &inputs = {"cuts1.root","cuts2.root","cuts3.root"},
    const std::vector<std::string> &outputs = {"cuts1_remade.root","cuts2_remade.root","cuts3_remade.root"},
    const std::vector<std::vector<double>> &parA_mat = {{a1,a2,ka},{a1,a2,ka},{a1,a2,ka}}
){
    for (int l=0;l<3;l++){
        TFile *f = new TFile(inputs[l].c_str(),"READ");
        TFile *fo1 = new TFile(outputs[l].c_str(),"RECREATE");
        std::vector<double> aparam = parA_mat[l];

        int num_sec=0;

        TIter next(f->GetListOfKeys());
        TKey *key;
        while ((key = (TKey*)next())) {
            TClass *cl = gROOT->GetClass(key->GetClassName());
            if (cl->InheritsFrom("TCutG")) {
                num_sec++;
            }
        }
        for (int i=1;i<=num_sec;i++){
            TCutG *cut = (TCutG*)f->Get(Form("cut_s%d",i));
            std::cout<<"original cut"<<"\n";
            cut->Print();
            int n = cut->GetN();
            TCutG *cut_new = new TCutG();
            cut_new->SetName(Form("cut_s%d",i));
            for (int j=0;j<n;j++){
                double x = cut->GetPointX(j);
                double y = cut->GetPointY(j);
                double y_new = aparam[0]*exp(aparam[2]*y)+aparam[1];
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