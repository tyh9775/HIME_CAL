#include "myConst3.h"
#include <TH2F.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>

struct lfit {
    int layer;
    double y_int;
	double y_int_err;
    double slope;
	double slope_err;
    double y_int2;
	double y_int_err2;
    double slope2;
	double slope_err2;	
};

std::vector<lfit> lin_par_load(const std::string &linp) {
    std::ifstream infile(linp);
    std::string line;
    std::vector<lfit> results;
    lfit current;
	bool v2=false;

    std::regex layer_regex(R"(Layer\s+(\d+))");
    std::regex values_regex(R"(y-int:\s*([-\d\.eE]+)\+/-([-\d\.eE]+),\s*slope:\s*([-\d\.eE]+)\+/-([-\d\.eE]+))");


    while (std::getline(infile, line)) {
        std::smatch match;
        if (std::regex_search(line, match, layer_regex)) {
            current = lfit();  // reset
            current.layer = std::stoi(match[1]);
            v2 = false;
        } else if (line.find("v2") != std::string::npos) {
            v2 = true;
        } else if (std::regex_search(line, match, values_regex)) {
            if (v2) {
                current.y_int2 = std::stod(match[1]);
                current.y_int_err2 = std::stod(match[2]);
                current.slope2 = std::stod(match[3]);
                current.slope_err2 = std::stod(match[4]);
                results.push_back(current);
            } else {
                current.y_int = std::stod(match[1]);
                current.y_int_err = std::stod(match[2]);
                current.slope = std::stod(match[3]);
                current.slope_err = std::stod(match[4]);
            }
        }
    }
    return results;
}

double E_to_tot(double E, std::vector<double> parA, std::vector<double> lp){
    double amp {(E-lp[0])/lp[1]};
    double tot {log((amp-parA[1])/parA[0])/parA[2]};
    return tot;
}

void muon_pts(
    const std::string &input = "input.root",
    const std::string &output = "output.root",
    const std::vector<std::vector<double>> parA_mat = {{a1,a2,ka},{a1,a2,ka},{a1,a2,ka}},
    const std::string &linp = "lin_par.txt"
){
    TFile *f = new TFile(input.c_str(),"READ");
    TFile *fo1 = new TFile(output.c_str(),"RECREATE");

    std::vector<int> hbases {hbase1,hbase2,hbase3};

    std::vector<lfit> lin_par_list = lin_par_load(linp);

	std::vector<std::vector<double>> lin_par {};
	std::vector<std::vector<double>> lin_par2 {};

	for (size_t i=0;i<lin_par_list.size();i++){
		std::vector<double> lp {lin_par_list[i].y_int,lin_par_list[i].slope};
		std::vector<double> lp2 {lin_par_list[i].y_int2,lin_par_list[i].slope2};
		lin_par.push_back(lp);
		lin_par2.push_back(lp2);
	}    

    for (int i=0;i<3;i++){
        TCanvas *ci = new TCanvas(Form("c%d",i),Form("Det%d",hbases[i]));
        TH2D *h = (TH2D*)f->Get(Form("hTotVsTDiff_mod_%d",hbases[i]))->Clone();
        ci->cd();
        h->Draw();
        std::vector<double> parA = parA_mat[i];
        std::vector<double> lp = lin_par[i];

        double tot = E_to_tot(8.0,parA,lp);
        double xmin = h->GetXaxis()->GetXmin();
        double xmax = h->GetXaxis()->GetXmax();

        ci->cd();
        TLine *line = new TLine(xmin,tot,xmax,tot);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw("SAME");

        fo1->cd();
        ci->Write();
    }


    fo1->Close();

    return;
}