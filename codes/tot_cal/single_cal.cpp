#include "myConst3.h"
#include <TFile.h>
#include <TGraph.h>
#include <TCanvas.h>

std::vector<double> cut_range(const TCutG *c){
    double xmax = c->GetPointX(0);
    double xmin = c->GetPointX(0);
    double ymax = c->GetPointY(0);
    double ymin = c->GetPointY(0);
    for (int i=1;i<c->GetN();i++){
        double xval = c->GetPointX(i);
        if (xval<xmin){
            xmin=xval;
        } else if (xval>xmax){
            xmax=xval;
        }
        double yval = c->GetPointY(i);
        if (yval<ymin){
            ymin=yval;
        } else if (yval>ymax){
            ymax=yval;
        }        
    }
    return {xmin,xmax,ymin,ymax};
}

TH2D *htmp_make(const TH2D *h,const TCutG *cut){
    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();
    TH2D *htmp = (TH2D*)h->Clone();
    for (int xi=1;xi<=nx;xi++){
        double xval = h->GetXaxis()->GetBinCenter(xi);
        for (int yi=1;yi<=ny;yi++){
            double yval = h->GetYaxis()->GetBinCenter(yi);
            if (!cut->IsInside(xval,yval)){
                htmp->SetBinContent(xi,yi,0);
                htmp->SetBinError(xi,yi,0);
            }
        }   
    }

    htmp->SetEntries(htmp->Integral(1,nx,1,ny));

    return htmp;
}


std::vector<std::vector<double>> data_load(const TH2D *h,const TCutG *cut){
    std::vector<double> x {};
    std::vector<double> y {};
    std::vector<double> err {};
    std::vector<double> counts {};
    
    TH2D *htmp = htmp_make(h,cut);

    int xmin_i = htmp->GetXaxis()->FindFixBin(cut_range(cut)[0]);
    int xmax_i = htmp->GetXaxis()->FindFixBin(cut_range(cut)[1]);
    int ymin_i = htmp->GetYaxis()->FindFixBin(cut_range(cut)[2]);
    int ymax_i = htmp->GetYaxis()->FindFixBin(cut_range(cut)[3]);

    int nybins = htmp->GetNbinsY();
    for (int xj=xmin_i;xj<xmax_i;xj++){
        TH1D *py_tmp = htmp->ProjectionY(Form("py%d",xj),xj-1,xj+1);
        double xval=htmp->GetXaxis()->GetBinCenter(xj);
        for (int yj=1;yj<=nybins;yj++){
            double yval=htmp->GetYaxis()->GetBinCenter(yj);
            if (!cut->IsInside(xval,yval)){
                py_tmp->SetBinContent(yj,0);
                py_tmp->SetBinError(yj,0);
            }
        }
        
        py_tmp->SetEntries(py_tmp->Integral(ymin_i,ymax_i));
        
        if (py_tmp->GetEntries()==0){
            continue;
        }
        int max_ybin = py_tmp->GetMaximumBin();
        double max_y = py_tmp->GetBinCenter(max_ybin);
        double max_count = py_tmp->GetBinContent(max_ybin);

        int y_low = py_tmp->FindFirstBinAbove(0);
        int y_high = py_tmp->FindLastBinAbove(0);
        double y_l = py_tmp->GetBinCenter(y_low);
        double y_h = py_tmp->GetBinCenter(y_high);

        TF1 *f_1 = new TF1("f_1","gaus",y_l,y_h);
        f_1->SetParameters(max_count,max_y,1);
        f_1->SetParLimits(1,y_l,y_h);
        f_1->SetParLimits(2,0.1,5.0);
        py_tmp->Fit(f_1,"QRB+");
        double max_val = f_1->GetParameter(0);
        double max_y_fit = f_1->GetParameter(1);
        double y_err = f_1->GetParError(1);
        //y_err=y_err/(max_val*num_ent);

        if (/* abs(y_err)<=25 && */ cut->IsInside(xval,max_y_fit)){
            x.push_back(xval);
            y.push_back(max_y_fit);
            err.push_back(y_err);
            counts.push_back(max_val);
        } 

    }

    return {x,y,err,counts};
}


std::vector<double> max_find(const TGraphErrors *g){
    int n = g->GetN();
    double x_max {0.0};
    double y_max {0.0};
    double err_max {0.0};

    for (int i=0;i<n;i++){
        double xi = g->GetPointX(i);
        double yi = g->GetPointY(i);
        double err_i = g->GetErrorY(i);
        /* if (err_i>2.5 || err_i<0.0){
            continue;
        } */
        if (yi>=y_max){
            y_max=yi;
            x_max=xi;
            err_max=err_i;
        }
    }
    
    return {x_max,y_max,err_max};
}


TH2D *cust_axis(TH2D *h, const TF1 *f){
    int nx = h->GetNbinsX();
    //int ny = h->GetNbinsY();
    int ny=310;
    std::vector<double> new_edges(ny+1);
    for (int j=0;j<=ny;j++){
        double yE = h->GetYaxis()->GetBinLowEdge(j+1);
        new_edges[j]=f->Eval(yE);
    }

    std::sort(new_edges.begin(),new_edges.end());

    TH2D *hnew = new TH2D("hnew",h->GetTitle(),nx,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),ny,new_edges.data());
    
    for (int i=1;i<=nx;i++) {
        for (int j=1;j<=ny;j++) {
            double x = h->GetXaxis()->GetBinCenter(i);
            double y = h->GetYaxis()->GetBinCenter(j);
            double z = h->GetBinContent(i, j);

            hnew->Fill(x, f->Eval(y), z);
        }
    }

    return hnew;
}

std::vector<double> mu_read(const std::string &tfile,int det){
    std::ifstream infile(tfile);
    std::string line;
    int di {0};
    double tot {0.0};
    double err {0.0};

    std::getline(infile,line);

    while (std::getline(infile,line)){
        if (line.empty()) {
            continue;
        }
        
        std::istringstream ss(line);
        std::string pm;
        ss>>di>>tot>>pm>>err;

        if (di==det){
            return {tot,err};
        }
    }
    return {};
}

TH1D *ht_1D(const TH1D *py,const TCutG *cut,const double xval){
    TH1D *py_tmp = (TH1D*)py->Clone();
    int n = py->GetNbinsX();
    for (int i=1;i<=n;i++){
        double yval = py->GetBinCenter(i);
        if (!cut->IsInside(xval,yval)){
            py_tmp->SetBinContent(i,0);
            py_tmp->SetBinError(i,0);
        }
    }
    py_tmp->SetEntries(py_tmp->Integral(1,n));
    return py_tmp;
}

double avg_find(std::vector<double> v){
    double sum {0};
    for (size_t i=0;i<v.size();i++){
        sum+=v[i];
    }
    return sum/v.size();
}

void single_cal_v3(
    const std::string &file = "input.root",
    const std::string &tfile = "input.txt",    
    const std::string &ofile = "output",
    const std::vector<std::string> &cuts = {"cut1.root","cut2.root","cut3.root"},
    bool limit = true,
    bool He = false,
    const double E_min = En_min //set negative if ignoring
){
    TFile *f = new TFile(file.c_str(),"READ");
    std::string oroot;
    oroot = ofile+".root";
    TFile *fo1 = new TFile(oroot.c_str(),"RECREATE");

    ofstream output;
    std::string output_file = ofile+".txt";
    output.open(output_file.c_str(),std::ios_base::trunc);
    std::vector<double> mu {};
    std::vector<double> mu_err {};
    std::vector<int> hbases {hbase1, hbase2, hbase3};    

    for (int i=0;i<3;i++){
        std::vector<double> mu_i=mu_read(tfile,hbases[i]);
        //std::cout<<mu_i<<"\n";
        mu.push_back(mu_i[0]);
        mu_err.push_back(mu_i[1]);
    }

    std::vector<std::vector<double>> tofs {L0_tof,L1_tof,L2_tof};

    for (int i=0;i<3;i++){
        std::vector<double> tot {mu[i]};
        std::vector<double> tot2 {mu[i]};
        std::vector<double> pnch_pts {cos_mu, pro_pt, deu_pt, tri_pt};
        int num_sec = 3;

        if (E_min>=0.0){
            pnch_pts.insert(pnch_pts.begin(),E_min);
        }

        if (He){
            pnch_pts.push_back(hel_pt);
            pnch_pts.push_back(alp_pt);
            num_sec=5;
        }

        std::vector<double> tot_err {mu_err[i]};
        std::vector<double> tot2_err {mu_err[i]};

        TH2D *h = (TH2D*)f->Get(Form("hTotVsTof_det%d",hbases[i]));
        TFile *fcut = new TFile(cuts[i].c_str(),"READ");
        std::vector<TGraphErrors*> graphs {};
        std::vector<TGraph*> mxpts {};
        std::vector<TGraphErrors*> mxpts2 {};

        if (E_min>=0.0){
            TH1D *py = (TH1D*)h->ProjectionY("py");
            int y_min_bin = py->FindFirstBinAbove(0);
            double tot_min = py->GetXaxis()->GetBinCenter(y_min_bin);
            double tot_min_err = py->GetXaxis()->GetBinWidth(y_min_bin);
            tot.insert(tot.begin(),tot_min);
            tot2.insert(tot2.begin(),tot_min);
            tot_err.insert(tot_err.begin(),tot_min_err);
            tot2_err.insert(tot2_err.begin(),tot_min_err);
        }
        for (int j=1;j<=num_sec;j++){
            TCutG *cut = (TCutG*)fcut->Get(Form("cut_s%d",j));
            std::vector<std::vector<double>> data = data_load(h,cut);
            std::vector<double> x_j {data[0]};
            std::vector<double> y_j {data[1]};
            std::vector<double> err_j {data[2]};
            std::vector<double> counts_j {data[3]};

            TGraphErrors *g_j = new TGraphErrors(x_j.size(),x_j.data(),y_j.data(),nullptr,err_j.data());
            g_j->SetName(Form("g_s%d_L%d",j,i));
            double tof_j = max_find(g_j)[0];
            double tot_j = max_find(g_j)[1];
            double e_j = max_find(g_j)[2];
            tot.push_back(tot_j);
            //tot_err.push_back(e_j);
            tot_err.push_back(avg_find(err_j));
            TGraphErrors *gc = (TGraphErrors*)g_j->Clone();
            graphs.push_back(gc);
            TGraph *mxpt = new TGraph(1,&tof_j,&tot_j);
            mxpt->SetName(Form("mxpt_s%d_L%d",j,i));
            mxpt->SetMarkerStyle(20);
            mxpt->SetMarkerColor(kRed);
            mxpts.push_back(mxpt);

            /* double tot2_j = tot_j+1.0;
            tot2.push_back(tot2_j);
            TGraph *mxpt2 = new TGraph(1,&tof_j,&tot2_j);
            mxpt2->SetName(Form("mxpt2_s%d_L%d",j+1,i));
            mxpt2->SetMarkerStyle(20);
            mxpt2->SetMarkerColor(kBlue);
            mxpts2.push_back(mxpt2);    */         

            double tof_tmp = tofs[i][j-1];
            int tof_tmp_bin = h->GetXaxis()->FindBin(tof_tmp);
            TH1D *py = (TH1D*)h->ProjectionY(Form("py_L%d_s%d",i,j),tof_tmp_bin-1,tof_tmp_bin+1);
            py->SetTitle(Form("Y Proj at ToF %.1f",tof_tmp));
            fo1->cd();
            py->Write();

            TH1D *py_cut = ht_1D(py,cut,tof_tmp);
            py_cut->SetName(Form("py_L%d_s%d_cut",i,j));
            py_cut->SetTitle(Form("Y Proj at ToF %.1f (cut)",tof_tmp));
            fo1->cd();
            py_cut->Write();

            int max_ybin = py_cut->GetMaximumBin();
            double max_y = py_cut->GetBinCenter(max_ybin);
            double max_count = py_cut->GetBinContent(max_ybin);

            int y_low = py_cut->FindFirstBinAbove(0);
            int y_high = py_cut->FindLastBinAbove(0);
            /* double y_l = cut_range(cut)[2];
            double y_h = cut_range(cut)[3]; */
            double y_l = py_cut->GetBinCenter(y_low);
            double y_h = py_cut->GetBinCenter(y_high);

            TF1 *f_1 = new TF1("f_1","gaus",y_l,y_h);
            f_1->SetParameters(max_count,max_y,1);
            f_1->SetParLimits(1,y_l,y_h);
            f_1->SetParLimits(2,0.1,5.0);
            py_cut->Fit(f_1,"QRB");
            double tot2_j = f_1->GetParameter(1);
            double tot2_err_j = f_1->GetParError(1);
            tot2.push_back(tot2_j);
            tot2_err.push_back(tot2_err_j);
            TGraphErrors *mxpt2 = new TGraphErrors(1,&tof_tmp,&tot2_j,nullptr,&tot2_err_j);
            mxpt2->SetName(Form("mxpt2_s%d_L%d",j,i));
            mxpt2->SetMarkerStyle(20);
            mxpt2->SetMarkerColor(kBlue);
            mxpts2.push_back(mxpt2);

            int tof_j_bin = h->GetXaxis()->FindBin(tof_j);
            TH1D *py2 = (TH1D*)h->ProjectionY(Form("py_L%d_s%d_f",i,j),tof_j_bin-1,tof_j_bin+1);
            py2->SetTitle(Form("Y Proj at ToF %.1f",tof_j));
            fo1->cd();
            py2->Write();

            TH1D *py2_cut = ht_1D(py2,cut,tof_j);
            py2_cut->SetName(Form("py_L%d_s%d_f_cut",i,j));
            py2_cut->SetTitle(Form("Y Proj at ToF %.1f (cut)",tof_j));
            fo1->cd();
            py2_cut->Write();

        }



        TCanvas *c0 = new TCanvas(Form("hTotVsTof_det%d",hbases[i]));
        c0->cd();
        h->Draw(); 

        for (size_t g=0; g< graphs.size(); g++) {
            graphs[g]->Draw("SAME");
            mxpts[g]->Draw("P SAME");
            mxpts2[g]->Draw("P SAME");

        }

        fo1->cd();
        c0->Write();

        TGraphErrors *g_E = new TGraphErrors(tot.size(),tot.data(),pnch_pts.data(),tot_err.data(),nullptr);
        g_E->SetName(Form("g_E_L%d",i));
        g_E->GetXaxis()->SetLimits(-5,35);
        g_E->GetYaxis()->SetRangeUser(-10,200);
        g_E->SetMarkerStyle(20);
        TF1 *f_E = new TF1(Form("f_E_L%d",i),"[0]*([1]*exp([2]*x)+[3])",tot[0],tot[tot.size()-1]);
        f_E->SetParameters(0.1,0.1,0.1,0.1);
        if (limit){
            f_E->SetParLimits(0,0.0,10.0);
            //f_E->FixParameter(0,1.0);
            f_E->SetParLimits(1,0.0,1.0);
            f_E->SetParLimits(2,0.0,1.0);
            f_E->SetParLimits(3,0.0,1.0);
            //f_E->SetParLimits(4,-5.0,5.0);
        }

        g_E->Fit(f_E,"R");

        fo1->cd();
        g_E->Write();

        output<<"det "<<hbases[i]<<"\n";
        output<<"a1: "<<f_E->GetParameter(1)<<"+/-"<<f_E->GetParError(1)<<", a2: "<<f_E->GetParameter(3)<<"+/-"<<f_E->GetParError(3)<<", ka: "<<f_E->GetParameter(2)<<"+/-"<<f_E->GetParError(2)<<"\n";
        output<<"slope: "<<f_E->GetParameter(0)<<"+/-"<<f_E->GetParError(0)<</* ", int: "<<f_E->GetParameter(4)<<"+/-"<<f_E->GetParError(4) <<*/"\n";
        output<<"\n";

        TCanvas *ccopy = new TCanvas(Form("hcopy_%d",hbases[i]));
        TH2D *hcopy = (TH2D*)h->Clone();
        TH2D *hnew = cust_axis(hcopy,f_E);
        hnew->SetName(Form("hcopy%d",hbases[i]));
        ccopy->cd();
        hnew->Draw();

        fo1->cd();
        ccopy->Write();


        TGraphErrors *g2_E = new TGraphErrors(tot2.size(),tot2.data(),pnch_pts.data(),tot2_err.data(),nullptr);
        g2_E->SetName(Form("g2_E_L%d",i));
        g2_E->GetXaxis()->SetLimits(-5,35);
        g2_E->GetYaxis()->SetRangeUser(-10,200);
        g2_E->SetMarkerStyle(20);
        TF1 *f2_E = new TF1(Form("f2_E_L%d",i),"[0]*([1]*exp([2]*x)+[3])",tot2[0],tot2[tot2.size()-1]);
        f2_E->SetParameters(0.1,0.1,0.1,0.1);
        if (limit){
            f2_E->SetParLimits(0,0.0,10.0);
            //f2_E->FixParameter(0,1.0);
            f2_E->SetParLimits(1,0.0,1.0);
            f2_E->SetParLimits(2,0.0,1.0);
            f2_E->SetParLimits(3,0.0,1.0);
            //f2_E->SetParLimits(4,-5.0,5.0);
        }
        g2_E->Fit(f2_E,"R");

        fo1->cd();
        g2_E->Write();

        output<<"det "<<hbases[i]<<" version 2"<<"\n";
        output<<"a1: "<<f2_E->GetParameter(1)<<"+/-"<<f2_E->GetParError(1)<<", a2: "<<f2_E->GetParameter(3)<<"+/-"<<f2_E->GetParError(3)<<", ka: "<<f2_E->GetParameter(2)<<"+/-"<<f2_E->GetParError(2)<<"\n";
        output<<"slope: "<<f2_E->GetParameter(0)<<"+/-"<<f2_E->GetParError(0)<</* ", int: "<<f2_E->GetParameter(4)<<"+/-"<<f2_E->GetParError(4)<< */"\n";
        output<<"\n";        

        TCanvas *ccopy2 = new TCanvas(Form("hcopy2_%d",hbases[i]));
        TH2D *hcopy2 = (TH2D*)h->Clone();
        TH2D *hnew2 = cust_axis(hcopy2,f2_E);
        hnew2->SetName(Form("hcopy%d_v2",hbases[i]));
        ccopy2->cd();
        hnew2->Draw();

        fo1->cd();
        ccopy2->Write();

    }

    fo1->Close();

    return;
}