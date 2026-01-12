#include "spirit.h"
#include "myConst3.h"
#include <TSpectrum.h>


const int nLayers = 3;
const int nModules = 72;
const double barWidthInMeter = 0.04; // 40 mm
const double C_LIGHT = 0.299792458;	 // m/ns

const double FRONT_DISTANCE = 4.8;	// m
const double LAYER_DISTANCE = 0.06; // m
const double barDepth = 2e-2;		// 2cm

unsigned int getLayer(const unsigned int &hitModule);

void getPositionCalibration(const std::string &filename, std::vector<double> &veff, std::vector<double> &offset);


std::vector<double> getModulePositions(const std::string &filename);

void loadTofOffset(const std::string &filename, std::vector<double> &offset);

//my calibration function/////////////////////////////////////////////
double amp_calc(const double *par, double ToT);

double f_cal(const double *param_b, double ToT, const std::vector<double> &parA);
//////////////////////////////////////////////////////////////////////

std::pair<std::vector<int>, std::vector<double*>> param_load(const std::string &pfile,int hbase,int hstart,int hstop);
//std::pair<std::vector<int>, std::vector<double*>> param_load(const std::string &pfile,int hbase);

void peak_find(TH1D *h, const std::vector<double> &yc_list);

//my multiplicity scaling function
double func(const double x, const double y, const double a, const double b, const double c);

int mult_count(const TH1D *h);

void shadow_mult_v2(
	// clang-format off
    const std::vector<int> runIds = {

        // forward angle shadow bar
        1148, 1149, 1150, 1151, 1152, 1153, 1154, 1155, 1156, 1158, 1160, 1161, 1162, 1163, 1164, 1165,
        1166, 1167, 1168, 1169, 1170, 1171, 1172, 1174, 1175, 1176, 1178, 1179, 1180,
        1181, 1182, 1183, 1184, 1185, 1186, 1189, 1190, 1191, 1192, 1195, 1196,
        1197, 1198, 1199, 1200, 1201, 1202, 1203, 1204, 1205, 1206, 1207, 1208, 1210, 
        
        // top forward-angled bar moved to bottom backward-angle
        1211, 1212, 1213, 1214, 1215, 1217, 1218, 1219, 1221, 1222, 1223, 1225, 1226, 1227, 1232, 1234, 1235, 1236, 1239, 1240, 1241, 1242, 1243, 1244, 1245, 
        
        // only bottom backward-angled bar
        1247, 1248, 1251, 1252, 1253, 1254, 1255, 1256, 1261, 1262, 1263, 1264, 1265, 1266, 1268, 1269, 1273, 1274, 1275, 1276, 1277, 1279, 1280, 1281, 1282, 1283, 1284, 1285, 1286, 1287, 1288, 1289, 1290, 1291, 1292, 1293, 1294, 1295, 1296, 1297, 1298, 1299,

        // remaining 
        // 1300, 1304, 1305, 1306, 1308, 1309, 1310, 1311, 1312, 1313, 1314, 1315, 1316, 1317, 1318, 1319, 1320, 1321, 1322, 1324, 1325, 1326, 1327, 1331,
    }, 
    const std::string &outputFile = "output1",
	const std::string &pfile = "parameters.txt",
    bool useWalkCorrection = true,
    bool useVetoWall = true,
	const std::vector<std::vector<double>> &parA_mat = {{a1,a2,ka},{a1,a2,ka},{a1,a2,ka}},
	const std::vector<std::vector<double>> &par_mult_mat = {{1.0,1.0,0.0},{1.0,1.0,0.0},{1.0,1.0,0.0}},

	// const std::vector<double> walkParams = {-0.20, -6.389},
    // const std::vector<double> walkParams = {0.0042, -0.305, -5.863},
    //const std::vector<double> walkParams = {7.8, 0.043, -13.549},
	const std::vector<std::vector<double>> walkParams = {{7.8, 0.043, -13.549}},
	const int wc_f_select = 0, //0 - exponential, 1 - polynomial
	//const double totThresh = 22.13,
	const double totThresh = my_tot_thresh_wc, //21.233
	bool cutVertex = true,

	const std::vector<double>& tofRange = {-1000,1000},

    std::vector<double> xRange = {-1000, 1000},
    std::vector<double> yRange = {-1000, 1000},

    const std::vector<int>& requiredLayers = {0,1,2},
    // const std::string& tofOffsetFile = "./hime/database/calibration/tof/tof_offset.json",
    const std::string& tofOffsetFile = "hime/database/calibration/tof/time_offsets_experiment_gaus.txt",
    const double timeSbtLeft = -110.8,

    const int hitPatternBins = 200,
    //const std::string &spiritDir="spirit",
    const std::string &spiritDir="hime_riken_bdc",
    
	//const std::string &tdiffOffsetFile = "hime/database/calibration/tdiff/tdiff_offset.json",
	//const std::string &velocityFile = "hime/database/velocity_marco.txt",


	const std::string &tdiffOffsetFile = "hime/database/calibration/tdiff/position-calib-muon.json",
	
	
	const std::string &modulePositionFilename = "hime/database/test.txt"
    //const std::string &modulePositionFilename = "hime/database/module_positions.dat"

	// clang-format on
) {

	auto spiritChain = getSpiritChain(runIds, spiritDir);
	auto spiritEntries = spiritChain->GetEntries();

	// load tdiff offset and velocity
	std::vector<double> veff(nModules, 0.0);
	std::vector<double> posOffset(nModules, 0.0);
	getPositionCalibration(tdiffOffsetFile, veff, posOffset);

	auto modulePositions = getModulePositions(modulePositionFilename);
	for (auto i = 0; i < nModules; i++) {
		modulePositions[i] *= 1000; // m -> mm
	}

	// load tof offset
	std::vector<double> tofOffset;
	loadTofOffset(tofOffsetFile, tofOffset); 

	/****************************************************************************************************************/

	auto hBdcXY = new TH2D("hBdcXY", "", 200, -100, 100, 200, -100, 100);
	auto hBdcXYCut = new TH2D("hBdcXYCut", "", 200, -100, 100, 200, -100, 100);
	auto hHitPattern = new TH2D("hHitPattern", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern = new TH2D("hBarHitPattern", "", hitPatternBins, -1000, 1000, 72, 0, 72);


	TH2D *hHitPatternLayer[3];
	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i] =
			new TH2D(Form("hHitPatternLayer%d", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer[i]->Sumw2();
	}
	TH1D *hBarHit[72];


	
	hBdcXY->Sumw2();
	hBdcXYCut->Sumw2();
	hHitPattern->Sumw2();
	hBarHitPattern->Sumw2();

	// for generating random positions
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> unif(0.0, 1.0);

	auto barWidth = barWidthInMeter * 1000; // mm
	/****************************************************************************************************************/


	//load my parameters
	std::pair<std::vector<int>, std::vector<double*>> par1 = param_load(pfile,hbase1,hstart1,hstop1);
	std::pair<std::vector<int>, std::vector<double*>> par2 = param_load(pfile,hbase2,hstart2,hstop2);
	std::pair<std::vector<int>, std::vector<double*>> par3 = param_load(pfile,hbase3,hstart3,hstop3);
	std::vector<std::vector<int>> cal_det_l = {par1.first,par2.first,par3.first};
	std::vector<std::vector<double*>> par_l = {par1.second,par2.second,par3.second};
	std::vector<int> hb_l = {hbase1,hbase2,hbase3};

	//my histograms

	TH1D *hMult_all = new TH1D("hMult_all","Multiplicity (All)",20,0,20);
    TH1D *hMult_L[3];
    for (int i = 0; i < 3; i++) {
        hMult_L[i] = new TH1D(Form("hMult_L%d", i), Form("Multiplicity(Layer %d)", i), 20, 0, 20);
        hMult_L[i]->Sumw2();
	}

	hMult_all->Sumw2();



	TH1D *hHits_all = new TH1D("hHits_all","Hits per Event (All)",20,0,20);
    TH1D *hHits_L[3];
    for (int i = 0; i < 3; i++) {
        hHits_L[i] = new TH1D(Form("hHits_L%d", i), Form("Hits per Event (Layer %d)", i), 20, 0, 20);
        hHits_L[i]->Sumw2();
	}
	hHits_all->Sumw2();


	TH2D *hHits_vs_Mult = new TH2D("hHits_vs_Mult","Hits vs Multiplicity",20,0,20,20,0,20);
	hHits_vs_Mult->Sumw2();
	TH2D *hHits_vs_Mult_L[3];
	for (int i = 0; i < 3; i++) {
		hHits_vs_Mult_L[i] = new TH2D(Form("hHits_vs_Mult_L%d", i), Form("Hits vs Multiplicity (Layer %d)", i),20,0,20,20,0,20);
		hHits_vs_Mult_L[i]->Sumw2();
	}

	TH2D *hHitPattern_test = (TH2D*)hHitPattern->Clone("hHitPattern_test"); //apply f(x,y) for comparison
	TH2D *hHitPattern_test_layer[3];
	for (int i = 0; i < 3; i++) {
		hHitPattern_test_layer[i] = (TH2D*)hHitPatternLayer[i]->Clone(Form("hHitPattern_test_L%d", i));
	}

    //layer exclusive events (hit only one layer within the event)
	auto hHitPattern_exc = new TH2D("hHitPattern_exc", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern_exc = new TH2D("hBarHitPattern_exc", "", hitPatternBins, -1000, 1000, 72, 0, 72);
	TH2D *hHitPatternLayer_exc[3];


	for (int i = 0; i < 3; i++) {
		hHitPatternLayer_exc[i] =
			new TH2D(Form("hHitPatternLayer%d_exc", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer_exc[i]->Sumw2();
	}
	TH1D *hBarHit_exc[72];

	hHitPattern_exc->Sumw2();
	hBarHitPattern_exc->Sumw2();
		

	for (auto iEvt = 0; iEvt < spiritEntries; iEvt++) {
		spiritChain->GetEntry(iEvt);
		hBdcXY->Fill(spirit.tbdc_x, spirit.tbdc_y);
	}
	auto xminBin = hBdcXY->GetXaxis()->FindBin(-20), xmaxBin = hBdcXY->GetXaxis()->FindBin(20);
	auto yminBin = hBdcXY->GetYaxis()->FindBin(-20), ymaxBin = hBdcXY->GetYaxis()->FindBin(20);
	auto hBdcX = (TH1D *)hBdcXY->ProjectionX("hBdcX", yminBin, ymaxBin);
	auto hBdcY = (TH1D *)hBdcXY->ProjectionY("hBdcY", xminBin, xmaxBin);
	auto bdcXmean = hBdcX->GetMean();
	auto bdcYmean = hBdcY->GetMean();

	auto bdcXsigma = hBdcX->GetRMS();
	auto bdcYsigma = hBdcY->GetRMS();
	std::cout << "BdcX mean: " << bdcXmean << " sigma: " << bdcXsigma << std::endl;
	std::cout << "BdcY mean: " << bdcYmean << " sigma: " << bdcYsigma << std::endl;

	// auto bdcXmean = 2.124;
	// auto bdcYmean = -5.569;
	// auto bdcXsigma = 3.154;
	// auto bdcYsigma = 3.154;

	/****************************************************************************************************************/

	for (auto iEvt = 0; iEvt < spiritEntries; iEvt++) {
		spiritChain->GetEntry(iEvt);

		auto vertexCut = spirit.tbdc_x > (bdcXmean - 2 * bdcXsigma) && spirit.tbdc_x < (bdcXmean + 2 * bdcXsigma) &&
						 spirit.tbdc_y > (bdcYmean - 2 * bdcYsigma) && spirit.tbdc_y < (bdcYmean + 2 * bdcYsigma);
		if (!vertexCut && cutVertex) {
			continue;
		}

		hBdcXYCut->Fill(spirit.tbdc_x, spirit.tbdc_y);
        
		auto vetoEvent = useVetoWall ? spirit.hime_veto_multi > 0 : false;
		if (vetoEvent) {
			continue;
		}

		hHits_all->Fill(spirit.hime_nHits);

		double sbtTime = 0.;
		double tpcTime = 0.;

		for (auto ihit = 0; ihit < spirit.hime_nHits; ihit++) {
			int moduleId = spirit.hime_moduleID[ihit];
			if (moduleId == 1) {
				sbtTime = spirit.hime_tofRaw[ihit];
			}

			// module 0 = channel 0 and 8 = tpc and veto
			if (moduleId == 0) {
				tpcTime = spirit.hime_tofRaw[ihit];
			}
		}
		if (sbtTime == 0.) {
			continue;
		}

		std::vector<int> layer_tmp {};
		std::vector<double> tof_tmp {};
		std::vector<double> tot_tmp {};
		std::vector<double> tot_c_tmp {};

		
		std::vector<double> x_tmp {};
		std::vector<double> y_tmp {};

		std::vector<int> mod_tmp {};

		int layer_mult[3] = {0,0,0};
		int n_hits[3] = {0,0,0};
		int full_mult = 0;

		for (auto ihit = 0; ihit < spirit.hime_nHits ; ihit++) {

			int moduleId = spirit.hime_moduleID[ihit];
			int layerId = getLayer(moduleId);
			
			std::vector<double> parA=parA_mat[layerId];
			
			if (layerId < 0 || layerId > 2) {
				std::cerr << "Invalid layerId: " << layerId << std::endl;
				continue;
			}

			if (std::find(requiredLayers.begin(), requiredLayers.end(), layerId) == requiredLayers.end()) {
				continue;
			}

			n_hits[layerId]++;

			auto tot = std::sqrt(spirit.hime_tot0[ihit] * spirit.hime_tot1[ihit]);
            // change the tot (calibration)
            auto tot_new = tot;

			auto cur_det = std::find(cal_det_l[layerId].begin(),cal_det_l[layerId].end(),moduleId);
			if (cur_det != cal_det_l[layerId].end()){
				int det_ind = std::distance(cal_det_l[layerId].begin(),cur_det);
				double* par_b = par_l[layerId][det_ind];
				tot_new = f_cal(par_b,tot,parA);								        
			}
           

			//auto tDiff = spirit.hime_tDiff[ihit] + posOffset[moduleId] / velocity[moduleId];
			auto tDiff = spirit.hime_tDiff[ihit];


			if (moduleId <= 1 || moduleId >= nModules) {
				continue;
			}

			//auto barPosition = tDiff * velocity[moduleId] * 10; // cm -> mm

			auto barPosition = 0.5 * veff[moduleId] * tDiff + posOffset[moduleId];
			
			double xHit = 0, yHit = 0;

			/* // extra position offset due to the module we used to gate the tdiff
			auto gatedMod = gatedModules[layerId][1];
			auto extra_offset = modulePositions[gatedMod]; */

			if (layerId == 0 || layerId == 2) {
				//xHit = barPosition + extra_offset;
				xHit = barPosition;
				yHit = modulePositions[moduleId] + barWidth * (unif(gen) - 0.5);
			} else {
				//yHit = barPosition + extra_offset;
				yHit = barPosition;
				xHit = modulePositions[moduleId] + barWidth * (unif(gen) - 0.5);
			}

			if (layerId == 0 || layerId == 2) {
				if (xHit < xRange[0] || xHit > xRange[1]) {
					continue;
				}
			}

			if (layerId == 1) {
				if (yHit < yRange[0] || yHit > yRange[1]) {
					continue;
				}
			}


			auto rHit = std::sqrt(xHit * xHit + yHit * yHit) * 1e-3; // mm -> m
			auto distanceCenter = FRONT_DISTANCE + LAYER_DISTANCE * layerId;
			auto distanceHit = std::sqrt(distanceCenter * distanceCenter + rHit * rHit);
			auto tof = spirit.hime_tofRaw[ihit] - (distanceHit - FRONT_DISTANCE) / C_LIGHT;

			// should not matter whether an extra common offset is added or not
			// added for consistency with Zibi's analysis
			tof -= timeSbtLeft - sbtTime;

			tof -= tofOffset[moduleId];

			// some constant offset applied in Zibi's analysis
			// to offset to 0, should not matter
			// tof-offset calibration only aligns prompt gamma
			tof += (11.8 + 2. * 4.8 / C_LIGHT);
			
			//correcting tof w/ new tot
			auto tof_new = tof;

			auto tof_mw = tof;

			if (useWalkCorrection && walkParams.size()==1) {
				if (tot < totThresh) {
					// tof -= (walkParams[0] * tot + walkParams[1]);
					// tof -= (walkParams[0] * tot * tot + walkParams[1] * tot + walkParams[2]);
					tof -= (walkParams[0][0] * std::exp(-walkParams[0][1] * tot) + walkParams[0][2]);
					tof_new -= (walkParams[0][0] * std::exp(-walkParams[0][1] * tot_new) + walkParams[0][2]);

				} else {
					// tof -= (walkParams[0] * totThresh + walkParams[1]);
					// tof -= (walkParams[0] * totThresh * totThresh + walkParams[1] * totThresh + walkParams[2]);
					tof -= (walkParams[0][0] * std::exp(-walkParams[0][1] * totThresh) + walkParams[0][2]);
					tof_new -= (walkParams[0][0] * std::exp(-walkParams[0][1] * totThresh) + walkParams[0][2]);
				}
			} else if (useWalkCorrection && walkParams.size()>1) {
				if (tot < totThresh) {
					// tof -= (walkParams[0] * tot + walkParams[1]);
					// tof -= (walkParams[0] * tot * tot + walkParams[1] * tot + walkParams[2]);
					if (wc_f_select==0){
						tof += walkParams[layerId][0]*exp(walkParams[layerId][1]*tot)+walkParams[layerId][2];
						tof_new += walkParams[layerId][0]*exp(walkParams[layerId][1]*tot_new)+walkParams[layerId][2];
						tof_mw += walkParams[3][0]*exp(walkParams[3][1]*tot_new)+walkParams[3][2];
					} else if (wc_f_select==1){
						for (size_t par_i=0;par_i<walkParams[layerId].size();par_i++){
							tof += walkParams[layerId][par_i]*pow(tot,par_i);
							tof_new += walkParams[layerId][par_i]*pow(tot_new,par_i);
							tof_mw += walkParams[3][par_i]*pow(tot_new,par_i);
						}
					}
				} else {
					// tof -= (walkParams[0] * totThresh + walkParams[1]);
					// tof -= (walkParams[0] * totThresh * totThresh + walkParams[1] * totThresh + walkParams[2]);
					if (wc_f_select==0){
						tof += walkParams[layerId][0]*exp(walkParams[layerId][1]*totThresh)+walkParams[layerId][2];
						tof_new += walkParams[layerId][0]*exp(walkParams[layerId][1]*totThresh)+walkParams[layerId][2];
						tof_mw += walkParams[3][0]*exp(walkParams[3][1]*totThresh)+walkParams[3][2];
					} else if (wc_f_select==1){
						for (size_t par_i=0;par_i<walkParams[layerId].size();par_i++){
							tof += walkParams[layerId][par_i]*pow(totThresh,par_i);
							tof_new += walkParams[layerId][par_i]*pow(totThresh,par_i);
							tof_mw += walkParams[3][par_i]*pow(totThresh,par_i);
						}
					}
				}
			}


			if (tof < tofRange[0] || tof > tofRange[1]) {
				continue;
			}
			if (tof_new < tofRange[0] || tof_new > tofRange[1]) {
				continue;
			}

			layer_mult[layerId]++;
			full_mult++;

			tot_c_tmp.push_back(tot_new);
			layer_tmp.push_back(layerId);
			tof_tmp.push_back(tof);
			tot_tmp.push_back(tot);
			
			x_tmp.push_back(xHit);
			y_tmp.push_back(yHit);

			mod_tmp.push_back(moduleId);

			double scl_test = func(xHit,yHit,par_mult_mat[layerId][0],par_mult_mat[layerId][1],par_mult_mat[layerId][2]);
			double prob_test = unif(gen);


			hHitPattern->Fill(xHit, yHit);
			hBarHitPattern->Fill(layerId == 0 || layerId == 2 ? xHit : yHit, moduleId);
			hHitPatternLayer[layerId]->Fill(xHit, yHit);
			if (prob_test<scl_test){
				hHitPattern_test->Fill(xHit, yHit);
				hHitPattern_test_layer[layerId]->Fill(xHit, yHit);
			}

		}
		hMult_all->Fill(full_mult);
		hMult_L[0]->Fill(layer_mult[0]);
		hMult_L[1]->Fill(layer_mult[1]);
		hMult_L[2]->Fill(layer_mult[2]);

		hHits_L[0]->Fill(n_hits[0]);
		hHits_L[1]->Fill(n_hits[1]);
		hHits_L[2]->Fill(n_hits[2]);

		hHits_vs_Mult->Fill(full_mult,spirit.hime_nHits);
		hHits_vs_Mult_L[0]->Fill(layer_mult[0],n_hits[0]);
		hHits_vs_Mult_L[1]->Fill(layer_mult[1],n_hits[1]);
		hHits_vs_Mult_L[2]->Fill(layer_mult[2],n_hits[2]);

		for (int tar_tmp=0;tar_tmp<3;tar_tmp++){
            if (!layer_tmp.empty() && std::all_of(layer_tmp.begin(),layer_tmp.end(),[tar_tmp](int n) {return n==tar_tmp;})){

				for (size_t ti=0;ti<tof_tmp.size();ti++){
					hHitPattern_exc->Fill(x_tmp[ti], y_tmp[ti]);
					hBarHitPattern_exc->Fill(layer_tmp[ti]== 0 || layer_tmp[ti]== 2 ? x_tmp[ti] : y_tmp[ti], mod_tmp[ti]);
					hHitPatternLayer_exc[layer_tmp[ti]]->Fill(x_tmp[ti], y_tmp[ti]);
			    }
		    } else{
				continue;
			}
            
        }

    }

	std::ofstream omult;
	omult.open((outputFile+".txt").c_str());
	omult<<std::setw(13)<<" "<<std::setw(20)<<" Including Background "<<std::setw(20)<<" Without Background "<<"\n";

	hBdcXY->Scale(1. / spiritEntries);

	std::cout<<"Total Multiplicity: "<<mult_count(hMult_all)<<"\n";
	for (int i = 0; i < 3; i++) {
		std::cout<<"Layer "<<i<<" Multiplicity: "<<mult_count(hMult_L[i])<<"\n";
	}

	double tot_mult = hHitPattern->Integral();
	std::cout<<"Including background:"<<"\n";
	std::cout<<"Total multiplicity: "<<tot_mult<<"\n";


	hHitPattern->Scale(1. / spiritEntries);
	std::cout<<"Average multiplicity: "<<hHitPattern->Integral()<<"\n";

	hBarHitPattern->Scale(1. / spiritEntries);

	hHitPattern_exc->Scale(1. / spiritEntries);
	hBarHitPattern_exc->Scale(1. / spiritEntries);

	
	double tot_mult_test = hHitPattern_test->Integral();

	std::cout<<"Without background:"<<"\n";
	std::cout<<"Total multiplicity: "<<tot_mult_test<<"\n";


	hHitPattern_test->Scale(1. / spiritEntries);

	std::cout<<"Average multiplicity: "<<hHitPattern_test->Integral()<<"\n";

	for (int i = 0; i < 3; i++) {
		double layer_mult = hHitPatternLayer[i]->Integral();
		double layer_mult_test = hHitPattern_test_layer[i]->Integral();
		omult<<"Layer "<<i<<"\n";
		omult<<std::setw(10)<<"Avg Mult: "<<std::setw(20)<<layer_mult/spiritEntries<<std::setw(20)<<layer_mult_test/spiritEntries<<"\n";
		omult<<std::setw(10)<<"Tot Mult: "<<std::setw(20)<<layer_mult<<std::setw(20)<<layer_mult_test<<"\n";

		hHitPatternLayer[i]->Scale(1. / spiritEntries);

		hHitPatternLayer_exc[i]->Scale(1. / spiritEntries);

		hHitPattern_test_layer[i]->Scale(1. / spiritEntries);

	}

	omult<<"All "<<"\n";
	omult<<std::setw(10)<<"Avg Mult: "<<std::setw(20)<<tot_mult/spiritEntries<<std::setw(20)<<tot_mult_test/spiritEntries<<"\n";
	omult<<std::setw(10)<<"Tot Mult: "<<std::setw(20)<<tot_mult<<std::setw(20)<<tot_mult_test<<"\n";

	for (int i = 0; i < 72; i++) {
		hBarHit[i] = (TH1D *)hBarHitPattern->ProjectionX(Form("hBarHit%d", i), i + 1, i + 1);
		hBarHit_exc[i] = (TH1D *)hBarHitPattern_exc->ProjectionX(Form("hBarHit_exc%d", i), i + 1, i + 1);		
	}

	TFile *ofile = new TFile((outputFile+".root").c_str(), "recreate");




	ofile->cd();

	hBdcXY->Scale(1.0,"width");
	hBdcXYCut->Scale(1.0,"width");
	hHitPattern->Scale(1.0,"width");
	hBarHitPattern->Scale(1.0,"width");
	hHitPattern_exc->Scale(1.0,"width");
	hBarHitPattern_exc->Scale(1.0,"width");

	hBdcXY->Write();
	hBdcXYCut->Write();
	hHitPattern->Write();
	hBarHitPattern->Write();


	for (int i=0;i<3;i++){
		hHitPatternLayer[i]->Scale(1.0,"width");
		hHitPatternLayer_exc[i]->Scale(1.0,"width");
		hHitPatternLayer[i]->Write();
		hHitPatternLayer_exc[i]->Write();
    }
	hHitPattern_exc->Write();
	hBarHitPattern_exc->Write();

	hHitPattern_test->Scale(1.0,"width");
	hHitPattern_test->Write();
	for (int i=0;i<3;i++){
		hHitPattern_test_layer[i]->Scale(1.0,"width");
		hHitPattern_test_layer[i]->Write();
	}

	hMult_all->Write();
	hMult_L[0]->Write();
	hMult_L[1]->Write();
	hMult_L[2]->Write();

	hHits_all->Write();
	hHits_L[0]->Write();
	hHits_L[1]->Write();
	hHits_L[2]->Write();

	hHits_vs_Mult->Write();
	hHits_vs_Mult_L[0]->Write();
	hHits_vs_Mult_L[1]->Write();
	hHits_vs_Mult_L[2]->Write();


	for (int i = 0; i < 72; i++) {
		hBarHit[i]->Scale(1.0,"width");
		hBarHit_exc[i]->Scale(1.0,"width");
		hBarHit[i]->Write();
		hBarHit_exc[i]->Write();

	}

	ofile->Close();


	return;
}



unsigned int getLayer(const unsigned int &hitModule) {
	if (hitModule >= 0 && hitModule <= 23) {
		return 0;
	}
	if (hitModule >= 24 && hitModule <= 47) {
		return 1;
	}
	if (hitModule >= 48 && hitModule <= 71) {
		return 2;
	}
	return -1;
}



/////////////////////////////
void getPositionCalibration(const std::string &filename, std::vector<double> &veff, std::vector<double> &offset) {
	std::ifstream infile(filename.c_str());

	assert(veff.size() == nModules);
	assert(offset.size() == nModules);

	for (int i = 0; i < nModules; i++) {
		veff[i] = 0.0;
		offset[i] = 0.0;
	}

	nlohmann::json j;
	infile >> j;
	infile.close();

	for (int i = 0; i < nModules; i++) {
		if (j.contains(std::to_string(i))) {
			if (!j[std::to_string(i)].contains("v") || !j[std::to_string(i)].contains("offset")) {
				continue; // skip if v or offset is not present
			}
			veff[i] = j[std::to_string(i)]["v"].get<double>();
			offset[i] = j[std::to_string(i)]["offset"].get<double>();
		} else {
			std::cerr << "Module " << i << " not found in calibration file." << std::endl;
		}
	}
}
////////////////////////////////




std::vector<double> getModulePositions(const std::string &filename) {
	std::vector<double> positions(nModules, 0);
	std::ifstream infile(filename.c_str());
	if (!infile.is_open()) {
		std::cerr << "Could not open file " << filename << std::endl;
		return positions;
	}

	infile.ignore(1000, '\n');
	//infile.ignore(1000, '\n');

	int id, layer;
	double pos;

	while (infile >> id >> layer >> pos) {
		positions[id] = pos;
	}
	return positions;
}

void loadTofOffset(const std::string &filename, std::vector<double> &offset) {
	offset.clear();
	offset.resize(nModules, 0);

	if (filename == "") {
		return;
	}

	if (!std::filesystem::exists(filename)) {
		std::cerr << "Could not open file " << filename << std::endl;
		return;
	}

	std::string extension = filename.substr(filename.find_last_of(".") + 1);

	if (extension == "json") {
		std::ifstream infile(filename.c_str());
		// load json file
		nlohmann::json j;
		infile >> j;

		for (auto i = 0; i < nModules; i++) {
			offset[i] = j[std::to_string(i)];
		}

	} else if (extension == "txt") {
		std::ifstream infile(filename.c_str());
		int id;
		double off;
		while (infile >> id >> off) {
			offset[id] = off;
		}
	}
	return;
}


//my calibration function/////////////////////////////////////////////
double amp_calc(const double *par, double ToT){
	double A {par[0]};
    double B {par[1]};
    double k {par[2]};

    return A*exp(k*ToT)+B;
}

//my calibration function/////////////////////////////////////////////
double f_cal(const double *param_b, double ToT, const std::vector<double> &parA){
    double Amp = amp_calc(param_b,ToT);
	double A = parA[0];
	double B = parA[1];
	double k = parA[2];
    double ToT_cal {log((Amp-B)/A)/k};
    return ToT_cal;
}
//////////////////////////////////////////////////////////////////////



std::pair<std::vector<int>, std::vector<double*>> param_load(const std::string &pfile,int hbase,int hstart, int hstop) {
	std::ifstream file(pfile);
	std::string line;
	std::vector<int> cal_det {}; //detectors that I calibrated
	std::vector<double*> parameters {}; //parameters of the detectors

	while (std::getline(file,line)){

		if (line.find("baseline") != std::string::npos || line.find("Det Num") != std::string::npos) {
            continue;
        }
	    std::istringstream iss(line);
		int det {};
		double b1 {};
		double b1_err {};
		double b2 {};
		double b2_err {};
		double kb {};
		double kb_err {};
	    

		/* if (!(iss >> det >> b1 >> std::ws && iss.ignore(10, '+') && iss.ignore(10, '/') && iss >> b1_err
                  >> b2 >> std::ws && iss.ignore(10, '+') && iss.ignore(10, '/') && iss >> b2_err
                  >> kb  >> std::ws && iss.ignore(10, '+') && iss.ignore(10, '/') && iss >> kb_err)) {
            continue;
        } */
		std::string dummy;
		if (!(iss >> det >> b1 >> dummy >> b1_err >> b2 >> dummy >> b2_err >> kb >> dummy >> kb_err)) {
			continue;
		}
		



	    if (det<hstart || det==hbase){ 
	        continue;
	    }
		if (det>hstop){
			break;
		}
	    cal_det.push_back(det);
	    
	    double* par = new double[3];
	    
	    par[0]=b1;
		par[1]=b2;
		par[2]=kb;
	    
	    parameters.push_back(par);

	}    

	file.close();
    return std::make_pair(cal_det,parameters);
}

void peak_find(TH1D *h, const std::vector<double> &yc_list){
	std::vector<TH1D*> hc_list {};
	if (yc_list.size()==3){
	    int bin1 = h->FindBin(yc_list[0]);
	    int bin2 = h->FindBin(yc_list[1]);
	    int bin3 = h->FindBin(yc_list[2]);
	    int n_b1 = bin2-bin1+1;
	    TH1D* hc1 =  new TH1D("hcl","lower copy",n_b1,h->GetBinLowEdge(bin1),h->GetBinLowEdge(bin2+1));
	    for (int j=0;j<n_b1;j++){
	        hc1->SetBinContent(j+1,h->GetBinContent(bin1+j));
	    }
	    int n_b2 = bin3-bin2+1;
	    TH1D* hc2 =  new TH1D("hcu","upper copy",n_b2,h->GetBinLowEdge(bin2),h->GetBinLowEdge(bin3+1));
	    for (int j=0;j<n_b2;j++){
	        hc2->SetBinContent(j+1,h->GetBinContent(bin2+j));
	    }
	    
	    hc_list.push_back(hc1);
	    hc_list.push_back(hc2);
	    		    
	} else {
	    int bin1 = h->FindBin(yc_list[0]);
	    int bin2 = h->FindBin(yc_list[1]);
	    int n_b1 = bin2-bin1+1;
	    TH1D* hc1 =  new TH1D("hcl","lower copy",n_b1,h->GetBinLowEdge(bin1),h->GetBinLowEdge(bin2+1));
	    for (int j=0;j<n_b1;j++){
	        hc1->SetBinContent(j+1,h->GetBinContent(bin1+j));
	    }
	    
	    hc_list.push_back(hc1);
	}
	
	for (size_t j=0;j<hc_list.size();j++){
        TSpectrum *spec = new TSpectrum(2);
	    int n_peaks = spec->Search(hc_list[j],4,"",0.3);
	    
	    Double_t *x_peaks = spec->GetPositionX();
	    
	    for (int xi=0;xi<n_peaks;xi++){
		    Double_t xp = x_peaks[xi];
		    Double_t a0 = h->GetBinContent(h->FindBin(xp));
		    Double_t sig = 1.0;
		    
		    TF1 *f1 = new TF1(Form("gaus_%d",xi),"gaus",xp-2.5,xp+2.5);
		    f1->SetParameters(a0,xp,sig);
		    f1->FixParameter(1,xp);
		    h->Fit(f1,"RQ+");
	    }		
	}
	
	for (size_t j=0;j<hc_list.size();j++){
	    hc_list[j]->Delete();
	}
		
    return;
}


double func(const double x, const double y, const double a, const double b, const double c){
	return a*x+b*y+c;
}

int mult_count(const TH1D *h){
	int mult = 0;
	for (int i=0;i<h->GetNbinsX();i++){
		mult+= h->GetBinContent(i+1)*i;
	}
	return mult;
}