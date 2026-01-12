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

double amp_calc(const double *par, double ToT);

//my calibration function/////////////////////////////////////////////
double f_cal(const double *param_b, double ToT,const std::vector<double> &parA);
//////////////////////////////////////////////////////////////////////

std::pair<std::vector<int>, std::vector<double*>> param_load(const std::string &pfile,int hbase,int hstart,int hstop);
//std::pair<std::vector<int>, std::vector<double*>> param_load(const std::string &pfile,int hbase);

void peak_find(TH1D *h, const std::vector<double> &yc_list);

void amp_tot( //simplified version that doesn't save individual detectors
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
    const std::string &outputFilename = "output.root",
	const std::string &pfile = "parameters.txt",
    bool useWalkCorrection = true,
    bool useVetoWall = true,
	const std::vector<std::vector<double>> &parA_mat = {{a1,a2,ka},{a1,a2,ka},{a1,a2,ka}},

	// const std::vector<double> walkParams = {-0.20, -6.389},
    // const std::vector<double> walkParams = {0.0042, -0.305, -5.863},
    //const std::vector<double> walkParams = {7.8, 0.043, -13.549},
	const std::vector<std::vector<double>> walkParams = {{7.8, 0.043, -13.549}},
	const int wc_f_select = 0, //0 - exponential, 1 - polynomial
	//const double totThresh = 22.13,
	const double totThresh = my_tot_thresh, //21.233
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
	//auto hTotVsTof = new TH2D("hTotVsTof", "", 2000, -100, 300, 150, 0, 75);
	auto hTotVsTof = new TH2D("hTotVsTof", "", 2000, -100, 300, 750, 0, 75); //my copy, finer binning in Tot
	auto hHitPattern = new TH2D("hHitPattern", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern = new TH2D("hBarHitPattern", "", hitPatternBins, -1000, 1000, 72, 0, 72);

    auto hAmpVsTof_L0 = new TH2D("hAmpVsTof_L0", "", 2000, -100, 300, 12000, 0, 6000);
    auto hAmpVsTof_L1 = new TH2D("hAmpVsTof_L1", "", 2000, -100, 300, 12000, 0, 6000);
    auto hAmpVsTof_L2 = new TH2D("hAmpVsTof_L2", "", 2000, -100, 300, 12000, 0, 6000);

    auto hAmpVsTof_L0_cal = new TH2D("hAmpVsTof_L0_cal", "", 2000, -100, 300, 12000, 0, 6000);
    auto hAmpVsTof_L1_cal = new TH2D("hAmpVsTof_L1_cal", "", 2000, -100, 300, 12000, 0, 6000);
    auto hAmpVsTof_L2_cal = new TH2D("hAmpVsTof_L2_cal", "", 2000, -100, 300, 12000, 0, 6000);

    hAmpVsTof_L0->Sumw2();
    hAmpVsTof_L1->Sumw2();
    hAmpVsTof_L2->Sumw2();
    hAmpVsTof_L0_cal->Sumw2();
    hAmpVsTof_L1_cal->Sumw2();
    hAmpVsTof_L2_cal->Sumw2();


	TH2D *hHitPatternLayer[3];
	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i] =
			new TH2D(Form("hHitPatternLayer%d", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer[i]->Sumw2();
	}


	hBdcXY->Sumw2();
	hBdcXYCut->Sumw2();
	hTotVsTof->Sumw2();
	hHitPattern->Sumw2();

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
	/* for (int i=0;i<3;i++){
		std::cout<<"layer "<<i<<"\n";
		for (size_t j=0;j<cal_det_l[i].size();j++){
			std::cout<<"det "<<cal_det_l[i][j]<<" param "<<par_l[i][j][0] <<", "<<par_l[i][j][1] <<", "<<par_l[i][j][2] <<"\n";
		}
	} */
	//my histograms
	//no shift
	TH2D *hTotVsTof_L0 = new TH2D("hTotVsTof_L0", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_L1 = new TH2D("hTotVsTof_L1", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_L2 = new TH2D("hTotVsTof_L2", "", 2000, -100, 300, 750, 0, 75);

	//shifted in Tot only
	TH2D *hTotVsTof_sh_t_L0 = new TH2D("hTotVsTof_sh_t_L0", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_sh_t_L1 = new TH2D("hTotVsTof_sh_t_L1", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_sh_t_L2 = new TH2D("hTotVsTof_sh_t_L2", "", 2000, -100, 300, 750, 0, 75);

	//hTotVsTof_sh_t_all->Sumw2();
    hTotVsTof_sh_t_L0->Sumw2();
	hTotVsTof_sh_t_L1->Sumw2();
	hTotVsTof_sh_t_L2->Sumw2();

	TH2D *hTotVsTof_my_wc_all = new TH2D("hTotVsTof_my_wc_all", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_my_wc_L0 = new TH2D("hTotVsTof_my_wc_L0", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_my_wc_L1 = new TH2D("hTotVsTof_my_wc_L1", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_my_wc_L2 = new TH2D("hTotVsTof_my_wc_L2", "", 2000, -100, 300, 750, 0, 75);

	hTotVsTof_my_wc_all->Sumw2();
    hTotVsTof_my_wc_L0->Sumw2();
	hTotVsTof_my_wc_L1->Sumw2();
	hTotVsTof_my_wc_L2->Sumw2();


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
			auto tot = std::sqrt(spirit.hime_tot0[ihit] * spirit.hime_tot1[ihit]);
            // change the tot (calibration)
            auto tot_new = tot;

			auto cur_det = std::find(cal_det_l[layerId].begin(),cal_det_l[layerId].end(),moduleId);
			double *parb = nullptr;
			if (cur_det != cal_det_l[layerId].end()){
				int det_ind = std::distance(cal_det_l[layerId].begin(),cur_det);
				double* par_b = par_l[layerId][det_ind];
				tot_new = f_cal(par_b,tot,parA);
				parb=par_b;
			}


			/* for (int ci=0;ci<3;ci++){
			    if (layerId == ci){
			        auto cur_det = std::find(cal_det_l[ci].begin(),cal_det_l[ci].end(),moduleId);
			        if (cur_det != cal_det_l[ci].end()){
                        int det_ind = std::distance(cal_det_l[ci].begin(),cur_det);
			            double* par_b = par_l[ci][det_ind];
			            tot_new = f_cal(par_b,tot);								        
			        }
			    }
			
			} */
            

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

			double *par = new double[3];
			par[0]=parA[0];
			par[1]=parA[1];
			par[2]=parA[2];
			double Amp = amp_calc(par,tot);
			double Amp_c = amp_calc(par,tot_new);

			mod_tmp.push_back(moduleId);

			hTotVsTof->Fill(tof, tot);

			hTotVsTof_my_wc_all->Fill(tof_mw,tot_new);

			if (layerId == 0){
    			hTotVsTof_L0->Fill(tof,tot);
			    hTotVsTof_sh_t_L0->Fill(tof,tot_new);
				hTotVsTof_my_wc_L0->Fill(tof_mw,tot_new);
				hAmpVsTof_L0->Fill(tof,Amp);
				hAmpVsTof_L0_cal->Fill(tof,Amp_c);
			} else if (layerId == 1){
    		    hTotVsTof_L1->Fill(tof,tot);
			    hTotVsTof_sh_t_L1->Fill(tof,tot_new);
				hTotVsTof_my_wc_L1->Fill(tof_mw,tot_new);
				hAmpVsTof_L1->Fill(tof,Amp);
				hAmpVsTof_L1_cal->Fill(tof,Amp_c);
			} else if (layerId == 2){
		        hTotVsTof_L2->Fill(tof,tot);
			    hTotVsTof_sh_t_L2->Fill(tof,tot_new);
				hTotVsTof_my_wc_L2->Fill(tof_mw,tot_new);
				hAmpVsTof_L2->Fill(tof,Amp);
				hAmpVsTof_L2_cal->Fill(tof,Amp_c);
			}
			hHitPattern->Fill(xHit, yHit);
			hBarHitPattern->Fill(layerId == 0 || layerId == 2 ? xHit : yHit, moduleId);
			hHitPatternLayer[layerId]->Fill(xHit, yHit);



		}

    }


	hBdcXY->Scale(1. / spiritEntries);
	hTotVsTof->Scale(1. / spiritEntries);
	hHitPattern->Scale(1. / spiritEntries);
	hBarHitPattern->Scale(1. / spiritEntries);


	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i]->Scale(1. / spiritEntries);

	}

    
    //for checking stats
    TList *p_og_l0 =  new TList();
    p_og_l0->SetName("yProj_L0");
    TList *p_og_l1 =  new TList();
    p_og_l1->SetName("yProj_L1");
    TList *p_og_l2 =  new TList();
    p_og_l2->SetName("yProj_L2");
    std::vector<TList*> p_og_l = {p_og_l0,p_og_l1,p_og_l2};

    TList *p_cal_l0 =  new TList();
    p_cal_l0->SetName("yProj_cal_L0");
    TList *p_cal_l1 =  new TList();
    p_cal_l1->SetName("yProj_cal_L1");
    TList *p_cal_l2 =  new TList();
    p_cal_l2->SetName("yProj_cal_L2");
    std::vector<TList*> p_cal_l = {p_cal_l0,p_cal_l1,p_cal_l2};

    if (!useVetoWall){
        std::vector<double> xcals {xc1,xc2,xc3,xc4};
        std::vector<TH2D*> hTotVsTof_L = {hTotVsTof_L0,hTotVsTof_L1,hTotVsTof_L2};
        std::vector<TH2D*> hTotVsTof_L_cal = {hTotVsTof_sh_t_L0,hTotVsTof_sh_t_L1,hTotVsTof_sh_t_L2};
        
        for (int i=0;i<3;i++){
            std::vector<std::vector<double>> yc_list {{18.0,27.0,32.0},{20.0,28.0,32.0},{20.0,28.0,30.0},{18.0,23.0,28.0}};

            for (size_t xci=0;xci<xcals.size();xci++){
                int xj = hTotVsTof_L[i]->GetXaxis()->FindBin(xcals[xci]);
                TH1D* p_og = hTotVsTof_L[i]->ProjectionY(Form("p_og_L%d_%zu",i,xci),xj,xj);
                p_og->SetTitle(Form("Y Proj of Layer %d at x=%f",i,xcals[xci]));
                TH1D* p_cal = hTotVsTof_L_cal[i]->ProjectionY(Form("p_cal_L%d_%zu",i,xci),xj,xj);
                p_cal->SetTitle(Form("Y Proj of Layer %d at x=%f (Calibrated)",i,xcals[xci]));
                
                peak_find(p_og,yc_list[xci]);
                peak_find(p_cal,yc_list[xci]);

                p_og_l[i]->Add(p_og);
                p_cal_l[i]->Add(p_cal);
            }
        }
    }

	auto ofile = new TFile(outputFilename.c_str(), "recreate");
	hBdcXY->Write();
	hBdcXYCut->Write();
	hTotVsTof->Write();
	hHitPattern->Write();
	hBarHitPattern->Write();
	
	for (int i=0;i<3;i++){
		hHitPatternLayer[i]->Write();
	}

	
	hTotVsTof_L0->Scale(1. / spiritEntries);
	hTotVsTof_L1->Scale(1. / spiritEntries);
	hTotVsTof_L2->Scale(1. / spiritEntries);

	hTotVsTof_L0->Write();
	hTotVsTof_L1->Write();
	hTotVsTof_L2->Write();

    hTotVsTof_sh_t_L0->Scale(1. / spiritEntries);
	hTotVsTof_sh_t_L1->Scale(1. / spiritEntries);
	hTotVsTof_sh_t_L2->Scale(1. / spiritEntries);

	hTotVsTof_sh_t_L0->Write();
	hTotVsTof_sh_t_L1->Write();
	hTotVsTof_sh_t_L2->Write();


    hAmpVsTof_L0->Scale(1. / spiritEntries); 
    hAmpVsTof_L1->Scale(1. / spiritEntries); 
    hAmpVsTof_L2->Scale(1. / spiritEntries); 
    hAmpVsTof_L0_cal->Scale(1. / spiritEntries); 
    hAmpVsTof_L1_cal->Scale(1. / spiritEntries); 
    hAmpVsTof_L2_cal->Scale(1. / spiritEntries); 

    hAmpVsTof_L0->Write(); 
    hAmpVsTof_L1->Write(); 
    hAmpVsTof_L2->Write();
    hAmpVsTof_L0_cal->Write();
    hAmpVsTof_L1_cal->Write();
    hAmpVsTof_L2_cal->Write();	



	if (useWalkCorrection && walkParams.size()>1){

		hTotVsTof_my_wc_all->Scale(1. / spiritEntries);
		hTotVsTof_my_wc_L0->Scale(1. / spiritEntries);
		hTotVsTof_my_wc_L1->Scale(1. / spiritEntries);
		hTotVsTof_my_wc_L2->Scale(1. / spiritEntries);

		hTotVsTof_my_wc_all->Write();
		hTotVsTof_my_wc_L0->Write();
		hTotVsTof_my_wc_L1->Write();
		hTotVsTof_my_wc_L2->Write();
	}
	



    if (!useVetoWall){
        for (int i=0;i<3;i++){
            p_og_l[i]->Write(Form("yProj_L%d",i),TObject::kSingleKey);
            p_cal_l[i]->Write(Form("yProj_cal_L%d",i),TObject::kSingleKey);
        }
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
