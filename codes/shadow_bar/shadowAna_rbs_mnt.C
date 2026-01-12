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

std::vector<lfit> lin_par_load(const std::string &linp);

void peak_find(TH1D *h, const std::vector<double> &yc_list);

double En_calc(const std::vector<double> &lin_par, double amp);

double n_mnt(double En);


void shadowAna_rbs_v3(
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
    const std::string &outputFilename = "output",
	const std::string &pfile = "parameters.txt",
	const std::string &linp = "lin_parameters.txt",
    bool useWalkCorrection = true,
    bool useVetoWall = true,
	const std::vector<std::vector<double>> &parA_mat = {{a1,a2,ka},{a1,a2,ka},{a1,a2,ka}},

	bool evt_norm = true,

	// const std::vector<double> walkParams = {-0.20, -6.389},
    // const std::vector<double> walkParams = {0.0042, -0.305, -5.863},
    //const std::vector<double> walkParams = {7.8, 0.043, -13.549},
	const std::vector<std::vector<double>> walkParams = {{7.8, 0.043, -13.549}},
	const int wc_f_select = 0, //0 - exponential, 1 - polynomial
	//const double totThresh = 22.13,
	const double totThresh = my_tot_thresh_wc, //21.233
	bool cutVertex = true, //true by default

	const std::vector<double>& tofRange = {-1000.,1000.},

    std::vector<double> xRange = {-1000., 1000.},
    std::vector<double> yRange = {-1000., 1000.},

    const std::vector<int>& requiredLayers = {0,1,2},
    // const std::string& tofOffsetFile = "./hime/database/calibration/tof/tof_offset.json",
    const std::string& tofOffsetFile = "hime/database/calibration/tof/time_offsets_experiment_gaus.txt",
    const double timeSbtLeft = -110.8,

    const int hitPatternBins = nbins_hitpattern, //default 200
    //const std::string &spiritDir="spirit",
    const std::string &spiritDir="hime_riken_bdc",
    
	//const std::string &tdiffOffsetFile = "hime/database/calibration/tdiff/tdiff_offset.json",
	//const std::string &velocityFile = "hime/database/velocity_marco.txt",


	const std::string &tdiffOffsetFile = "hime/database/calibration/tdiff/position-calib-muon.json",
	
	
	const std::string &modulePositionFilename = "hime/database/module_positions_new.txt"
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

	double barDepthmm = barDepth * 1000; // mm
	auto barWidth = barWidthInMeter * 1000; // mm

	auto hBdcXY = new TH2D("hBdcXY", "", 200, -100, 100, 200, -100, 100);
	auto hBdcXYCut = new TH2D("hBdcXYCut", "", 200, -100, 100, 200, -100, 100);

	TH2D *hBdcXYCut_random = new TH2D("hBdcXYCut_random", "", 200, -100, 100, 200, -100, 100);

	auto hHitPattern = new TH2D("hHitPattern", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
	auto hBarHitPattern = new TH2D("hBarHitPattern", "", hitPatternBins, -1000, 1000, 72, 0, 72);

	TH3D *hHitPattern_mnt= new TH3D("hHitPattern_mnt", "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000, 250, 0, 500);
	TH3D *hBarHitPattern_mnt = new TH3D("hBarHitPattern_mnt", "", hitPatternBins, -1000, 1000, 72, 0, 72, 250, 0, 500);

	TH2D *hHitPatternLayer[3];
	TH3D *hHitPatternLayer_mnt[3];

	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i] =
			new TH2D(Form("hHitPatternLayer%d", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000);
		hHitPatternLayer[i]->Sumw2();
		hHitPatternLayer_mnt[i] = new TH3D(
			Form("hHitPatternLayer%d_mnt", i), "", hitPatternBins, -1000, 1000, hitPatternBins, -1000, 1000, 250, 0, 500);
		hHitPatternLayer_mnt[i]->Sumw2();
	}

	TH1D *hBarHit[72];
	TH2D *hBarHit_mnt[72];

	TH1D *hBarHit_all = new TH1D("hBarHit_all", "", hitPatternBins, -1000, 1000);
	TH2D *hBarHit_all_mnt = new TH2D("hBarHit_all_mnt", "", hitPatternBins, -1000, 1000, 250, 0, 500);

	TH2D *hTotVsTofBar[72];
	for (int i = 0; i < 72; i++) {
		//hTotVsTofBar[i] = new TH2D(Form("hTotVsTofBar%d", i), "", 2000, -100, 300, 150, 0, 75);
		hTotVsTofBar[i] = new TH2D(Form("hTotVsTofBar%d", i), "", 2000, -100, 300, 750, 0, 75);
		hTotVsTofBar[i]->Sumw2();
	}

    std::vector<TH2D*> hTotVsTofBar_sh_t {}; //just tot shift
	for (int i = 0; i < 72; i++) {
		TH2D *hTotVsTofBar_i = new TH2D(Form("hTotVsTofBar_sh_t%d", i), "", 2000, -100, 300, 750, 0, 75);
		hTotVsTofBar_i->Sumw2();
		hTotVsTofBar_sh_t.push_back(hTotVsTofBar_i);
	}
	
    std::vector<TH2D*> hTotVsTofBar_sh_f {}; //also tof shift from tot shift
	for (int i = 0; i < 72; i++) {
		TH2D *hTotVsTofBar_i = new TH2D(Form("hTotVsTofBar_sh_f%d", i), "", 2000, -100, 300, 750, 0, 75);
		hTotVsTofBar_i->Sumw2();
		hTotVsTofBar_sh_f.push_back(hTotVsTofBar_i);
	}

	hBdcXY->Sumw2();

	hBdcXYCut->Sumw2();
	hBdcXYCut_random->Sumw2();

	//hTotVsTof->Sumw2();
	hHitPattern->Sumw2();
	hHitPattern_mnt->Sumw2();
	hBarHitPattern->Sumw2();
	hBarHitPattern_mnt->Sumw2();
	hBarHit_all->Sumw2();
	hBarHit_all_mnt->Sumw2();


	TH2D *hHitPattern_rnd = (TH2D*)hHitPattern->Clone("hHitPattern_rnd");
	TH2D *hBarHitPattern_rnd = (TH2D*)hBarHitPattern->Clone("hBarHitPattern_rnd");

	TH2D *hHitPatternLayer_rnd[3];
	for (int i = 0; i < 3; i++) {
		hHitPatternLayer_rnd[i] = (TH2D*)hHitPatternLayer[i]->Clone(Form("hHitPatternLayer%d_rnd", i));
	}

	TH1D *hBarHit_all_rnd = new TH1D("hBarHit_all_rnd", "", hitPatternBins, -1000, 1000);
	hBarHit_all_rnd->Sumw2();
	

	// for generating random positions
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> unif(0.0, 1.0);


	/****************************************************************************************************************/

	//load my parameters
	std::pair<std::vector<int>, std::vector<double*>> par1 = param_load(pfile,hbase1,hstart1,hstop1);
	std::pair<std::vector<int>, std::vector<double*>> par2 = param_load(pfile,hbase2,hstart2,hstop2);
	std::pair<std::vector<int>, std::vector<double*>> par3 = param_load(pfile,hbase3,hstart3,hstop3);
	std::vector<std::vector<int>> cal_det_l = {par1.first,par2.first,par3.first};
	std::vector<std::vector<double*>> par_l = {par1.second,par2.second,par3.second};
	std::vector<int> hb_l = {hbase1,hbase2,hbase3};
	
	std::vector<lfit> lin_par_list = lin_par_load(linp);
	std::vector<std::vector<double>> lin_par {};
	std::vector<std::vector<double>> lin_par2 {};
	for (size_t i=0;i<lin_par_list.size();i++){
		std::vector<double> lp {lin_par_list[i].y_int,lin_par_list[i].slope};
		std::vector<double> lp2 {lin_par_list[i].y_int2,lin_par_list[i].slope2};
		lin_par.push_back(lp);
		lin_par2.push_back(lp2);
	}

	

	//my histograms
	//no shift
	TH2D *hTotVsTof_L0 = new TH2D("hTotVsTof_L0", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_L1 = new TH2D("hTotVsTof_L1", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_L2 = new TH2D("hTotVsTof_L2", "", 2000, -100, 300, 750, 0, 75);

	//shifted in Tot only
	TH2D *hTotVsTof_sh_t_all = new TH2D("hTotVsTof_sh_t_all", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_sh_t_L0 = new TH2D("hTotVsTof_sh_t_L0", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_sh_t_L1 = new TH2D("hTotVsTof_sh_t_L1", "", 2000, -100, 300, 750, 0, 75);
	TH2D *hTotVsTof_sh_t_L2 = new TH2D("hTotVsTof_sh_t_L2", "", 2000, -100, 300, 750, 0, 75);

	hTotVsTof_sh_t_all->Sumw2();
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

    
    //layer exclusive events (hit only one layer within the event)
    std::vector<TH2D*> hTotVsTof_exc_L {};
    std::vector<TH2D*> hTotVsTof_exc_L_cal {};
    for (int i=0;i<3;i++){
	    TH2D *hTotVsTof_exc_i = new TH2D(Form("hTotVsTof_L%d_exc",i), "", 2000, -100, 300, 750, 0, 75); 
	    TH2D *hTotVsTof_exc_i_cal = new TH2D(Form("hTotVsTof_L%d_exc_cal",i), "", 2000, -100, 300, 750, 0, 75);
	    hTotVsTof_exc_i->Sumw2();
		hTotVsTof_exc_i_cal->Sumw2();
		hTotVsTof_exc_L.push_back(hTotVsTof_exc_i);
	    hTotVsTof_exc_L_cal.push_back(hTotVsTof_exc_i_cal);
    }

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

		double rand_num = unif(gen);
		if (rand_num > 0.47) { 
			hBdcXYCut_random->Fill(spirit.tbdc_x, spirit.tbdc_y);
		}


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

			auto barPosition = 0.5 * veff[moduleId] * tDiff + posOffset[moduleId];
			
			double xHit = 0, yHit = 0;

			/* // extra position offset due to the module we used to gate the tdiff
			auto gatedMod = gatedModules[layerId][1];
			auto extra_offset = modulePositions[gatedMod]; */

			if (layerId == 0 || layerId == 2) {
				xHit = barPosition;
				yHit = modulePositions[moduleId] + barWidth * (unif(gen) - 0.5);
			} else {
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

			double amp = parA[0]*exp(parA[2]*tot_new)+parA[1];
			double En = En_calc(lin_par[layerId], amp);
			double mnt = n_mnt(En);

			tot_c_tmp.push_back(tot_new);
			layer_tmp.push_back(layerId);
			tof_tmp.push_back(tof);
			tot_tmp.push_back(tot);
			
			x_tmp.push_back(xHit);
			y_tmp.push_back(yHit);

			mod_tmp.push_back(moduleId);
			hTotVsTofBar[moduleId]->Fill(tof, tot);
			hTotVsTofBar_sh_t[moduleId]->Fill(tof, tot_new);
			hTotVsTofBar_sh_f[moduleId]->Fill(tof_new, tot_new);
			hTotVsTof_sh_t_all->Fill(tof,tot_new);
			hTotVsTof_my_wc_all->Fill(tof_mw,tot_new);
			

			if (layerId == 0){
				hTotVsTof_L0->Fill(tof,tot);
				hTotVsTof_sh_t_L0->Fill(tof,tot_new);
				hTotVsTof_my_wc_L0->Fill(tof_mw,tot_new);					
		    } else if (layerId == 1){
				hTotVsTof_L1->Fill(tof,tot);
				hTotVsTof_sh_t_L1->Fill(tof,tot_new);
				hTotVsTof_my_wc_L1->Fill(tof_mw,tot_new);
			} else if (layerId == 2){
				hTotVsTof_L2->Fill(tof,tot);
				hTotVsTof_sh_t_L2->Fill(tof,tot_new);
				hTotVsTof_my_wc_L2->Fill(tof_mw,tot_new);
		    }

			hHitPattern->Fill(xHit, yHit);
			hBarHitPattern->Fill(layerId == 0 || layerId == 2 ?	 xHit : yHit, moduleId);
			hHitPatternLayer[layerId]->Fill(xHit, yHit);

			hBarHit_all->Fill(layerId == 0 || layerId == 2 ? xHit : yHit);

			hHitPattern_mnt->Fill(xHit, yHit, mnt);
			hBarHitPattern_mnt->Fill(layerId == 0 || layerId == 2 ? xHit : yHit, moduleId, mnt);
			hHitPatternLayer_mnt[layerId]->Fill(xHit, yHit, mnt);

			hBarHit_all_mnt->Fill(layerId == 0 || layerId == 2 ? xHit : yHit, mnt);
			
			if (rand_num > 0.47) {
				hHitPattern_rnd->Fill(xHit, yHit);
				hBarHitPattern_rnd->Fill(layerId == 0 || layerId == 2 ? xHit : yHit, moduleId);
				hHitPatternLayer_rnd[layerId]->Fill(xHit, yHit);

				hBarHit_all_rnd->Fill(layerId == 0 || layerId == 2 ? xHit : yHit);
			}



		}

		for (int tar_tmp=0;tar_tmp<3;tar_tmp++){
            if (!layer_tmp.empty() && std::all_of(layer_tmp.begin(),layer_tmp.end(),[tar_tmp](int n) {return n==tar_tmp;})){

				for (size_t ti=0;ti<tof_tmp.size();ti++){
				//for (int ti=0;ti<spirit.hime_nHits;ti++){
					hTotVsTof_exc_L[tar_tmp]->Fill(tof_tmp[ti],tot_tmp[ti]);
					hTotVsTof_exc_L_cal[tar_tmp]->Fill(tof_tmp[ti],tot_c_tmp[ti]); 
					hHitPattern_exc->Fill(x_tmp[ti], y_tmp[ti]);
					hBarHitPattern_exc->Fill(layer_tmp[ti]== 0 || layer_tmp[ti]== 2 ? x_tmp[ti] : y_tmp[ti], mod_tmp[ti]);
					hHitPatternLayer_exc[layer_tmp[ti]]->Fill(x_tmp[ti], y_tmp[ti]);	
				}
		    } else{
				continue;
			}
            
        }

    }

	std::cout<<spiritEntries<<"\n";

	int rnd_ent = spiritEntries;
	if (!evt_norm){
		spiritEntries = hBarHitPattern->GetEntries();
		std::cout<<"Normalized by number of hits: "<<spiritEntries<<"\n";
		rnd_ent = hBarHitPattern_rnd->GetEntries();
	}


	hBdcXY->Scale(1. / spiritEntries);
	

	hBdcXYCut->Scale(1. / spiritEntries);


	hBdcXYCut_random->Scale(1. / rnd_ent);




	//hTotVsTof->Scale(1. / spiritEntries);	

	hHitPattern->Scale(1.0,"width");
	hHitPattern_exc->Scale(1.0,"width");	
	hHitPattern->Scale(1. / spiritEntries);
	hHitPattern_exc->Scale(1. / spiritEntries);


	hHitPattern_mnt->Scale(1.0/(hHitPattern_mnt->GetYaxis()->GetBinWidth(1)*hHitPattern_mnt->GetXaxis()->GetBinWidth(1)));
	hHitPattern_mnt->Scale(1. / spiritEntries);
	

	hHitPattern_rnd->Scale(1.0,"width");
	hHitPattern_rnd->Scale(1. / rnd_ent);
	

	for (int i = 0; i < 3; i++) {
		hHitPatternLayer[i]->Scale(1.0,"width");
		hHitPatternLayer_exc[i]->Scale(1.0,"width");
		hHitPatternLayer[i]->Scale(1. / spiritEntries);
		hHitPatternLayer_exc[i]->Scale(1. / spiritEntries);
		hHitPatternLayer_rnd[i]->Scale(1.0,"width");
		hHitPatternLayer_rnd[i]->Scale(1. / rnd_ent);

		hTotVsTof_exc_L[i]->Scale(1. / spiritEntries);
		hTotVsTof_exc_L_cal[i]->Scale(1. / spiritEntries);

		hHitPatternLayer_mnt[i]->Scale(1.0/(hHitPatternLayer_mnt[i]->GetYaxis()->GetBinWidth(1)*hHitPatternLayer_mnt[i]->GetXaxis()->GetBinWidth(1)));
		hHitPatternLayer_mnt[i]->Scale(1. / spiritEntries);

	}

	TH1D *hBarHit_rb[72];
	TH1D *hBarHit_rb2[72];
	TH1D *hBarHit_rb3[72];

	TH1D *hBarHit_rnd[72];
	TH1D *hBarHit_rnd_rb[72];
	TH1D *hBarHit_rnd_rb2[72];
	TH1D *hBarHit_rnd_rb3[72];

	TH2D *hBarHit_mnt_rb[72];


	TString htitlex = "Normalized #frac{dN}{dx} vs x; x [mm]; #frac{dN}{dx} [mm^{-1}]";
	TString htitley = "Normalized #frac{dN}{dy} vs y; y [mm]; #frac{dN}{dy} [mm^{-1}]";
	TString htitleAx = "Normalized #frac{dN}{dA} vs x; x [mm]; #frac{dN}{dA} [mm^{-2}]";
	TString htitleAy = "Normalized #frac{dN}{dA} vs y; y [mm]; #frac{dN}{dA} [mm^{-2}]";

	//bin contents are still counts here, make them density later in the code
	for (int i = 0; i < 72; i++) {
		hBarHit[i] = (TH1D*)hBarHitPattern->ProjectionX(Form("hBarHit%d", i), i + 1, i + 1);
		hBarHit[i]->SetTitle(htitlex);

		if (i>=24 && i<48){
			hBarHit[i]->SetTitle(htitley);
		}	

		hBarHit_rb[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb", i));
		hBarHit_rb2[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb2", i));
		hBarHit_rb3[i] = (TH1D*)hBarHit[i]->Clone(Form("hBarHit%d_rb3", i));

		hBarHit_rb[i]->Rebin(rbf);
		hBarHit_rb2[i]->Rebin(rbf*2);
		hBarHit_rb3[i]->Rebin(rbf*4);

		hBarHit_exc[i] = (TH1D *)hBarHitPattern_exc->ProjectionX(Form("hBarHit_exc%d", i), i + 1, i + 1);


		hBarHit_rnd[i] = (TH1D*)hBarHitPattern_rnd->ProjectionX(Form("hBarHit%d_rnd", i), i + 1, i + 1);
		hBarHit_rnd[i]->SetTitle(htitlex);
		if (i>=24 && i<48){
			hBarHit_rnd[i]->SetTitle(htitley);
		}

		hBarHit_rnd_rb[i] = (TH1D*)hBarHit_rnd[i]->Clone(Form("hBarHit%d_rnd_rb", i));
		hBarHit_rnd_rb2[i] = (TH1D*)hBarHit_rnd[i]->Clone(Form("hBarHit%d_rnd_rb2", i));
		hBarHit_rnd_rb3[i] = (TH1D*)hBarHit_rnd[i]->Clone(Form("hBarHit%d_rnd_rb3", i));

		hBarHit_rnd_rb[i]->Rebin(rbf);
		hBarHit_rnd_rb2[i]->Rebin(rbf*2);
		hBarHit_rnd_rb3[i]->Rebin(rbf*4);

		hBarHitPattern_mnt->GetYaxis()->SetRange(i+1,i+1);
		hBarHit_mnt[i] = (TH2D*)hBarHitPattern_mnt->Project3D("zx");
		hBarHit_mnt[i]->SetName(Form("hBarHit_mnt%d", i));
		hBarHit_mnt[i]->SetTitle(Form("Det %d Momentum vs Position; Position [mm]; Momentum [MeV/c]", i));
		hBarHit_mnt_rb[i] = (TH2D*)hBarHit_mnt[i]->Clone(Form("hBarHit_mnt%d_rb", i));
		hBarHit_mnt_rb[i]->RebinX(2*rbf);

		hTotVsTofBar[i]->Scale(1. / spiritEntries);
		hTotVsTofBar_sh_t[i]->Scale(1. / spiritEntries);
		hTotVsTofBar_sh_f[i]->Scale(1. / spiritEntries);

	}
	
	//now make density plots
	hBarHitPattern->Scale(1.0,"width");
	hBarHitPattern_exc->Scale(1.0,"width");
	
	hBarHitPattern->Scale(1. / spiritEntries);
	hBarHitPattern_exc->Scale(1. / spiritEntries);

	hBarHitPattern_mnt->Scale(1.0/(hBarHitPattern_mnt->GetYaxis()->GetBinWidth(1)*hBarHitPattern_mnt->GetXaxis()->GetBinWidth(1)));
	hBarHitPattern_mnt->Scale(1. / spiritEntries);

	hBarHitPattern_rnd->Scale(1.0,"width");
	//hBarHitPattern_rnd->Scale(1.0,"width");	
	hBarHitPattern_rnd->Scale(1. / spiritEntries);

	hBarHit_all->SetTitle(htitlex);

	TH1D *hBarHit_all_rb = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb");
	TH1D *hBarHit_all_rb2 = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb2");
	TH1D *hBarHit_all_rb3 = (TH1D*)hBarHit_all->Clone("hBarHit_all_rb3");

	hBarHit_all_rb->Rebin(rbf);
	hBarHit_all_rb2->Rebin(rbf*2);
	hBarHit_all_rb3->Rebin(rbf*4);

	hBarHit_all->Scale(1.0,"width");
	hBarHit_all_rb->Scale(1.0,"width");
	hBarHit_all_rb2->Scale(1.0,"width");
	hBarHit_all_rb3->Scale(1.0,"width");

	hBarHit_all->Scale(1. / spiritEntries);
	hBarHit_all_rb->Scale(1. / spiritEntries);
	hBarHit_all_rb2->Scale(1. / spiritEntries);
	hBarHit_all_rb3->Scale(1. / spiritEntries);


	hBarHit_all_rnd->SetTitle("Normalized #frac{dN}{dx} vs x of All Det (random); x [mm]; #frac{dN}{dx} [mm^{-1}]");

	TH1D *hBarHit_all_rnd_rb = (TH1D*)hBarHit_all_rnd->Clone("hBarHit_all_rnd_rb");
	TH1D *hBarHit_all_rnd_rb2 = (TH1D*)hBarHit_all_rnd->Clone("hBarHit_all_rnd_rb2");
	TH1D *hBarHit_all_rnd_rb3 = (TH1D*)hBarHit_all_rnd->Clone("hBarHit_all_rnd_rb3");

	hBarHit_all_rnd_rb->Rebin(rbf);
	hBarHit_all_rnd_rb2->Rebin(rbf*2);
	hBarHit_all_rnd_rb3->Rebin(rbf*4);

	hBarHit_all_rnd->Scale(1.0,"width");
	hBarHit_all_rnd_rb->Scale(1.0,"width");
	hBarHit_all_rnd_rb2->Scale(1.0,"width");
	hBarHit_all_rnd_rb3->Scale(1.0,"width");

	hBarHit_all_rnd->Scale(1. / spiritEntries);
	hBarHit_all_rnd_rb->Scale(1. / spiritEntries);
	hBarHit_all_rnd_rb2->Scale(1. / spiritEntries);
	hBarHit_all_rnd_rb3->Scale(1. / spiritEntries);

	hBarHit_all_mnt->Scale(1.0,"width");
	hBarHit_all_mnt->Scale(1. / spiritEntries);


	auto ofile1 = new TFile((outputFilename+"-tot.root").c_str(), "recreate");
	auto ofile2 = new TFile((outputFilename+"-hit.root").c_str(), "recreate");
	auto ofile3 = new TFile((outputFilename+"-mnt.root").c_str(), "recreate");

	ofile1->cd();

    hTotVsTof_sh_t_L0->Scale(1. / spiritEntries);
	hTotVsTof_sh_t_L1->Scale(1. / spiritEntries);
	hTotVsTof_sh_t_L2->Scale(1. / spiritEntries);
	hTotVsTof_L0->Scale(1. / spiritEntries);
	hTotVsTof_L1->Scale(1. / spiritEntries);
	hTotVsTof_L2->Scale(1. / spiritEntries);

	hTotVsTof_L0->Write();
	hTotVsTof_L1->Write();
	hTotVsTof_L2->Write();
	hTotVsTof_sh_t_L0->Write();
	hTotVsTof_sh_t_L1->Write();
	hTotVsTof_sh_t_L2->Write();

    for (int i=0;i<3;i++){
        hTotVsTof_exc_L[i]->Write();
        hTotVsTof_exc_L_cal[i]->Write();
    }

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
	
	for (int i = 0; i < 72; i++) {
		hTotVsTofBar[i]->Write();
		hTotVsTofBar_sh_t[i]->Write();
	}



	ofile2->cd();
	
	hBdcXY->Write();

	hBdcXYCut->Write();

	hBdcXYCut_random->Write();


	hHitPattern->Write();
	hBarHitPattern->Write();

	for (int i=0;i<3;i++){
		hHitPatternLayer[i]->Write();
		hHitPatternLayer_exc[i]->Write();
		hHitPatternLayer_rnd[i]->Write();
		
    }
	hHitPattern_exc->Write();
	hBarHitPattern_exc->Write();

	hHitPattern_rnd->Write();
	hBarHitPattern_rnd->Write();

	hBarHit_all->Write();
	hBarHit_all_rb->Write();
	hBarHit_all_rb2->Write();
	hBarHit_all_rb3->Write();

	hBarHit_all_rnd->Write();
	hBarHit_all_rnd_rb->Write();
	hBarHit_all_rnd_rb2->Write();
	hBarHit_all_rnd_rb3->Write();


	for (int i = 0; i < 72; i++) {
		hBarHit[i]->Scale(1.0,"width");
		hBarHit_exc[i]->Scale(1.0,"width");

		hBarHit[i]->Scale(1. / spiritEntries);
		hBarHit_exc[i]->Scale(1. / spiritEntries);

		hBarHit_rb[i]->Scale(1.0,"width");
		hBarHit_rb2[i]->Scale(1.0,"width");
		hBarHit_rb3[i]->Scale(1.0,"width");

		hBarHit_rb[i]->Scale(1. / spiritEntries);
		hBarHit_rb2[i]->Scale(1. / spiritEntries);
		hBarHit_rb3[i]->Scale(1. / spiritEntries);


		hBarHit[i]->Write();
		hBarHit_rb[i]->Write();
		hBarHit_rb2[i]->Write();
		hBarHit_rb3[i]->Write();


		hBarHit_exc[i]->Write();

		hBarHit_rnd[i]->Scale(1.0,"width");

		hBarHit_rnd[i]->Scale(1. / rnd_ent);
		hBarHit_rnd[i]->Write();

		hBarHit_rnd_rb[i]->Scale(1.0,"width");
		hBarHit_rnd_rb2[i]->Scale(1.0,"width");
		hBarHit_rnd_rb3[i]->Scale(1.0,"width");

		hBarHit_rnd_rb[i]->Scale(1. / rnd_ent);
		hBarHit_rnd_rb2[i]->Scale(1. / rnd_ent);
		hBarHit_rnd_rb3[i]->Scale(1. / rnd_ent);

		hBarHit_rnd_rb[i]->Write();
		hBarHit_rnd_rb2[i]->Write();
		hBarHit_rnd_rb3[i]->Write();


	}

	ofile3->cd();

	hHitPattern_mnt->Write();
	hBarHitPattern_mnt->Write();
	hBarHit_all_mnt->Write();

	hHitPatternLayer_mnt[0]->Write();
	hHitPatternLayer_mnt[1]->Write();
	hHitPatternLayer_mnt[2]->Write();
	
	for (int i = 0; i < 72; i++) {
		double xwidth = hBarHit_mnt[i]->GetXaxis()->GetBinWidth(1);
		hBarHit_mnt[i]->Scale(1.0/xwidth);
		hBarHit_mnt[i]->Scale(1. / spiritEntries);
		hBarHit_mnt[i]->Write();

		double xwidth_rb = hBarHit_mnt_rb[i]->GetXaxis()->GetBinWidth(1);
		hBarHit_mnt_rb[i]->Scale(1.0/xwidth_rb);
		hBarHit_mnt_rb[i]->Scale(1. / spiritEntries);
		hBarHit_mnt_rb[i]->Write();
	}

	ofile1->Close();
	ofile2->Close();
	ofile3->Close();

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


double En_calc(const std::vector<double> &lin_par, double amp){
	return lin_par[0]+amp*lin_par[1];
}

double n_mnt(const double En){
	double E_total = En + neu_m;
	return pow(E_total*E_total-neu_m*neu_m, 0.5); 
}
