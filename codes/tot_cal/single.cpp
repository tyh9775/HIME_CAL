#include "../spirit.h"
#include "../myConst.h"
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



void single(
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
	const std::vector<int> mod_check = {14, 38, 62},
    bool useWalkCorrection = true,
    bool useVetoWall = true,

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
    const std::string& tofOffsetFile = "codes/time_offsets_experiment_gaus.txt",
    const double timeSbtLeft = -110.8,

    const int hitPatternBins = 200,
    //const std::string &spiritDir="spirit",
    const std::string &spiritDir="hime_riken_bdc",
    
	//const std::string &tdiffOffsetFile = "hime/database/calibration/tdiff/tdiff_offset.json",
	//const std::string &velocityFile = "hime/database/velocity_marco.txt",


	const std::string &tdiffOffsetFile = "codes/position-calib-muon.json",


	const std::string &modulePositionFilename = "codes/positions.txt"
    //const std::string &modulePositionFilename = "hime/database/module_positions.dat"

	// clang-format on
) {

	auto spiritChain = getSpiritChain(runIds, spiritDir);
	auto spiritEntries = spiritChain->GetEntries();
	std::cout << "Total entries: " << spiritEntries << "\n";

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

	auto hBdcXY = new TH2D("hBdcXY", "", 200, -100, 100, 200, -100, 100);
	auto hBdcXYCut = new TH2D("hBdcXYCut", "", 200, -100, 100, 200, -100, 100);
	hBdcXY->Sumw2();
	hBdcXYCut->Sumw2();

	/****************************************************************************************************************/

	std::vector<TH2D*> hTotVsTof_det {};
	std::vector<std::vector<TH2D*>> hAmpVsTof_det_par {};
	for (size_t i=0;i<mod_check.size();i++){
		TH2D *hTotVsTof_i = new TH2D(Form("hTotVsTof_det%d",mod_check[i]), "", 2000, -100, 300, 750, 0, 75);
		hTotVsTof_i->Sumw2();
		hTotVsTof_det.push_back(hTotVsTof_i);

	}

	// for generating random positions
	static std::random_device rd;
	static std::mt19937 gen(rd());
	static std::uniform_real_distribution<> unif(0.0, 1.0);

	auto barWidth = barWidthInMeter * 1000; // mm
	/****************************************************************************************************************/

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
			
		for (auto ihit = 0; ihit < spirit.hime_nHits ; ihit++) {

			int moduleId = spirit.hime_moduleID[ihit];
			int layerId = getLayer(moduleId);

			if (std::find(mod_check.begin(), mod_check.end(), moduleId) == mod_check.end()) {
				continue;
			}

			if (layerId < 0 || layerId > 2) {
				std::cerr << "Invalid layerId: " << layerId << std::endl;
				continue;
			}

			if (std::find(requiredLayers.begin(), requiredLayers.end(), layerId) == requiredLayers.end()) {
				continue;
			}
			auto tot = std::sqrt(spirit.hime_tot0[ihit] * spirit.hime_tot1[ihit]);
                       

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

			double tot_new=tot;

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

			int mod_pos = std::find(mod_check.begin(), mod_check.end(), moduleId) - mod_check.begin();

			hTotVsTof_det[mod_pos]->Fill(tof, tot);
			
		}

    }

	for (size_t i=0;i<mod_check.size();i++){
		hTotVsTof_det[i]->Scale(1. / spiritEntries);
	}

	auto ofile = new TFile(outputFilename.c_str(), "recreate");
	std::vector<std::vector<double>> pt_tofs {L0_tof, L1_tof, L2_tof};
	std::vector<std::vector<double>> pt_tots {L0_tot, L1_tot, L2_tot};

	for (size_t i=0;i<mod_check.size();i++){
		ofile->cd();
		hTotVsTof_det[i]->Write();
		TGraph *g1 = new TGraph(pt_tofs[i].size(), pt_tofs[i].data(), pt_tots[i].data());
		TCanvas *c1 = new TCanvas(Form("c_tot_tof_det%d",mod_check[i]), "");
		c1->cd();
		gPad->SetLogz();
		hTotVsTof_det[i]->Draw();
		g1->SetMarkerStyle(20);
		g1->SetMarkerColor(kRed);
		g1->Draw("P SAME");
		c1->Write();
	}
	ofile->Close();

	return;
}

std::vector<double> getOffset(const std::string &filename) {
	std::vector<double> offset(nModules, 0);
	std::ifstream infile(filename.c_str());
	if (!infile.is_open()) {
		std::cerr << "Could not open file " << filename << std::endl;
		return offset;
	}

	std::string line;
	while (std::getline(infile, line)) {
		if (line[0] == '#') {
			continue;
		}
		std::stringstream ss(line);
		int detid;
		ss >> detid >> offset[detid];
	}
	return offset;
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



