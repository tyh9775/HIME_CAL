#ifndef SPIRIT_HH
#define SPIRIT_HH

#include "TChain.h"
#include <array>
#include <string>
#include <vector>

constexpr int HIME_MAXHITS = 256;
constexpr int SAMURAI_MAXHITS = 256;
struct spiritData {
	unsigned int hime_nHits;
	std::array<double, HIME_MAXHITS> hime_tofRaw;
	std::array<double, HIME_MAXHITS> hime_tDiff;
	std::array<double, HIME_MAXHITS> hime_tSum;
	std::array<double, HIME_MAXHITS> hime_tot0;
	std::array<double, HIME_MAXHITS> hime_tot1;
	std::array<unsigned int, HIME_MAXHITS> hime_moduleID;

	ULong64_t hime_slowScaler;
	ULong64_t hime_fastScaler;
	ULong64_t hime_eventNumber;

	unsigned int run, event;
	ULong64_t lupots;
	unsigned int kyoto_multi;
	std::array<unsigned int, SAMURAI_MAXHITS> kyoto_bar;
	unsigned int hime_veto_multi;
	std::array<unsigned int, SAMURAI_MAXHITS> hime_veto_bar;
	std::array<double, SAMURAI_MAXHITS> hime_veto_tof;
	std::array<double, SAMURAI_MAXHITS> hime_veto_charge;
	std::array<double, SAMURAI_MAXHITS> hime_veto_tdiff;
	std::array<double, SAMURAI_MAXHITS> hime_veto_x;

	Double_t tbdc_x;
	Double_t tbdc_y;
	Double_t tbdc_a;
	Double_t tbdc_b;

	void reset();
};

spiritData spirit;

void spiritData::reset() {
	hime_nHits = 0;
	hime_tofRaw.fill(0);
	hime_tDiff.fill(0);
	hime_tSum.fill(0);
	hime_tot0.fill(0);
	hime_tot1.fill(0);
	hime_moduleID.fill(0);

	hime_slowScaler = 0;
	hime_fastScaler = 0;
	hime_eventNumber = 0;

	run = 0;
	event = 0;
	lupots = 0;
	kyoto_multi = 0;
	kyoto_bar.fill(0);
	hime_veto_multi = 0;
	hime_veto_bar.fill(0);
	hime_veto_tof.fill(0);
	hime_veto_charge.fill(0);
	hime_veto_tdiff.fill(0);
	hime_veto_x.fill(0);

	tbdc_x = 0;
	tbdc_y = 0;
	tbdc_a = 0;
	tbdc_b = 0;

	return;
}

TChain *getSpiritChain(const std::vector<int> &runIds, const std::string &spiritDir = "spirit") {
	auto chain = new TChain("spirit", "spirit");

	for (auto idx : runIds) {
		std::string fname = Form("%s/%04d.root", spiritDir.c_str(), idx);
		//std::string fname = Form("%s/data%04d.root", spiritDir.c_str(), idx);
		if (runIds.size() == 1){
			std::cout<<fname<<"\n";
		}
		chain->AddFile(fname.c_str());
	}

	chain->SetBranchAddress("hime_nHits", &spirit.hime_nHits);
	chain->SetBranchAddress("hime_tofRaw", &spirit.hime_tofRaw[0]);
	chain->SetBranchAddress("hime_tDiff", &spirit.hime_tDiff[0]);
	chain->SetBranchAddress("hime_tSum", &spirit.hime_tSum[0]);
	chain->SetBranchAddress("hime_tot0", &spirit.hime_tot0[0]);
	chain->SetBranchAddress("hime_tot1", &spirit.hime_tot1[0]);
	chain->SetBranchAddress("hime_moduleID", &spirit.hime_moduleID[0]);

	chain->SetBranchAddress("hime_slowScaler", &spirit.hime_slowScaler);
	chain->SetBranchAddress("hime_fastScaler", &spirit.hime_fastScaler);
	chain->SetBranchAddress("hime_eventNumber", &spirit.hime_eventNumber);

	chain->SetBranchAddress("runNumber", &spirit.run);
	chain->SetBranchAddress("eventNumber", &spirit.event);
	chain->SetBranchAddress("lupoTimeStamp", &spirit.lupots);
	chain->SetBranchAddress("kyotoMulti", &spirit.kyoto_multi);
	chain->SetBranchAddress("kyotoBarId", &spirit.kyoto_bar[0]);
	chain->SetBranchAddress("vetoMulti", &spirit.hime_veto_multi);
	chain->SetBranchAddress("vetoBarId", &spirit.hime_veto_bar[0]);
	chain->SetBranchAddress("vetoTof", &spirit.hime_veto_tof[0]);
	chain->SetBranchAddress("vetoTot", &spirit.hime_veto_charge[0]);
	chain->SetBranchAddress("vetoTdiff", &spirit.hime_veto_tdiff[0]);

	chain->SetBranchAddress("tbdc_x", &spirit.tbdc_x);
	chain->SetBranchAddress("tbdc_y", &spirit.tbdc_y);
	chain->SetBranchAddress("tbdc_a", &spirit.tbdc_a);
	chain->SetBranchAddress("tbdc_b", &spirit.tbdc_b);

	return chain;
}

#endif
