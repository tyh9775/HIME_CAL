{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect18{ 16.73517192625534, 18.79498571053063, 18.88202009578169, 16.73517192625534, 16.44505730875179, 16.73517192625534 };
   std::vector<Double_t> cutg_vect19{ 31.55973686981751, 31.5402631853168, 26.96394732765069, 26.92499995864928, 30.9755263347963, 31.55973686981751 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect18.data(), cutg_vect19.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
