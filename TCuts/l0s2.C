{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect2{ 31.33008593798198, 33.75415474201173, 33.6381088950103, 31.22693407398071, 31.23982805698087, 31.33008593798198 };
   std::vector<Double_t> cutg_vect3{ 27.78815791691213, 27.94078949813388, 24.29289470693391, 24.35394733942261, 27.66605265193472, 27.78815791691213 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect2.data(), cutg_vect3.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
