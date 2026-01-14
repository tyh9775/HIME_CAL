{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect12{ 23.81547277065951, 26.77249287202913, 26.73810891736204, 23.67793695199115, 23.47163322398862, 23.81547277065951 };
   std::vector<Double_t> cutg_vect13{ 27.42710531273563, 27.37763162778789, 22.50447366043533, 22.62815787280469, 26.33868424388531, 27.42710531273563 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect12.data(), cutg_vect13.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
