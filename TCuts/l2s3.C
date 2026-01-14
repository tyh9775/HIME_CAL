{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect14{ 28.12492842226795, 31.3842407917524, 31.5704872128658, 28.0783668169896, 27.65931236948446, 28.12492842226795 };
   std::vector<Double_t> cutg_vect15{ 28.39473690610183, 28.26000006198883, 20.78210521371741, 20.8494736357739, 27.11473688702835, 28.39473690610183 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect14.data(), cutg_vect15.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
