{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect0{ 19.54383947216508, 19.54383947216508, 14.60830931266006, 14.60830931266006, 17.54169044519606, 19.54383947216508 };
   std::vector<Double_t> cutg_vect1{ 23.97807376494967, 20.97167003162618, 19.79390362063346, 22.67633194227352, 23.4821721182159, 23.97807376494967 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
