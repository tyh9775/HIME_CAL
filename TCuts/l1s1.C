{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect0{ 19.19484231292262, 23.3982807709742, 23.3982807709742, 19.01969904383714, 19.25322340261778, 19.19484231292262 };
   std::vector<Double_t> cutg_vect1{ 25.42289472824257, 25.42289472824257, 21.74552625239287, 21.71868414672974, 24.64447366401161, 25.42289472824257 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
