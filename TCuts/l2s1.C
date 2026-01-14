{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect10{ 18.00458443192152, 20.23954148528228, 20.03323775727975, 17.38567324791392, 17.28252138391266, 18.00458443192152 };
   std::vector<Double_t> cutg_vect11{ 26.38815792883306, 26.41289477130693, 21.02026311200308, 20.99552626952921, 25.10184212019177, 26.38815792883306 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect10.data(), cutg_vect11.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
