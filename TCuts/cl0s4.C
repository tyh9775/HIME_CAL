{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect10{ 21.5925501640171, 21.5925501640171, 14.20429790492642, 14.20429790492642, 17.95644695797247, 21.5925501640171 };
   std::vector<Double_t> cutg_vect11{ 31.06513161344552, 28.04539472634267, 25.42828942418686, 27.91118419802699, 29.09671053148218, 31.06513161344552 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect10.data(), cutg_vect11.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
