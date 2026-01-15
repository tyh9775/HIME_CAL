{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect0{ 16.3522204612774, 16.3522204612774, 23.91740682042581, 23.91740682042581, 22.76339534191165, 16.3522204612774 };
   std::vector<Double_t> cutg_vect1{ 23.76907894765272, 20.71644732321759, 22.7769736697113, 25.82960529414643, 25.48618423639748, 23.76907894765272 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
