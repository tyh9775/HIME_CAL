{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect8{ 19.10257879510992, 21.51375361613952, 21.39197711002691, 19.10257879510992, 18.81031518043967, 19.10257879510992 };
   std::vector<Double_t> cutg_vect9{ 30.46552635252868, 30.46552635252868, 25.84710523107726, 25.98921049635269, 30.02500003017485, 30.46552635252868 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect8.data(), cutg_vect9.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
