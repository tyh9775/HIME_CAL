{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect6{ 16.93495698630553, 18.55601718238098, 18.62263609454847, 18.60042979049264, 16.13553004029572, 16.06891112812824, 16.93495698630553 };
   std::vector<Double_t> cutg_vect7{ 30.04000003039837, 30.02736845126278, 29.04210527868647, 26.12421049836435, 26.12421049836435, 29.11789475350004, 30.04000003039837 };
   auto *cutg = new TCutG("CUTG", 7, cutg_vect6.data(), cutg_vect7.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
