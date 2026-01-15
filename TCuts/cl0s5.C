{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect12{ 24.6484241350546, 24.6484241350546, 28.16848149409781, 28.16848149409781, 25.38338116606363, 24.6484241350546 };
   std::vector<Double_t> cutg_vect13{ 31.53486846255041, 29.61118422335896, 27.08355260674695, 29.09671053148218, 31.5796053053223, 31.53486846255041 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect12.data(), cutg_vect13.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
