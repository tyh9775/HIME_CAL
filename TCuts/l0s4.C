{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect0{ 19.16475637494066, 21.89398277664082, 21.97994266330855, 19.42263603494382, 19.27220623327531, 19.16475637494066 };
   std::vector<Double_t> cutg_vect1{ 31.24986844489253, 31.20013160204613, 27.90092102656827, 27.86776313133733, 30.43750001173467, 31.24986844489253 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
