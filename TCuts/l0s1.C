{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect0{ 23.94054438191917, 27.25358168057095, 27.11031520279141, 23.51074494858056, 23.61819480691521, 23.94054438191917 };
   std::vector<Double_t> cutg_vect1{ 27.07960528893102, 27.10855265778343, 23.31644733811876, 23.35986839139736, 27.09407897335723, 27.07960528893102 };
   TCutG *cutg = new TCutG("CUTG", 6, cutg_vect0.data(), cutg_vect1.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
