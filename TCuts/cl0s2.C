{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect4{ 33.21361039734546, 34.49584537347231, 40.45823801246216, 40.45823801246216, 39.24011478514166, 33.59828089018352, 32.0595989188313, 32.0595989188313, 33.21361039734546 };
   std::vector<Double_t> cutg_vect5{ 27.73750005941838, 27.73750005941838, 23.99802631948535, 21.55592101993725, 21.63223681054813, 26.0967105612845, 25.67697371292467, 27.50855268758575, 27.73750005941838 };
   TCutG *cutg = new TCutG("CUTG", 9, cutg_vect4.data(), cutg_vect5.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
