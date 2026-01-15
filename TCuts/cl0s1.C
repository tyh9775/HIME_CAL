{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect2{ 23.98151856923215, 25.32786529416535, 27.25121775835562, 33.59828089018352, 34.23939837824695, 34.23939837824695, 31.54670492838056, 26.35365327506683, 23.98151856923215, 23.98151856923215 };
   std::vector<Double_t> cutg_vect3{ 25.90592108475731, 26.59276320025521, 26.55460530494977, 22.43355261196235, 21.67039470585357, 20.10592099833057, 20.22039468424688, 24.91381580681589, 24.26513158662343, 25.90592108475731 };
   TCutG *cutg = new TCutG("CUTG", 10, cutg_vect2.data(), cutg_vect3.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
