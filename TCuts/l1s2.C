{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect2{ 26.49247852481773, 29.06124647140482, 29.35315191988062, 26.72600288359838, 26.37571634542741, 26.49247852481773 };
   std::vector<Double_t> cutg_vect3{ 26.46973684910489, 26.44289474344176, 22.63131573927638, 22.79236837325519, 25.77184210186334, 26.46973684910489 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect2.data(), cutg_vect3.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
