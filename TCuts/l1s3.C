{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect4{ 31.39649005921126, 33.32306601915157, 34.37392563366446, 34.54906890274995, 31.57163332829674, 31.04620352104029, 31.39649005921126 };
   std::vector<Double_t> cutg_vect5{ 27.24815791333585, 27.16763159634644, 26.57710527175744, 20.8060525541831, 20.96710518816192, 26.71131580007312, 27.24815791333585 };
   auto *cutg = new TCutG("CUTG", 7, cutg_vect4.data(), cutg_vect5.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
