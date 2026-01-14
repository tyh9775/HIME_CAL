{
//========= Macro generated from object: CUTG/Graph
//========= by ROOT version6.36.04
   
   std::vector<Double_t> cutg_vect16{ 14.87843837423256, 16.61912607925392, 16.73517192625534, 14.26919767747508, 14.15315183047365, 14.87843837423256 };
   std::vector<Double_t> cutg_vect17{ 30.99500001929701, 31.01447370379772, 26.74973679814292, 26.6913157446408, 30.1186842167652, 30.99500001929701 };
   auto *cutg = new TCutG("CUTG", 6, cutg_vect16.data(), cutg_vect17.data());
   cutg->SetVarX("");
   cutg->SetVarY("");
   cutg->SetTitle("Graph");
   cutg->SetFillStyle(1000);
   cutg->Draw();
}
