#include "spirit.h"


void get_num_entries(
    const std::vector<std::vector<int>> runIds = {{123,124,125},{123,124,125}}, 
    const std::string &ofile = "num_entries.txt",
    //const std::string &spiritDir = "hime/spirit"
    const std::string &spiritDir="hime_riken_bdc"

){
    std::ofstream output;
    output.open(ofile.c_str(),std::ios_base::trunc);  
    for (size_t i=0;i<runIds.size();i++){
        auto spiritChain = getSpiritChain(runIds[i], spiritDir);
        auto spiritEntries = spiritChain->GetEntries();
        std::cout<<spiritEntries<<"\n";
        output<<spiritEntries<<"\n";
    }
    output.close();
}