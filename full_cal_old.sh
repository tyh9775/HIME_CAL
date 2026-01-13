#!/bin/bash

set -e

all=1 #switch between using just june data (0) or all data (1)
original=0 #run Kin's original code
mu_st=0
limit=true
single=0 #run the single detector estimation
load=0 #for loading given data for calibration
cal=0 #for applying the calibration parameters
recut=0 #for automatically remaking tcutg files 
amp_e=0 #calibrating amplitude to energy
ecal=0 #creating E vs tot graphs
vis=0
muon=0
veto=0

Hel=false
E_min=-1 # neg = don't use; zero = use val in myConst; pos = use given value
read=1

#function to pull out parameters from txt file
get_params() {
    local det_line="$1"
    local file="$2"

    local line
    line=$(awk -v target="$det_line" '$0 == target {getline; print}' "$file")

    if [[ -z $line ]]
    then
        echo "Line not found: $det_line"
        return 1
    fi

    a1=$(echo "$line" | sed -E 's/a1: ([^+]+)\+.*a2: ([^+]+)\+.*ka: ([^+]+)\+.*/\1/')
    a2=$(echo "$line" | sed -E 's/a1: ([^+]+)\+.*a2: ([^+]+)\+.*ka: ([^+]+)\+.*/\2/')
    ka=$(echo "$line" | sed -E 's/a1: ([^+]+)\+.*a2: ([^+]+)\+.*ka: ([^+]+)\+.*/\3/')

    echo "{$a1,$a2,$ka}"
}

get_double() {
    local header="$1"
    local name="$2"
    grep -E "double[[:space:]]+${name}[[:space:]]*\{" "$header" \
      | sed -E 's/.*\{([^}]+)\}.*/\1/'
}

cut_file1="./TCuts/cuts_L0.root"
cut_file2="./TCuts/cuts_L1.root"
cut_file3="./TCuts/cuts_L2.root"

result_dir="sb_results"
mkdir -p $result_dir
amp_dir_0="amp_results_v3"
mkdir -p $amp_dir_0

amp_dir=$amp_dir_0

En=$E_min
vis_res="${amp_dir_0}/vis"
if [ $E_min -lt 0 ]
then
    vis_res="${vis_res}_no_E"
    amp_dir="${amp_dir}/no_E"
elif [ $E_min == 0 ]
then
    En=$(get_double myConst3.h En_min)
    vis_res="${vis_res}_min_E"
    amp_dir="${amp_dir}/min_E"
else
    vis_res="${vis_res}_min_E_cust"
    amp_dir="${amp_dir}/min_E_cust"
fi

if [ $Hel == true ]
then
    vis_res="${vis_res}_He"
    amp_dir="${amp_dir}_He"
else 
    vis_res="${vis_res}_no_He"
    amp_dir="${amp_dir}_no_He"
fi

vis_res_2="${vis_res}_v2"

mkdir -p $amp_dir

og_file="./hime/results/gain_matching/june-walk.root"
#og_file="./sb_results/gain_matching/june-walk.root" #if new loaded data is needed, use reload.sh
sb_res="${result_dir}/june"
amp_res="${amp_dir}/june"
amp_res_0="${amp_dir_0}/june"

if [ $all == 1 ]
then
    sb_res="${result_dir}/all"
    amp_res="${amp_dir}/all"
    amp_res_0="${amp_dir_0}/all"
    og_file="./hime/results/gain_matching/all-walk.root"
    #og_file="./sb_results/gain_matching/all-walk.root"
fi

tt_cut0="./TCuts/tt_L0.root"
tt_cut1="./TCuts/tt_L1.root"
tt_cut2="./TCuts/tt_L2.root"

E_cut0="./TCuts/E_L0.root"
E_cut1="./TCuts/E_L1.root"
E_cut2="./TCuts/E_L2.root"

ar_sc="${amp_res}-single-cal"
ar_par="${amp_res}-parameters"
ar_par_2="${amp_res}-parameters-v2"
ar_clb="${amp_res}-clb.root"
ar_clb_2="${amp_res}-clb_2.root"
ar_a_t="${amp_res}-a-t.root"
ar_a_t_2="${amp_res}-a-t-2.root"
ar_ae="${amp_res}-amp-E"
ar_ae_2="${amp_res}-amp-E-v2"
cut_dir="${amp_dir}/cuts"
mkdir -p $cut_dir
ar_cut0="${cut_dir}/amp_L0.root"
ar_cut1="${cut_dir}/amp_L1.root"
ar_cut2="${cut_dir}/amp_L2.root"
ar_cut0_v2="${cut_dir}/amp_L0_v2.root"
ar_cut1_v2="${cut_dir}/amp_L1_v2.root"
ar_cut2_v2="${cut_dir}/amp_L2_v2.root"
ar_sp="${amp_res}_simple"
ar_sp_2="${amp_res}_simple_v2"
ar_sv="${amp_res}_simple_vis"
ar_sv_2="${amp_res}_simple_vis_2"
ar_mu="${amp_res}_mu"
ar_mu_2="${amp_res}_mu_2"

if [ "$limit" == true ]
then
    ar_sc="${amp_res}-single-cal-lim"
    ar_par="${amp_res}-parameters-lim"
    ar_par_2="${amp_res}-parameters-v2-lim"
    ar_clb="${amp_res}-clb-lim.root"
    ar_clb_2="${amp_res}-clb_2-lim.root"
    ar_a_t="${amp_res}-a-t-lim.root"
    ar_a_t_2="${amp_res}-a-t-2-lim.root"
    ar_ae="${amp_res}-amp-E-lim"
    ar_ae_2="${amp_res}-amp-E-v2-lim"
    ar_cut0="${cut_dir}/amp_L0_lim.root"
    ar_cut1="${cut_dir}/amp_L1_lim.root"
    ar_cut2="${cut_dir}/amp_L2_lim.root"
    ar_cut0_v2="${cut_dir}/amp_L0_v2_lim.root"
    ar_cut1_v2="${cut_dir}/amp_L1_v2_lim.root"
    ar_cut2_v2="${cut_dir}/amp_L2_v2_lim.root"
    ar_sp="${amp_res}_simple_lim"
    ar_sp_2="${amp_res}_simple_v2_lim"
    ar_sv="${amp_res}_simple_vis_lim"
    ar_sv_2="${amp_res}_simple_vis_2_lim"
    ar_mu="${amp_res}_mu_lim"
    ar_mu_2="${amp_res}_mu_2_lim"

    vis_res="${vis_res}_lim"
    vis_res_2="${vis_res_2}_lim"
fi


runs23=( 1148 1149 1150 1151 1152 1153 1154 1155 1156 1158 1160 1161 1162 1163 1164 1165 1166 1167 1168 1169 1170 1171 1172 1174 1175 1176 1178 1179 1180 1181 1182 1183 1184 1185 1186 1189 1190 1191 1192 1195 1196 1197 1198 1199 1200 1201 1202 1203 1204 1205 1206 1207 1208 1210 )

runs34=( 1211 1212 1213 1214 1215 1217 1218 1219 1221 1222 1223 1225 1226 1227 1232 1234 1235 1236 1239 1240 1241 1242 1243 1244 1245 )

run34_nov=( 1644  1645  1646  1647  1648  1649  1650  1651  1652  1653  1654  1655  1656  1657  1658  1659  1660  1661  1662  1663  1664  1665  1666  1667  1668  1669  1670  1671  1672  1673  1675  1676  1677  1679  1680  1681  1682  1683  1684  1685  1686  1687  1688  1689  1690  1691  1692  1693  1694  1695  1696  1697  1698  1699  1700  1701  1702  1703  1704  1705  1706  1707  1708  1709  1710  1711  1721  1729  1730  1731  1732  1733  1734  1735  1736  1737  1738  1739  1740  1741  1742  1743  1744  1745  1746  1747  1748  1750  1751  1752  1753  1754  1755  1756  1757  1758  1759  1760  1761  1762  1763  1764  1765  1766  1767  1768  1769  1770  1771  1772  1773  1774  1775  1776  1777  1778  1779  1780  1781  1782  1783  1784  1785  1786  1788  1790  1791  1792  1793  1794  1800  1801  1802  1803  1804  1805  1806  1807  1808  1809  1810  1811  1812  1813  1814  1815  1816  1817  1818  1819  1822  1823  1824  1825  1826  1827  1828  1829  1830  1831  1832  1833  1834  1835  1836  1837  1838  1839  1840  1841  1842  1843  1844  1845  1846  1847  1848  1849  1850  1851  1852  1853  1854  1855  1856  1857  1858  1859  1860  1861  1862  1863  1864  1865  1866  1867  1868  1869 )

# shadow bar 4, at bottom, backward angle
runs4a=(1247 1248 1251 1252 1253 1254 1255 1256 1261 1262 1263 1264 1265 1266 1268 1269 1273 1274 1275 1276 1277 1279 1280 1281 1282 1283 1284 1285 1286 1287 1288 1289 1290 1291 1292 1293 1294 1295 1296 1297 1298 1299)

# shadow bar 4, at bottom, forward angle, but target switched to Xe132
runs4b=(1300 1304 1305 1306 1308 1309 1310 1311 1312 1313 1314 1315 1316 1317 1318 1319 1320 1321 1322 1324 1325 1326 1327 1331)

#every_run=( ${runs23[@]} ${runs34[@]} ${runs4a[@]} ${runs4b[@]} ${run34_nov[@]} )
runs_june=( ${runs23[@]} ${runs34[@]} ${runs4a[@]} ${runs4b[@]} )
runs_nov=( ${run34_nov[@]} )
runs_all=( ${runs_june[@]} ${runs_nov[@]} )

runs_j_2=( ${runs34[@]} ${runs4a[@]} ${runs4b[@]} ) #june run without sb 2
runs_j_3=( ${runs4a[@]} ${runs4b[@]} ) #june run without sb 3
runs_j_4=( ${runs23[@]} ) #june run without sb 4

runs_a_2=( ${runs34[@]} ${runs4a[@]} ${runs4b[@]} ${run34_nov[@]} ) #all run w/o sb 2
runs_a_3=( ${runs4a[@]} ${runs4b[@]} ) #all run w/o sb 3
runs_a_4=( ${runs23[@]} ) #all run w/o sb4

runs_j_2_inc=( ${runs23[@]} ) #june run including sb 2
runs_j_3_inc=( ${runs23[@]} ${runs34[@]} ) #june run inc sb 3
runs_j_4_inc=( ${runs34[@]} ${runs4a[@]} ${runs4b[@]} ) #june run inc sb 4

runs_a_2_inc=( ${runs23[@]} )
runs_a_3_inc=( ${runs23[@]} ${runs34[@]} ${run34_nov[@]} )
runs_a_4_inc=( ${runs34[@]} ${runs4a[@]} ${runs4b[@]} ${run34_nov[@]} ) 


runs_june=$( IFS=, ; echo "${runs_june[*]}" )
runs_nov=$( IFS=, ; echo "${runs_nov[*]}" )
runs_all=$( IFS=, ; echo "${runs_all[*]}" )

runs_j_2=$( IFS=, ; echo "${runs_j_2[*]}" )
runs_j_3=$( IFS=, ; echo "${runs_j_3[*]}" )
runs_j_4=$( IFS=, ; echo "${runs_j_4[*]}" )

runs_a_2=$( IFS=, ; echo "${runs_a_2[*]}" )
runs_a_3=$( IFS=, ; echo "${runs_a_3[*]}" )
runs_a_4=$( IFS=, ; echo "${runs_a_4[*]}" )

runs_j_2_inc=$( IFS=, ; echo "${runs_j_2_inc[*]}" )
runs_j_3_inc=$( IFS=, ; echo "${runs_j_3_inc[*]}" )
runs_j_4_inc=$( IFS=, ; echo "${runs_j_4_inc[*]}" )

runs_a_2_inc=$( IFS=, ; echo "${runs_a_2_inc[*]}" )
runs_a_3_inc=$( IFS=, ; echo "${runs_a_3_inc[*]}" )
runs_a_4_inc=$( IFS=, ; echo "${runs_a_4_inc[*]}" )



runs=$runs_june
runs_2=$runs_j_2 #runs without sb 2
runs_3=$runs_j_3 #runs without sb 3
runs_4=$runs_j_4  #runs without sb 4

runs_2_inc=$runs_j_2_inc #runs inc sb 2
runs_3_inc=$runs_j_3_inc #runs inc sb 3
runs_4_inc=$runs_j_4_inc  #runs inc sb 4

if [ $all == 1 ]
then
    runs=$runs_all
    runs_2=$runs_a_2 
    runs_3=$runs_a_3 
    runs_4=$runs_a_4 

    runs_2_inc=$runs_a_2_inc 
    runs_3_inc=$runs_a_3_inc 
    runs_4_inc=$runs_a_4_inc 
fi

if [ $original == 1 ] 
then
    root -l -b -q "shadowAna.C({$runs}, \"${sb_res}-shadow-og.root\", true, false)"
    root -l -b -q "shadowAna.C({$runs}, \"${sb_res}-shadow-vw-og.root\", true, true)"
fi

if [ $mu_st == 1 ]
then
    root -b -q -l "muon_load.cpp(\"ma.root\",\"mu_ld\")"
fi


if [ $single == 1 ] 
then
    #root -b -q -l "single_v2.cpp({$runs},\"${amp_res_0}-single.root\",{14,38,62},true,false)"
    root -b -q -l "single_cal_v3.cpp(\"${amp_res_0}-single.root\",\"mu_ld.txt\",\"$ar_sc\",{\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},$limit,$Hel,$En)"
fi

a_pars=()   
for det in "det 14" "det 38" "det 62"
do
    params=$(get_params "$det" "${ar_sc}.txt")
    if [[ $? -eq 0 ]]
    then
        a_pars+=("$params")
    fi
done
a_mat="{${a_pars[0]},${a_pars[1]},${a_pars[2]}}"

a_pars_2=()   
for det in "det 14 version 2" "det 38 version 2" "det 62 version 2"
do
    params=$(get_params "$det" "${ar_sc}.txt")
    if [[ $? -eq 0 ]]
    then
        a_pars_2+=("$params")
    fi
done
a_mat_2="{${a_pars_2[0]},${a_pars_2[1]},${a_pars_2[2]}}"


if [ $load == 1 ]
then   
    #root -b -q -l "mult_loader_v3.cpp(\"${sb_res}-shadow-og.root\",\"ml.root\",{\"${cut_file1}\",\"${cut_file2}\",\"${cut_file3}\"})"
    root -b -q -l "ToT_cal_v5.cpp(\"ml.root\",\"${ar_clb}\",\"$ar_par\",true,$a_mat)"
    root -b -q -l "ToT_cal_v5.cpp(\"ml.root\",\"${ar_clb_2}\",\"$ar_par_2\",true,$a_mat_2)"
fi

if [ $cal == 1 ]
then
    root -b -q -l "amp_tot.C({$runs},\"$ar_a_t\",\"${ar_par}.txt\",true,false,$a_mat)"
    root -b -q -l "amp_tot.C({$runs},\"$ar_a_t_2\",\"${ar_par_2}.txt\",true,false,$a_mat_2)"
fi

if [ $recut == 1 ]
then
    root -l -b -q "amp_cut_remaker.cpp({\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$ar_cut0\",\"$ar_cut1\",\"$ar_cut2\"},$a_mat)"
    root -l -b -q "amp_cut_remaker.cpp({\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$ar_cut0_v2\",\"$ar_cut1_v2\",\"$ar_cut2_v2\"},$a_mat_2)"
fi

if [ $amp_e == 1 ]
then 
    root -l -b -q "amp_E_cal_v2.cpp(\"$ar_a_t\",\"mu_ld.txt\",\"$ar_ae\",{\"$ar_cut0\",\"$ar_cut1\",\"$ar_cut2\"},$a_mat,true,true,$En)" #set offset to 0
    root -l -b -q "amp_E_cal_v2.cpp(\"$ar_a_t\",\"mu_ld.txt\",\"${ar_ae}-os\",{\"$ar_cut0\",\"$ar_cut1\",\"$ar_cut2\"},$a_mat,false,true,$En)" 

    root -l -b -q "amp_E_cal_v2.cpp(\"$ar_a_t_2\",\"mu_ld.txt\",\"$ar_ae_2\",{\"$ar_cut0_v2\",\"$ar_cut1_v2\",\"$ar_cut2_v2\"},$a_mat_2,true,true,$En)" 
    root -l -b -q "amp_E_cal_v2.cpp(\"$ar_a_t_2\",\"mu_ld.txt\",\"${ar_ae_2}-os\",{\"$ar_cut0_v2\",\"$ar_cut1_v2\",\"$ar_cut2_v2\"},$a_mat_2,false,true,$En)" 
fi

if [ $ecal == 1 ]
then 
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp}.root\",\"${ar_ae}.txt\",\"${ar_par}.txt\", true, false,$a_mat)"
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp}_os.root\",\"${ar_ae}-os.txt\",\"${ar_par}.txt\", true, false,$a_mat)"

    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp_2}.root\",\"${ar_ae_2}.txt\",\"${ar_par_2}.txt\", true, false,$a_mat_2)"
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp_2}_os.root\",\"${ar_ae_2}-os.txt\",\"${ar_par_2}.txt\", true, false,$a_mat_2)"
fi

if [ $vis == 1 ]
then
    root -l -b -q "visualize_simple.cpp(\"${ar_sp}.root\",\"${ar_sv}.root\")"
    root -l -b -q "visualize_simple.cpp(\"${ar_sp}_os.root\",\"${ar_sv}_os.root\")"
    
    root -l -b -q "visualize_simple.cpp(\"${ar_sp_2}.root\",\"${ar_sv_2}.root\")"
    root -l -b -q "visualize_simple.cpp(\"${ar_sp_2}_os.root\",\"${ar_sv_2}_os.root\")"
fi

if [ $read == 1 ]
then
    root -b -q -l "reader.cpp(\"${ar_sp}.root\",\"${vis_res}.root\",{\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$E_cut0\",\"$E_cut1\",\"$E_cut2\"})"
    root -b -q -l "reader.cpp(\"${ar_sp_2}.root\",\"${vis_res_2}.root\",{\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$E_cut0\",\"$E_cut1\",\"$E_cut2\"})"

    root -b -q -l "reader.cpp(\"${ar_sp}_os.root\",\"${vis_res}_os.root\",{\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$E_cut0\",\"$E_cut1\",\"$E_cut2\"})"
    root -b -q -l "reader.cpp(\"${ar_sp_2}_os.root\",\"${vis_res_2}_os.root\",{\"$tt_cut0\",\"$tt_cut1\",\"$tt_cut2\"},{\"$E_cut0\",\"$E_cut1\",\"$E_cut2\"})"
fi


if [ $muon == 1 ]
then
    root -b -q -l "muon_pts.cpp(\"ma.root\",\"${ar_mu}.root\",$a_mat,\"${ar_ae}.txt\")"
    root -b -q -l "muon_pts.cpp(\"ma.root\",\"${ar_mu}-os.root\",$a_mat,\"${ar_ae}-os.txt\")"

    root -b -q -l "muon_pts.cpp(\"ma.root\",\"${ar_mu_2}.root\",$a_mat_2,\"${ar_ae_2}.txt\")"
    root -b -q -l "muon_pts.cpp(\"ma.root\",\"${ar_mu_2}-os.root\",$a_mat_2,\"${ar_ae_2}-os.txt\")"
fi


if [ $veto == 1 ]
then
    #root -l -b -q "tot_simple_v2.cpp({$runs},\"${amp_res}_simple_vw.root\",\"${amp_res}-amp-E-z.txt\",\"parameters_amp.txt\", true, true,{$a1, $a2, $ka})"
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp}_vw.root\",\"${ar_ae}.txt\",\"${ar_par}.txt\", true, true,$a_mat)"
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp}_os_vw.root\",\"${ar_ae}-os.txt\",\"${ar_par}.txt\", true, true,$a_mat)"

    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp_2}_vw.root\",\"${ar_ae_2}.txt\",\"${ar_par_2}.txt\", true, true,$a_mat_2)"
    root -l -b -q "tot_simple_v3.cpp({$runs},\"${ar_sp_2}_os_vw.root\",\"${ar_ae_2}-os.txt\",\"${ar_par_2}.txt\", true, true,$a_mat_2)"
fi