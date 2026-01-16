#!/bin/bash

#this is the shell for running the full calibration process
#Ultimately want to get Energy from Time over Threshold (ToT)
#Energy=c*Amplitude+d
#Amplitude=A*exp(k*ToT)+B


set -e

#switches for running different parts of the code
original=false #run Kin's original code
mu_st=false #load and combine cosmic muon data
single=false #run the single detector estimation
load=true #for loading given data for calibration
cal=false #for applying the calibration parameters
recut=false #for automatically remaking tcutg files 
amp_e=false #calibrating amplitude to energy
ecal=false #creating E vs tot graphs

read=1

#switches for different output options
Hel=false


#functions####################################################################################
#function to pull out parameters from txt file

get_double() {
    local header="$1"
    local name="$2"
    grep -E "double[[:space:]]+${name}[[:space:]]*\{" "$header" \
      | sed -E 's/.*\{([^}]+)\}.*/\1/'
}
##############################################################################################

#paths and filenames
tot_code_dir="codes/tot_cal"
sb_code_dir="codes/shadow_bar"

mu_dir="cosmic_muons"
mu_pattern="${mu_dir}/cosmics_layers_1_2_3_thr0x200_volt2024-06-25_00"

result_dir="tot_cal_results"
mkdir -p $result_dir

if [ $Hel == true ]
then
    res_dir="${result_dir}/He"
    mkdir -p "$res_dir"
else 
    res_dir="${result_dir}/no_He"
    mkdir -p "$res_dir"
fi

s_est="${res_dir}/single_estimate"



runs23=( 1148 1149 1150 1151 1152 1153 1154 1155 1156 1158 1160 1161 1162 1163 1164 1165 1166 1167 1168 1169 1170 1171 1172 1174 1175 1176 1178 1179 1180 1181 1182 1183 1184 1185 1186 1189 1190 1191 1192 1195 1196 1197 1198 1199 1200 1201 1202 1203 1204 1205 1206 1207 1208 1210 )

runs34=( 1211 1212 1213 1214 1215 1217 1218 1219 1221 1222 1223 1225 1226 1227 1232 1234 1235 1236 1239 1240 1241 1242 1243 1244 1245 )

run34_nov=(  1644  1645  1646  1647  1648  1649  1650  1651  1652  1653  1654  1655  1656  1657  1658  1659  1660  1661  1662  1663  1664  1665  1666  1667  1668  1669  1670  1671  1672  1673  1675  1676  1677  1679  1680  1681  1682  1683  1684  1685  1686  1687  1688  1689  1690  1691  1692  1693  1694  1695  1696  1697  1698  1699  1700  1701  1702  1703  1704  1705  1706  1707  1708  1709  1710  1711  1721  1729  1730  1731  1732  1733  1734  1735  1736  1737  1738  1739  1740  1741  1742  1743  1744  1745  1746  1747  1748  1750  1751  1752  1753  1754  1755  1756  1757  1758  1759  1760  1761  1762  1763  1764  1765  1766  1767  1768  1769  1770  1771  1772  1773  1774  1775  1776  1777  1778  1779  1780  1781  1782  1783  1784  1785  1786  1788  1790  1791  1792  1793  1794  1800  1801  1802  1803  1804  1805  1806  1807  1808  1809  1810  1811  1812  1813  1814  1815  1816  1817  1818  1819  1822  1823  1824  1825  1826  1827  1828  1829  1830  1831  1832  1833  1834  1835  1836  1837  1838  1839  1840  1841  1842  1843  1844  1845  1846  1847  1848  1849  1850  1851  1852  1853  1854  1855  1856  1857  1858  1859  1860  1861  1862  1863  1864  1865  1866  1867  1868  1869  )

# shadow bar 4, at bottom, backward angle
runs4a=(1247 1248 1251 1252 1253 1254 1255 1256 1261 1262 1263 1264 1265 1266 1268 1269 1273 1274 1275 1276 1277 1279 1280 1281 1282 1283 1284 1285 1286 1287 1288 1289 1290 1291 1292 1293 1294 1295 1296 1297 1298 1299)

# shadow bar 4, at bottom, forward angle, but target switched to Xe132
runs4b=(1300 1304 1305 1306 1308 1309 1310 1311 1312 1313 1314 1315 1316 1317 1318 1319 1320 1321 1322 1324 1325 1326 1327 1331)

runs_test=( 1870 1871 1872 1873 1874 1875 1876 1877 1878 1879 1880 1881 1882 1883 1884 1885 1886 1887 1888 1889 1890 1891 1892 1893 1894 1895 1896 1897 1898 1899 1900 1901 1902 1903 1904 1905 1906 1907 1908 1909 1910 1911 1912 1913 1914 1915 1916 )
#runs_t=(${runs_test[@]})
#runs_t=$( IFS=, ; echo "${runs_t[*]}" )

runs_23=( ${runs23[@]} )
runs_34=( ${runs34[@]} )
#runs_34n=( ${run34_nov[@]} )
runs_34n=( ${run34_nov[@]} ${runs_test[@]} )
runs_4=( ${runs4a[@]} ${runs4b[@]} )

runs_23=$( IFS=, ; echo "${runs_23[*]}" )
runs_34=$( IFS=, ; echo "${runs_34[*]}" )
runs_34n=$( IFS=, ; echo "${runs_34n[*]}" )
runs_4=$( IFS=, ; echo "${runs_4[*]}" )

runs_34t=( ${runs34[@]} ${run34_nov[@]} ${runs_test[@]} )
runs_34t=$( IFS=, ; echo "${runs_34t[*]}" )

#runs_all=( ${runs23[@]} ${runs34[@]} ${run34_nov[@]} ${runs4a[@]} ${runs4b[@]} )
runs_all=( ${runs23[@]} ${runs34[@]} ${run34_nov[@]} ${runs_test[@]} ${runs4a[@]} ${runs4b[@]} )

runs_all=$( IFS=, ; echo "${runs_all[*]}" )



if [ $original == true ] 
then
    root -l -b -q "${sb_code_dir}/shadowAna.C({$runs_all}, \"${result_dir}/og.root\", true, false)"
    root -l -b -q "${sb_code_dir}/shadowAna.C({$runs_all}, \"${result_dir}/vw-og.root\", true, true)"
fi

if [ $mu_st == true ]
then
    root -b -q -l "${tot_code_dir}/muon_add.cpp(\"${mu_pattern}\",\"${mu_dir}/ma.root\")"
    root -b -q -l "${tot_code_dir}/muon_load.cpp(\"${mu_dir}/ma.root\",\"${mu_dir}/mu_ld\")"
fi


if [ $single == true ] 
then
    #root -b -q -l "${tot_code_dir}/single.cpp({$runs_all},\"${result_dir}/single.root\",{14,38,62},true,false)"
    root -b -q -l "${tot_code_dir}/single_estimate.cpp(\"${result_dir}/single.root\",\"${mu_dir}/mu_ld.txt\",\"$s_est\",\"./TCuts/l\",$Hel)"
fi


if [ $load == true ]
then
    root -b -q -l "${tot_code_dir}/mult_loader.cpp(\"${result_dir}/og.root\",\"${res_dir}/ml.root\",\"./TCuts/cal_l\",false,true)"
    #root -b -q -l "ToT_cal.cpp(\"ml.root\",\"${ar_clb}\",\"$ar_par\",true,$a_mat)"

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