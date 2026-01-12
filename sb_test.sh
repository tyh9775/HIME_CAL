#!/bin/bash

set -e

par_type=0 #choose which folder to pull parameters from; 0 = no E no He, 1 = no E He
use_veto=true #switch for veto detectors; use to filter out charged particles

lim=true #were the initial parameters limited
offset=false #did the linear params have an offset

evt_norm=true

initialize=true
num_check=true
scl_check=true
combine=true
compare=true
percent_grid=true
mult_check=true



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

get_mpar() {

    declare -n out=$1
    local file=$2

    out=()

    local a b c
    local rows=()
    while IFS=',' read -r a b c; do
        a=$(echo "$a" | xargs)
        b=$(echo "$b" | xargs)
        c=$(echo "$c" | xargs)
        rows+=("{$a,$b,$c}")
    done < "$file"

    out=$( IFS=, ; echo "${rows[*]}" )
}


amp_dir_0="amp_results_v3" #input directory
result_dir="sb_test_results" #output directory

if [ $evt_norm == true ]
then
    result_dir="${result_dir}/evt_norm"
elif [ $evt_norm == false ]
then
    result_dir="${result_dir}/hit_norm"
fi

if [ $par_type == 0 ]
then
    amp_dir="${amp_dir_0}/no_E_no_He"
    result_dir="${result_dir}/no_E_no_He"
elif [ $par_type == 1 ]
then
    amp_dir="${amp_dir_0}/no_E_He"
    result_dir="${result_dir}/no_E_He"
fi

mkdir -p $result_dir


amp_res="${amp_dir}/all"
amp_res_0="${amp_dir_0}/all"

par_name="${amp_res}-parameters"
par_name_2="${amp_res}-parameters-v2"

lin_par_name="${amp_res}-amp-E"
lin_par_name_2="${amp_res}-amp-E-v2"

ar_sc="${amp_res}-single-cal"

if [ $lim == true ]
then
    par_name="${par_name}-lim"
    par_name_2="${par_name_2}-lim"
    lin_par_name="${lin_par_name}-lim"
    lin_par_name_2="${lin_par_name_2}-lim"
    ar_sc="${ar_sc}-lim"
fi

if [ $offset == true ]
then
    lin_par_name="${lin_par_name}-os"
    lin_par_name_2="${lin_par_name_2}-os"
fi


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

#root -b -q -l "shadowAna_v4.C({$runs_all},{\"${result_dir}/runs-all-tot-test.root\",\"${result_dir}/runs-all-hit-test.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},1,false)"
#root -b -q -l "shadowAna_v4.C({$runs_34t},{\"${result_dir}/runs-34-tot-test.root\",\"${result_dir}/runs-34-hit-test.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},1,false)"

if [ $initialize == true ]
then
    root -b -q -l "shadowAna_rbs.C({$runs_23},{\"${result_dir}/runs23-tot.root\",\"${result_dir}/runs23-hit.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},$evt_norm)"

    #root -b -q -l "shadowAna_rbs.C({$runs_34},{\"${result_dir}/runs34-tot.root\",\"${result_dir}/runs34-hit.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},$evt_norm)"
    #root -b -q -l "shadowAna_rbs.C({$runs_34n},{\"${result_dir}/runs34n-tot.root\",\"${result_dir}/runs34n-hit.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},$evt_norm)"
    #root -b -q -l "shadowAna_rbs.C({$runs_4},{\"${result_dir}/runs4-tot.root\",\"${result_dir}/runs4-hit.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},$evt_norm)"

    #root -b -q -l "shadowAna_rbs.C({$runs_all},{\"${result_dir}/runs-all-tot.root\",\"${result_dir}/runs-all-hit.root\"},\"${par_name}.txt\",true,${use_veto},${a_mat},$evt_norm)"
fi


if [ $num_check == true ]
then
    root -b -q -l "get_num_entries.cpp({{$runs_23},{$runs_34},{$runs_34n},{$runs_4}},\"${result_dir}/num_entries.txt\")"
    #root -b -q -l "get_num_entries.cpp({{$runs_23},{$runs_34},{$runs_34n},{$runs_t},{$runs_4}},\"num_entries_test.txt\")"
    if [ $evt_norm == false ]
    then
        root -b -q -l "get_num_hits.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\",\"${result_dir}/runs4-hit.root\"},\"${result_dir}/num_hits.txt\")"
    fi
fi

num_entries=($(cat "${result_dir}/num_entries.txt"))
if [ $evt_norm == false ]
then
    num_entries=""
    num_entries=($(cat "${result_dir}/num_hits.txt"))
fi
num_list="{$(printf "%s," "${num_entries[@]}" | sed 's/,$//')}"
num_list2="{$(printf "%s," "${num_entries[@]:1}" | sed 's/,$//')}"
num_list3="{$(printf "%s," "${num_entries[@]::${#num_entries[@]}-1}" | sed 's/,$//')}"

num_list34_comp="{$(printf "%s," "${num_entries[@]:1:2}" | sed 's/,$//')}"
#echo $num_list
#echo $num_list2
#echo $num_list3
if [ $scl_check == true ]
then
    root -b -q -l "scale_check_v2.cpp({\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\"},\"${result_dir}/runs34_comp\")"
fi

scl_list=($(cat ${result_dir}/runs34_comp.txt))
scl_ftr=${scl_list[0]}
scl_ftr_rnd=${scl_list[1]}

#echo $scl_ftr
#echo $scl_ftr_rnd

if [ $combine == true ]
then 
    #with no scaling
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\",\"${result_dir}/runs4-hit.root\"},\"${result_dir}/noSB2-hit-no-scl.root\",${use_veto},$num_list2,1.0)"
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\"},\"${result_dir}/SB3-hit-no-scl.root\",${use_veto},$num_list3,1.0)"

    #with default scaling factor (1.23)
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\",\"${result_dir}/runs4-hit.root\"},\"${result_dir}/noSB2-hit.root\",${use_veto},$num_list2)" #no sb2 is the same as inc sb 4
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\"},\"${result_dir}/SB3-hit.root\",${use_veto},$num_list3)" 

    #with new scaling factor
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\",\"${result_dir}/runs4-hit.root\"},\"${result_dir}/noSB2-hit-new.root\",${use_veto},$num_list2,$scl_ftr)" 
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\"},\"${result_dir}/SB3-hit-new.root\",${use_veto},$num_list3,$scl_ftr)" 

    #combine all
    root -b -q -l "sb_combiner_v2.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/runs34-hit.root\",\"${result_dir}/runs34n-hit.root\",\"${result_dir}/runs4-hit.root\"},\"${result_dir}/all-cmb-hit.root\",${use_veto},$num_list,$scl_ftr)"

fi

if [ $compare == true ]
then
    root -b -q -l "shadow_comp_v4.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/SB3-hit-no-scl.root\",\"${result_dir}/noSB2-hit-no-scl.root\"},{\"${result_dir}/noSB2-hit-no-scl.root\",\"${result_dir}/runs4-hit.root\",\"${result_dir}/runs23-hit.root\"},\"${result_dir}/comp-ns\")"
    root -b -q -l "visualize_comp_v3.cpp(\"${result_dir}/comp-ns.root\",\"${result_dir}/comp-vis-ns.root\")"

    root -b -q -l "shadow_comp_v4.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/SB3-hit.root\",\"${result_dir}/noSB2-hit.root\"},{\"${result_dir}/noSB2-hit.root\",\"${result_dir}/runs4-hit.root\",\"${result_dir}/runs23-hit.root\"},\"${result_dir}/comp\")"
    root -b -q -l "visualize_comp_v3.cpp(\"${result_dir}/comp.root\",\"${result_dir}/comp-vis.root\")"

    root -b -q -l "shadow_comp_v4.cpp({\"${result_dir}/runs23-hit.root\",\"${result_dir}/SB3-hit-new.root\",\"${result_dir}/noSB2-hit-new.root\"},{\"${result_dir}/noSB2-hit-new.root\",\"${result_dir}/runs4-hit.root\",\"${result_dir}/runs23-hit.root\"},\"${result_dir}/comp-new\")"
    root -b -q -l "visualize_comp_v3.cpp(\"${result_dir}/comp-new.root\",\"${result_dir}/comp-vis-new.root\")"

fi

if [ $percent_grid == true ]
then
    mkdir -p "${result_dir}/grid"
    #root -b -q -l "grid_make_v2.cpp(\"${result_dir}/all-cmb-hit.root\",\"${result_dir}/grid/ogb\",\"${result_dir}/comp-new_res/ogb.txt\",$num_list)"
    #root -b -q -l "grid_make_v2.cpp(\"${result_dir}/all-cmb-hit.root\",\"${result_dir}/grid/rb\",\"${result_dir}/comp-new_res/rb.txt\",$num_list)"
    #root -b -q -l "grid_make_v2.cpp(\"${result_dir}/all-cmb-hit.root\",\"${result_dir}/grid/rb2\",\"${result_dir}/comp-new_res/rb2.txt\",$num_list)"
    #root -b -q -l "grid_make_v2.cpp(\"${result_dir}/all-cmb-hit.root\",\"${result_dir}/grid/rb3\",\"${result_dir}/comp-new_res/rb3.txt\",$num_list)"

    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid/ogb\",\"${result_dir}/comp-new_res/ogb.txt\",$num_list)"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid/rb\",\"${result_dir}/comp-new_res/rb.txt\",$num_list)"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid/rb2\",\"${result_dir}/comp-new_res/rb2.txt\",$num_list)"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid/rb3\",\"${result_dir}/comp-new_res/rb3.txt\",$num_list)"

    mkdir -p "${result_dir}/grid-rnd"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid-rnd/ogb\",\"${result_dir}/comp-new_res/ogb_rnd.txt\",$num_list,\"hHitPattern_rnd\")"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid-rnd/rb\",\"${result_dir}/comp-new_res/rb_rnd.txt\",$num_list,\"hHitPattern_rnd\")"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid-rnd/rb2\",\"${result_dir}/comp-new_res/rb2_rnd.txt\",$num_list,\"hHitPattern_rnd\")"
    root -b -q -l "grid_make_v2.cpp(\"${result_dir}/runs-all-hit.root\",\"${result_dir}/grid-rnd/rb3\",\"${result_dir}/comp-new_res/rb3_rnd.txt\",$num_list,\"hHitPattern_rnd\")"

fi

if [ $mult_check == true ]
then
    mkdir -p "${result_dir}/mult"
    get_mpar par_mult_ogb "${result_dir}/grid/ogb.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult/ogb\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_ogb}})"

    get_mpar par_mult_rb "${result_dir}/grid/rb.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult/rb\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb}})"

    get_mpar par_mult_rb2 "${result_dir}/grid/rb2.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult/rb2\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb2}})"

    get_mpar par_mult_rb3 "${result_dir}/grid/rb3.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult/rb3\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb3}})"

        
    mkdir -p "${result_dir}/mult-rnd"
    get_mpar par_mult_ogb_rnd "${result_dir}/grid-rnd/ogb.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult-rnd/ogb\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_ogb_rnd}})"

    get_mpar par_mult_rb_rnd "${result_dir}/grid-rnd/rb.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult-rnd/rb\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb_rnd}})"

    get_mpar par_mult_rb2_rnd "${result_dir}/grid-rnd/rb2.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult-rnd/rb2\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb2_rnd}})"

    get_mpar par_mult_rb3_rnd "${result_dir}/grid-rnd/rb3.txt"
    root -b -q -l "shadow_mult_v2.cpp({$runs_all},\"${result_dir}/mult-rnd/rb3\",\"${par_name}.txt\",true,${use_veto},${a_mat},{${par_mult_rb3_rnd}})"

fi

