#ifndef MYCONST_H
#define MYCONST_H
#include <cmath>
#include <iostream>
#include <list>
#include <vector>
#include <unordered_map>
#include <tuple>
#include <string>
#include <fstream>
#include <regex>
#include <algorithm>


//particle masses in MeV/c^2 (from NIST reference)
double pro_m {938.272};
double neu_m {939.565};
double deu_m {1875.612};
double tri_m {2808.921};
double hel_m {2808.391};
double alp_m {3727.379};

//distance to the detector in meters
double det_dist {4.0};

//time offset in ns
double t_os {16.0};

//punch through energies in BC400 scintillator (from Lisecute++) in MeV
double pro_pt {47.971};
double deu_pt {65.526};
double tri_pt {78.672};
double hel_pt {168.557};
double alp_pt {191.353};

std::vector<double> pt_ens {pro_pt,deu_pt,tri_pt,hel_pt,alp_pt};

double pro_pt_err {1.9639};
double deu_pt_err {1.7098};
double tri_pt_err {1.5816};
double hel_pt_err {2.3814};
double alp_pt_err {2.2796};

std::vector<double> pt_errs {pro_pt_err,deu_pt_err,tri_pt_err,hel_pt_err,alp_pt_err};

double cos_mu {8.0};
double mu_tot0 {21.0};
double mu_tot1 {19.5};
double mu_tot2 {20.5};

//BC408 can detect x rays of energies <100keV (around 10% E res around here)
double En_min {0.1};
double En_min_err {0.01};

double num_ent = 13215485.0;

//punch though points (estimated by eye)
std::vector<double> pro_0 {25.4,25.7};
std::vector<double> deu_0 {32.8,26.4};
std::vector<double> tri_0 {37.8,27.1};
std::vector<double> hel_0 {21.0,30.0};
std::vector<double> alp_0 {23.8,30.8};

std::vector<double> pro_1 {21.4,24.2};
std::vector<double> deu_1 {27.8,25.2};
std::vector<double> tri_1 {32.4,25.9};
std::vector<double> hel_1 {17.8,28.4};
std::vector<double> alp_1 {20.0,29.2};

std::vector<double> pro_2 {19.0,24.9};
std::vector<double> deu_2 {25.2,25.9};
std::vector<double> tri_2 {29.4,26.6};
std::vector<double> hel_2 {15.4,29.3};
std::vector<double> alp_2 {17.6,30.4};

std::vector<double> L0_tof {pro_0[0],deu_0[0],tri_0[0],hel_0[0],alp_0[0]};
std::vector<double> L1_tof {pro_1[0],deu_1[0],tri_1[0],hel_1[0],alp_1[0]};
std::vector<double> L2_tof {pro_2[0],deu_2[0],tri_2[0],hel_2[0],alp_2[0]};

std::vector<double> L0_tot {pro_0[1],deu_0[1],tri_0[1],hel_0[1],alp_0[1]};
std::vector<double> L1_tot {pro_1[1],deu_1[1],tri_1[1],hel_1[1],alp_1[1]};
std::vector<double> L2_tot {pro_2[1],deu_2[1],tri_2[1],hel_2[1],alp_2[1]};


/* std::vector<double> L0_tof {pro_tof_0,deu_tof_0,tri_tof_0,hel_tof_0,alp_tof_0};
std::vector<double> L1_tof {pro_tof_1,deu_tof_1,tri_tof_1,hel_tof_1,alp_tof_1};
std::vector<double> L2_tof {pro_tof_2,deu_tof_2,tri_tof_2,hel_tof_2,alp_tof_2};
 */

// (m/ns)^2 = (J/kg)/(10^-9)^2 = (6241509074461 MeV / 6.0221366516752E+26 u) * (10^18) = 10364.27673 MeV/u
//conversion factor
double cnv_fctr {10364.27673};

////////
//limit loaded histograms
///////

int hstart1 {2};
int hstop1 {23};
int hbase1 {14};
int hstart2 {24};
int hstop2 {47};
int hbase2 {38};
int hstart3 {48};
int hstop3 {71};
int hbase3 {62};


std::vector<int> hstarts {hstart1,hstart2,hstart3};
std::vector<int> hstops {hstop1,hstop2,hstop3};
std::vector<int> hbases {hbase1,hbase2,hbase3};


//sections
//{1,2,3,4,5,6,7};
std::vector<int> sec_omit {1,6,7}; //change which sections to omit
std::vector<int> sec_omit_r {1};
std::vector<int> sec_omit_e {};

int max_deg {2}; //for fitting polynomials
int min_deg {1};

//initial guess for a parameters (adjust to the parameters of the baseline detector)
double a1 {0.5}; 
double a2 {0.29};
double ka {0.25};

double E_A {0.02};
double E_B {-0.01};
double E_k {0.30};

/* double a1 {0.5}; 
double a2 {0.29};
double ka {0.15}; */

//for checking the improvement of statistics w/ the calibration method
double xc1 {20.0};
double xc2 {25.0};
double xc3 {28.0};
double xc4 {32.0};

double yc1 {20.0};
double yc2 {24.0};
double yc3 {26.0};
double yc4 {30.0};

//detectors blocked completely by shadow bars
std::vector<int> sb2 {6,41,54};
std::vector<int> sb3 {17,41,65};
std::vector<int> sb4 {17,30,65};

//positions of the shadow bars
std::vector<double> sb2_a {-184.0,181.0};
std::vector<double> sb2_b {-184.0,274.0};
std::vector<double> sb2_c {-278.0,274.0};
std::vector<double> sb2_d {-277.0,181.0};

double sb2_xmin {-277.5};
double sb2_xmax {-184.0};
double sb2_ymin {181.0};
double sb2_ymax {274.0};

std::vector<double> sb3_a {-181.0,-276.0};
std::vector<double> sb3_b {-181.0,-184.0};
std::vector<double> sb3_c {-280.0,-184.0};
std::vector<double> sb3_d {-280.0,-277.0};

double sb3_xmin {-280.0};
double sb3_xmax {-181.0};
double sb3_ymin {-276.5};
double sb3_ymax {-184.0};

std::vector<double> sb4_a {277.0,-279.0};
std::vector<double> sb4_b {276.0,-181.0};
std::vector<double> sb4_c {184.0,-181.0};
std::vector<double> sb4_d {185.0,-280.0};

double sb4_xmin {184.0};
double sb4_xmax {276.5};
double sb4_ymin {-279.5};
double sb4_ymax {-181.0};

//vw edges
std::vector<int> v0l0 {2,7};
std::vector<int> v1l0 {7,16};
std::vector<int> v2l0 {16,23};

std::vector<double> v0l1 {155.0,500.0};
std::vector<double> v1l1 {-185.0,195.0};
std::vector<double> v2l1 {-500.0,-165.0};

std::vector<int> v0l2 {48,55};
std::vector<int> v1l2 {55,64};
std::vector<int> v2l2 {64,71};


std::vector<int> v0l1d {26,30};
std::vector<int> v1l1d {30,39};
std::vector<int> v2l1d {39,56};

//tot threshold for wc
double my_tot_thresh_wc = 21.233;

//tot threshold for bg
double my_tot_thresh = 18.0;

//rebin factor
int rbf {2};


int nbins_hitpattern {200}; //default 200


double z_target_jun = -21.86;
double z_target_nov = -32.8;

double z_sig_jun = 2.15;
double z_sig_nov = 2.5;


std::vector<int> nov_runs = {1644,1645,1646,1647,1648,1649,1650,1651,1652,1653,1654,1655,1656,1657,1658,1659,1660,1661,1662,1663,1664,1665,1666,1667,1668,1669,1670,1671,1672,1673,1675,1676,1677,1679,1680,1681,1682,1683,1684,1685,1686,1687,1688,1689,1690,1691,1692,1693,1694,1695,1696,1697,1698,1699,1700,1701,1702,1703,1704,1705,1706,1707,1708,1709,1710,1711,1721,1729,1730,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,1745,1746,1747,1748,1750,1751,1752,1753,1754,1755,1756,1757,1758,1759,1760,1761,1762,1763,1764,1765,1766,1767,1768,1769,1770,1771,1772,1773,1774,1775,1776,1777,1778,1779,1780,1781,1782,1783,1784,1785,1786,1788,1790,1791,1792,1793,1794,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,};


#endif
