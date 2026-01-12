#!/bin/bash

root -b -q -l "cut_load.cpp({\"E_s1_L0.C\",\"E_s2_L0.C\",\"E_s3_L0.C\",\"E_s4_L0.C\",\"E_s5_L0.C\"},\"E_L0.root\")"
root -b -q -l "cut_load.cpp({\"E_s1_L1.C\",\"E_s2_L1.C\",\"E_s3_L1.C\",\"E_s4_L1.C\",\"E_s5_L1.C\"},\"E_L1.root\")"
root -b -q -l "cut_load.cpp({\"E_s1_L2.C\",\"E_s2_L2.C\",\"E_s3_L2.C\",\"E_s4_L2.C\",\"E_s5_L2.C\"},\"E_L2.root\")"

root -b -q -l "cut_load.cpp({\"tt_s1_L0.C\",\"tt_s2_L0.C\",\"tt_s3_L0.C\",\"tt_s4_L0.C\",\"tt_s5_L0.C\"},\"tt_L0.root\")"
root -b -q -l "cut_load.cpp({\"tt_s1_L1.C\",\"tt_s2_L1.C\",\"tt_s3_L1.C\",\"tt_s4_L1.C\",\"tt_s5_L1.C\"},\"tt_L1.root\")"
root -b -q -l "cut_load.cpp({\"tt_s1_L2.C\",\"tt_s2_L2.C\",\"tt_s3_L2.C\",\"tt_s4_L2.C\",\"tt_s5_L2.C\"},\"tt_L2.root\")"
