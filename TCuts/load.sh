#!/bin/bash

root -b -q -l "cut_load.cpp({\"l0s1.C\",\"l0s2.C\",\"l0s3.C\",\"l0s4.C\",\"l0s5.C\"},\"l0.root\")"
root -b -q -l "cut_load.cpp({\"l1s1.C\",\"l1s2.C\",\"l1s3.C\",\"l1s4.C\",\"l1s5.C\"},\"l1.root\")"
root -b -q -l "cut_load.cpp({\"l2s1.C\",\"l2s2.C\",\"l2s3.C\",\"l2s4.C\",\"l2s5.C\"},\"l2.root\")"

root -b -q -l "cut_load.cpp({\"cl0s0.C\",\"cl0s1.C\",\"cl0s2.C\",\"cl0s3.C\",\"cl0s4.C\",\"cl0s5.C\"},\"cal_l0.root\")"
