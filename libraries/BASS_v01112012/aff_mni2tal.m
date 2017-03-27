function outpoint = aff_mni2tal(inpoint)
Tfrm = [0.88 0 0 -0.8;0 0.97 0 -3.32; 0 0.05 0.88 -0.44;0 0 0 1];
tmp = Tfrm * [inpoint 1] ';
outpoint = tmp(1:3)';