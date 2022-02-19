function L1 = sos_construct_sosdec
pvar x1 x2
L_au = 6;
[L1,L1_Q] = sosdecvar('L1_w',monomials([x1;x2],0:L_au/2));
end