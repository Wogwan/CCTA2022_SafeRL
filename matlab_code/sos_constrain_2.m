function pconstr_22 = sos_constrain_2(hdot,L1,h)
pconstr_22 = hdot - L1*h >= 0;
end