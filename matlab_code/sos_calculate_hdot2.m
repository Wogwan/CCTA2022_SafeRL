function hdot = sos_calculate_hdot2(h,f1,f21,u2)
pvar x1 x2
hdot1 = jacobian(h, x1)*(f1);
hdot2 = jacobian(h, x2)*(f21+u2);
hdot = hdot1 + hdot2;
end