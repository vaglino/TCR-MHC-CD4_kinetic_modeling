L1 = 3.700760241587129;
L2 = 1.914718101631063;
L3 = 2.4008921528083533;
L4 = 36.00905301970478;


% Model 2 vs 1
[ p2, Fstat2, df12, df22 ] = ftest(97, 1,3, L1,L2)

% Model 3 vs 1
[ p3, Fstat3, df13, df13 ] = ftest(97, 1,3, L1,L3)

% Model 4 vs 1
[ p4, Fstat4, df14, df24 ] = ftest(97, 1,3, L1,L4)
