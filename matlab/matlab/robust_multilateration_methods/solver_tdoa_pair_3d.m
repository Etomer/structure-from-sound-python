function sols = solver_tdoa_pair_3d(data)
[C0,C1] = setup_elimination_template(data);
C1 = C0 \ C1;
RR = [-C1(end-2:end,:);eye(8)];
AM_ind = [9,6,1,8,2,10,11,3];
AM = RR(AM_ind,:);
[V,D] = eig(AM);
V = V ./ (ones(size(V,1),1)*V(1,:));
sols(1,:) = V(2,:);
sols(2,:) = V(4,:);
sols(3,:) = diag(D).';

% Action =  z
% Quotient ring basis (V) = 1,x,x*z,y,y*z,z,z^2,z^3,
% Available monomials (RR*V) = x*z^2,y*z^2,z^4,1,x,x*z,y,y*z,z,z^2,z^3,
function [coeffs] = compute_coeffs(data)
coeffs(1) = 4*data(1)^2 - 8*data(1)*data(4) + 4*data(4)^2 - 4*data(19);
coeffs(2) = 8*data(1)*data(2) - 8*data(2)*data(4) - 8*data(1)*data(5) + 8*data(4)*data(5);
coeffs(3) = 4*data(2)^2 - 8*data(2)*data(5) + 4*data(5)^2 - 4*data(19);
coeffs(4) = 8*data(1)*data(3) - 8*data(3)*data(4) - 8*data(1)*data(6) + 8*data(4)*data(6);
coeffs(5) = 8*data(2)*data(3) - 8*data(3)*data(5) - 8*data(2)*data(6) + 8*data(5)*data(6);
coeffs(6) = 4*data(3)^2 - 8*data(3)*data(6) + 4*data(6)^2 - 4*data(19);
coeffs(7) = -4*data(1)^3 - 4*data(1)*data(2)^2 - 4*data(1)*data(3)^2 + 4*data(1)^2*data(4) + 4*data(2)^2*data(4) + 4*data(3)^2*data(4) + 4*data(1)*data(4)^2 - 4*data(4)^3 + 4*data(1)*data(5)^2 - 4*data(4)*data(5)^2 + 4*data(1)*data(6)^2 - 4*data(4)*data(6)^2 + 4*data(1)*data(19) + 4*data(4)*data(19);
coeffs(8) = -4*data(1)^2*data(2) - 4*data(2)^3 - 4*data(2)*data(3)^2 + 4*data(2)*data(4)^2 + 4*data(1)^2*data(5) + 4*data(2)^2*data(5) + 4*data(3)^2*data(5) - 4*data(4)^2*data(5) + 4*data(2)*data(5)^2 - 4*data(5)^3 + 4*data(2)*data(6)^2 - 4*data(5)*data(6)^2 + 4*data(2)*data(19) + 4*data(5)*data(19);
coeffs(9) = -4*data(1)^2*data(3) - 4*data(2)^2*data(3) - 4*data(3)^3 + 4*data(3)*data(4)^2 + 4*data(3)*data(5)^2 + 4*data(1)^2*data(6) + 4*data(2)^2*data(6) + 4*data(3)^2*data(6) - 4*data(4)^2*data(6) - 4*data(5)^2*data(6) + 4*data(3)*data(6)^2 - 4*data(6)^3 + 4*data(3)*data(19) + 4*data(6)*data(19);
coeffs(10) = data(1)^4 + 2*data(1)^2*data(2)^2 + data(2)^4 + 2*data(1)^2*data(3)^2 + 2*data(2)^2*data(3)^2 + data(3)^4 - 2*data(1)^2*data(4)^2 - 2*data(2)^2*data(4)^2 - 2*data(3)^2*data(4)^2 + data(4)^4 - 2*data(1)^2*data(5)^2 - 2*data(2)^2*data(5)^2 - 2*data(3)^2*data(5)^2 + 2*data(4)^2*data(5)^2 + data(5)^4 - 2*data(1)^2*data(6)^2 - 2*data(2)^2*data(6)^2 - 2*data(3)^2*data(6)^2 + 2*data(4)^2*data(6)^2 + 2*data(5)^2*data(6)^2 + data(6)^4 - 2*data(1)^2*data(19) - 2*data(2)^2*data(19) - 2*data(3)^2*data(19) - 2*data(4)^2*data(19) - 2*data(5)^2*data(19) - 2*data(6)^2*data(19) + data(19)^2;
coeffs(11) = 4*data(7)^2 - 8*data(7)*data(10) + 4*data(10)^2 - 4*data(20);
coeffs(12) = 8*data(7)*data(8) - 8*data(8)*data(10) - 8*data(7)*data(11) + 8*data(10)*data(11);
coeffs(13) = 4*data(8)^2 - 8*data(8)*data(11) + 4*data(11)^2 - 4*data(20);
coeffs(14) = 8*data(7)*data(9) - 8*data(9)*data(10) - 8*data(7)*data(12) + 8*data(10)*data(12);
coeffs(15) = 8*data(8)*data(9) - 8*data(9)*data(11) - 8*data(8)*data(12) + 8*data(11)*data(12);
coeffs(16) = 4*data(9)^2 - 8*data(9)*data(12) + 4*data(12)^2 - 4*data(20);
coeffs(17) = -4*data(7)^3 - 4*data(7)*data(8)^2 - 4*data(7)*data(9)^2 + 4*data(7)^2*data(10) + 4*data(8)^2*data(10) + 4*data(9)^2*data(10) + 4*data(7)*data(10)^2 - 4*data(10)^3 + 4*data(7)*data(11)^2 - 4*data(10)*data(11)^2 + 4*data(7)*data(12)^2 - 4*data(10)*data(12)^2 + 4*data(7)*data(20) + 4*data(10)*data(20);
coeffs(18) = -4*data(7)^2*data(8) - 4*data(8)^3 - 4*data(8)*data(9)^2 + 4*data(8)*data(10)^2 + 4*data(7)^2*data(11) + 4*data(8)^2*data(11) + 4*data(9)^2*data(11) - 4*data(10)^2*data(11) + 4*data(8)*data(11)^2 - 4*data(11)^3 + 4*data(8)*data(12)^2 - 4*data(11)*data(12)^2 + 4*data(8)*data(20) + 4*data(11)*data(20);
coeffs(19) = -4*data(7)^2*data(9) - 4*data(8)^2*data(9) - 4*data(9)^3 + 4*data(9)*data(10)^2 + 4*data(9)*data(11)^2 + 4*data(7)^2*data(12) + 4*data(8)^2*data(12) + 4*data(9)^2*data(12) - 4*data(10)^2*data(12) - 4*data(11)^2*data(12) + 4*data(9)*data(12)^2 - 4*data(12)^3 + 4*data(9)*data(20) + 4*data(12)*data(20);
coeffs(20) = data(7)^4 + 2*data(7)^2*data(8)^2 + data(8)^4 + 2*data(7)^2*data(9)^2 + 2*data(8)^2*data(9)^2 + data(9)^4 - 2*data(7)^2*data(10)^2 - 2*data(8)^2*data(10)^2 - 2*data(9)^2*data(10)^2 + data(10)^4 - 2*data(7)^2*data(11)^2 - 2*data(8)^2*data(11)^2 - 2*data(9)^2*data(11)^2 + 2*data(10)^2*data(11)^2 + data(11)^4 - 2*data(7)^2*data(12)^2 - 2*data(8)^2*data(12)^2 - 2*data(9)^2*data(12)^2 + 2*data(10)^2*data(12)^2 + 2*data(11)^2*data(12)^2 + data(12)^4 - 2*data(7)^2*data(20) - 2*data(8)^2*data(20) - 2*data(9)^2*data(20) - 2*data(10)^2*data(20) - 2*data(11)^2*data(20) - 2*data(12)^2*data(20) + data(20)^2;
coeffs(21) = 4*data(13)^2 - 8*data(13)*data(16) + 4*data(16)^2 - 4*data(21);
coeffs(22) = 8*data(13)*data(14) - 8*data(14)*data(16) - 8*data(13)*data(17) + 8*data(16)*data(17);
coeffs(23) = 4*data(14)^2 - 8*data(14)*data(17) + 4*data(17)^2 - 4*data(21);
coeffs(24) = 8*data(13)*data(15) - 8*data(15)*data(16) - 8*data(13)*data(18) + 8*data(16)*data(18);
coeffs(25) = 8*data(14)*data(15) - 8*data(15)*data(17) - 8*data(14)*data(18) + 8*data(17)*data(18);
coeffs(26) = 4*data(15)^2 - 8*data(15)*data(18) + 4*data(18)^2 - 4*data(21);
coeffs(27) = -4*data(13)^3 - 4*data(13)*data(14)^2 - 4*data(13)*data(15)^2 + 4*data(13)^2*data(16) + 4*data(14)^2*data(16) + 4*data(15)^2*data(16) + 4*data(13)*data(16)^2 - 4*data(16)^3 + 4*data(13)*data(17)^2 - 4*data(16)*data(17)^2 + 4*data(13)*data(18)^2 - 4*data(16)*data(18)^2 + 4*data(13)*data(21) + 4*data(16)*data(21);
coeffs(28) = -4*data(13)^2*data(14) - 4*data(14)^3 - 4*data(14)*data(15)^2 + 4*data(14)*data(16)^2 + 4*data(13)^2*data(17) + 4*data(14)^2*data(17) + 4*data(15)^2*data(17) - 4*data(16)^2*data(17) + 4*data(14)*data(17)^2 - 4*data(17)^3 + 4*data(14)*data(18)^2 - 4*data(17)*data(18)^2 + 4*data(14)*data(21) + 4*data(17)*data(21);
coeffs(29) = -4*data(13)^2*data(15) - 4*data(14)^2*data(15) - 4*data(15)^3 + 4*data(15)*data(16)^2 + 4*data(15)*data(17)^2 + 4*data(13)^2*data(18) + 4*data(14)^2*data(18) + 4*data(15)^2*data(18) - 4*data(16)^2*data(18) - 4*data(17)^2*data(18) + 4*data(15)*data(18)^2 - 4*data(18)^3 + 4*data(15)*data(21) + 4*data(18)*data(21);
coeffs(30) = data(13)^4 + 2*data(13)^2*data(14)^2 + data(14)^4 + 2*data(13)^2*data(15)^2 + 2*data(14)^2*data(15)^2 + data(15)^4 - 2*data(13)^2*data(16)^2 - 2*data(14)^2*data(16)^2 - 2*data(15)^2*data(16)^2 + data(16)^4 - 2*data(13)^2*data(17)^2 - 2*data(14)^2*data(17)^2 - 2*data(15)^2*data(17)^2 + 2*data(16)^2*data(17)^2 + data(17)^4 - 2*data(13)^2*data(18)^2 - 2*data(14)^2*data(18)^2 - 2*data(15)^2*data(18)^2 + 2*data(16)^2*data(18)^2 + 2*data(17)^2*data(18)^2 + data(18)^4 - 2*data(13)^2*data(21) - 2*data(14)^2*data(21) - 2*data(15)^2*data(21) - 2*data(16)^2*data(21) - 2*data(17)^2*data(21) - 2*data(18)^2*data(21) + data(21)^2;
function [C0,C1] = setup_elimination_template(data)
[coeffs] = compute_coeffs(data);
coeffs0_ind = [1,11,2,1,11,12,21,3,2,12,13,22,3,13,23,1,11,21,4,14,2,1,11,12,21,22,5,4,14,15,3,2,12,13,22,24,23,5,15,3,13,23,25,4,14,11,1,21,24,6,...
16,5,4,14,15,12,2,22,24,25,6,16,5,15,13,3,23,25,26,6,16,14,4,24,26,6,16,15,5,25,26,1,11,21,7,17,2,12,1,21,11,22,8,7,17,18,27,3,13,2,...
22,12,23,8,18,28,3,23,13,7,17,4,14,21,11,1,24,27,9,19,8,7,17,18,27,5,15,4,24,22,12,2,14,25,28,9,19,8,18,28,29,5,25,23,13,3,15,7,17,21,...
11,1,27,10,20,8,18,7,27,22,12,2,17,28,10,20,30,8,28,23,13,3,18,9,19,17,7,27,6,16,24,14,4,26,29,9,19,18,8,28,29,6,26,25,15,5,16,16,6,26];
coeffs1_ind = [30,20,10,10,20,27,17,7,30,10,20,9,19,27,17,7,24,14,4,29,30,10,30,28,18,8,20,10,20,30,9,29,28,18,8,25,15,5,19,30,20,10,29,19,9,20,10,30,29,19,...
9,26,16,6,19,9,29,26,16,6];
C0_ind = [1,4,27,28,29,30,39,53,54,55,56,65,80,81,91,109,112,130,131,134,135,136,137,138,142,156,157,158,159,160,161,162,163,164,168,169,182,184,185,188,189,194,195,213,216,217,218,219,234,235,...
238,239,240,241,242,243,244,245,246,260,262,263,266,267,269,270,271,272,273,291,294,295,296,297,312,318,319,321,322,323,324,352,353,363,365,368,378,379,380,381,388,389,391,392,393,394,403,404,405,406,...
407,414,415,418,419,429,432,433,440,447,450,456,457,460,461,462,467,468,469,472,473,474,475,476,480,482,483,484,485,486,487,488,492,493,494,496,497,500,501,506,507,510,511,512,513,514,518,534,535,541,...
542,543,545,547,550,560,561,562,563,567,568,569,570,571,574,575,585,588,589,593,594,595,596,603,606,607,608,609,612,613,616,617,618,623,624,630,631,633,634,635,636,640,641,642,643,644,648,659,660,661];
C1_ind = [21,22,23,40,41,47,48,49,51,57,60,66,67,70,71,72,73,74,75,77,78,94,95,99,100,101,102,110,111,116,120,121,122,123,124,125,126,127,128,148,149,150,151,152,153,165,166,167,174,175,...
176,177,178,179,191,192,193,200,201,202];
C0 = zeros(26,26);
C1 = zeros(26,8);
C0(C0_ind) = coeffs(coeffs0_ind);
C1(C1_ind) = coeffs(coeffs1_ind);

