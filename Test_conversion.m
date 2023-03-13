CNV = VHDL_Conversion;
fft_lib = FFT;

L = 6;
x1 = 1;
x2 = -1;
x3 = 4;

b1 = [1, 0, 0, 0];
b2 = [1, 1, 1, 1];
b3 = [1, 0, 0, 1];

y1a = CNV.bin2unsigned(b1);
y1b = CNV.bin2signed_signal(b1);
y1c = CNV.bin2signed_2s(b1);
%display([' unsigned ' num2str(y1a) ' signed signal ' num2str(y1b) ' signed 2s ' num2str(y1c)]);

y2a = CNV.bin2unsigned(b2);
y2b = CNV.bin2signed_signal(b2);
y2c = CNV.bin2signed_2s(b2);
%display([' unsigned ' num2str(y2a) ' signed signal ' num2str(y2b) ' signed 2s ' num2str(y2c)]);

y3a = CNV.bin2unsigned(b3);
y3b = CNV.bin2signed_signal(b3);
y3c = CNV.bin2signed_2s(b3);
%display([' unsigned ' num2str(y3a) ' signed signal ' num2str(y3b) ' signed 2s ' num2str(y3c)]);

y4 = CNV.unsigned2bin(x1, L);
y5 = CNV.unsigned2bin(x3, L);
%display([' 1 ' CNV.bits2str(y4) ' 4 ' CNV.bits2str(y5)]);

y6 = CNV.signed2bin_signal(x2, L);
y7 = CNV.signed2bin_2s(x2, L);
%display([' -1, signal ' CNV.bits2str(y6) ' -1, 2s ' CNV.bits2str(y7)]);

xd = 1.62;
vxd = CNV.Double2FxP(xd, 2, 6);
xd_e = CNV.FxP2Double(vxd, 2, 6);

%display(['Double: x = ' num2str(xd) ' Via FXP: x = ' num2str(xd_e)]);
display(vxd);

vx1 = [b1; b2; b3];
% vy1 = CNV.bin2unsigned_vector(vx1);
% vy2 = CNV.bin2signed_signal_vector(vx1);
% vy3 = CNV.bin2signed_2s_vector(vx1);
vy1x = CNV.bin2unsigned(vx1);
vy2x = CNV.bin2signed_signal(vx1);
vy3x = CNV.bin2signed_2s(vx1);

vx2 = [x1; x3; x1];
% vy4 = CNV.unsigned2bin_vector(vx2, L);
vy4x = CNV.unsigned2bin(vx2, L);

vx3 = [x1; x2; x3];
% vy5 = CNV.signed2bin_signal_vector(vx3, L);
% vy6 = CNV.signed2bin_2s_vector(vx3, L);
vy5x = CNV.signed2bin_signal(vx3, L);
vy6x = CNV.signed2bin_2s(vx3, L);

vx4 = [xd; xd; xd];
% vy7 = CNV.Double2FxP_vector(vx4, 2, 6);
% vy8 = CNV.FxP2Double_vector(vy7, 2, 6); 
vy7x = CNV.Double2FxP(vx4, 2, 6);
vy8x = CNV.FxP2Double(vy7x, 2, 6); 

vy9 = CNV.lshift(vx1, 1);
vy10 = CNV.rshift(vx1, 1);
vy11 = CNV.invert(vx1);
vy12 = CNV.sum_signed(vx1, vx1);
vy13 = CNV.minus(vx1, vx1);
vy14 = CNV.twoscomplement_vector(vx1);
% vy15 = CNV.bits2str_vector(vx1);
% 
% btest1 = [1, 1, 1, 0, 0, 0, 0];
% btest2 = [0, 1, 1, 1, 0, 0, 0];
% vytest = CNV.sum_unsigned(btest1, btest2);
% 
% vy16 = CNV.binarymultiplication(vx1, vx1);
% 
y11 = CNV.bin2signed_signal(vy11);
y11x = -CNV.bin2signed_signal(vx1);
y12 = CNV.bin2signed_signal(vy12);
y12x = CNV.bin2signed_signal(vx1) + CNV.bin2signed_signal(vx1);
y13 = CNV.bin2signed_signal(vy13);
y13x = CNV.bin2signed_signal(vx1) - CNV.bin2signed_signal(vx1);

% y16 = CNV.bin2signed_signal(vy16);
% y16x = CNV.bin2signed_signal(vx1) .* CNV.bin2signed_signal(vx1);
a0 = 0;
a1 = 3;
a2 = -1;
a3 = 1.5;
ba0 = CNV.signed2bin_signal(a0, 8);
ba1 = CNV.signed2bin_signal(a1, 8);
ba2 = CNV.signed2bin_signal(a2, 8);
fa0 = CNV.Double2FxP(a0, 4, 4);
fa1 = CNV.Double2FxP(a1, 4, 4);
fa2 = CNV.Double2FxP(a2, 4, 4);
fa3 = CNV.Double2FxP(a3, 4, 4);
fa1x = CNV.Double2FxP(a1, 16, 16);
fa2x = CNV.Double2FxP(a2, 16, 16);
fa3x = CNV.Double2FxP(a3, 16, 16);

cai(1, :, 1) = ba0;
cai(1, :, 2) = ba2;
ca1(1, :, 1) = ba1;
ca1(1, :, 2) = ba0;
ca2(1, :, 1) = ba2;
ca2(1, :, 2) = ba0;
cf1(1, :, 1) = fa1;
cf1(1, :, 2) = ba0; 
cf2(1, :, 1) = fa2;
cf2(1, :, 2) = ba0;
cf3(1, :, 1) = fa3;
cf3(1, :, 2) = ba0;
cai1(1, :, 1) = ba1;
cai1(1, :, 2) = ba1;
cai2(1, :, 1) = ba2;
cai2(1, :, 2) = ba2;
cai3(1, :, 1) = fa2;
cai3(1, :, 2) = fa2;
cai4(1, :, 1) = fa3;
cai4(1, :, 2) = fa3;

cfft1(1, 1, :, 1) = fa3x;
cfft1(1, 1, :, 2) = fa3x;
cfft2(1, 1, :, 1) = fa2x;
cfft2(1, 1, :, 2) = fa2x;

vyb12 = CNV.minus(ca1, ca2);
vyb12ss = CNV.sum_signedx(ca1, ca2);
vyb12su = CNV.sum_unsigned(ca1, ca2);
vyf12 = CNV.minus(cf1, cf2);
vyf13 = CNV.minus(cf1, cf3);
vyf12x = CNV.minusx(cf1, cf2);
vyf13x = CNV.minusx(cf1, cf3);
vymul3 = CNV.binarymultiplication(fa1x, fa2x);
vymul = CNV.complex_binary_multiplication(ca1, cai);
vymul1 = CNV.complex_binary_multiplication(cai1, cai2);
vymul2 = CNV.complex_binary_multiplication(cai3, cai4);
vysumfft = fft_lib.sum_fft(cfft1, cfft2);
vyminusfft = fft_lib.minus_fft(cfft1, cfft2);
vymulfft = fft_lib.complex_binary_multiplication_fft(cfft1, cfft2);
vymulfft2 = fft_lib.complex_binary_multiplication_fft(cfft1, cfft1);

yb12 = CNV.bin2signed_signal(vyb12(1, :, 1));
yb12ss = CNV.bin2signed_signal(vyb12ss(1, :, 1));
yb12su = CNV.bin2signed_signal(vyb12su(1, :, 1));
yf12 = CNV.FxP2Double(vyf12(1, :, 1), 5, 4);
yf13 = CNV.FxP2Double(vyf13(1, :, 1), 5, 4);
yf12x = CNV.FxP2Double(vyf12x(1, :, 1), 4, 4);
yf13x = CNV.FxP2Double(vyf13x(1, :, 1), 4, 4);
ymul = CNV.bin2signed_signal(vymul(1, :, 2));
ymul1i = CNV.bin2signed_signal(vymul1(1, :, 2));
ymul1r = CNV.bin2signed_signal(vymul1(1, :, 1));
vymul2x = CNV.reduce_size(vymul2, 8);
ymul2i = CNV.FxP2Double(vymul2(1, :, 2), 8, 7);
vymul3x = CNV.reduce_size(vymul3, 8);
ymul3 = CNV.FxP2Double(vymul3x(1, :), 4, 4);
ymul2r = CNV.FxP2Double(vymul2(1, :, 1), 8, 7);
ymul2xi = CNV.FxP2Double(vymul2x(1, :, 2), 4, 4);
ymul2xr = CNV.FxP2Double(vymul2x(1, :, 1), 4, 4);

ysumfftx(1, :, :) = vysumfft(1, 1, :, :);
yminusfftx(1, :, :) = vyminusfft(1, 1, :, :);
ymulfftx(1, :, :) = vymulfft(1, 1, :, :);
ymulfft2x(1, :, :) = vymulfft2(1, 1, :, :);

ysumfftr = CNV.FxP2Double(ysumfftx(1, :, 1), 16, 16);
ysumffti = CNV.FxP2Double(ysumfftx(1, :, 2), 16, 16);
yminusfftr = CNV.FxP2Double(yminusfftx(1, :, 1), 16, 16);
yminusffti = CNV.FxP2Double(yminusfftx(1, :, 2), 16, 16);
ymulfftr = CNV.FxP2Double(ymulfftx(1, :, 1), 16, 16);
ymulffti = CNV.FxP2Double(ymulfftx(1, :, 2), 16, 16);
ymulfft2r = CNV.FxP2Double(ymulfft2x(1, :, 1), 16, 16);
ymulfft2i = CNV.FxP2Double(ymulfft2x(1, :, 2), 16, 16);







z1 = [0, 0, 0, 0, 1, 1, 1, 1, 1];
z1x = CNV.reduce_size(z1, 4);