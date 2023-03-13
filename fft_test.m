CNV = VHDL_Conversion;
fft_lib = FFT;
hold on

N = 64;
L = 31;
dx = L/N;
step = 1;
x(1:L+1, 1) = ((-L/2 + 1):step:(L/2 + 1))';
x(1:L+1, 2) = zeros(L+1, 1);
f = cos(x).*exp((-power(x, 2)/25));

%fft from matlab
fft_mat = fft(f(:, 1));
ifft_mat = ifft(fft_mat);

%fft using traditional i as imaginary number
fft_ft2_test = fft_lib.ft2(f(:, 1));
ifft_ft2_test = fft_lib.ift2(fft_ft2_test);

%fft using one dimension for real and one for imaginary
fft_ft2_r_test = fft_lib.ft2_r(f);
ifft_ft2_r_test = fft_lib.ift2_r(fft_ft2_r_test);

%convert the to fixed point
vx(:, :, 1) = CNV.Double2FxP(f(:, 1), (L+1)/2, (L+1)/2);
vx(:, :, 2) = CNV.Double2FxP(f(:, 2), (L+1)/2, (L+1)/2);

%fft using fixed point and one dimension for real and one for imaginary
fft_ft2_r_bin_test_bin = fft_lib.ft2_r_bin(vx);
ifft_ft2_r_bin_test_bin = fft_lib.ift2_r_bin(fft_ft2_r_bin_test_bin);
[mw, nw, Lw, dw] = size(fft_ft2_r_bin_test_bin);
fft_ft2_r_bin_test_bink = zeros(mw, Lw);
for k = 1:nw
    fft_ft2_r_bin_test_bink(:, :) = fft_ft2_r_bin_test_bin(:, k, :, 1);
    fft_ft2_r_bin_test = CNV.FxP2Double(fft_ft2_r_bin_test_bink, (L+1)/2, (L+1)/2);
end

[miw, niw, Liw, diw] = size(ifft_ft2_r_bin_test_bin);
ifft_ft2_r_bin_test_bink = zeros(miw, Liw);
for k = 1:niw
    ifft_ft2_r_bin_test_bink(:, :) = ifft_ft2_r_bin_test_bin(:, k, :, 1);
    ifft_ft2_r_bin_test = CNV.FxP2Double(ifft_ft2_r_bin_test_bink, (L+1)/2, (L+1)/2);
    ifft_ft2_r_bin_test = ifft_ft2_r_bin_test/(L+1);
end


%plot(x(:, 1),fft_ft2_r_test(:, 1));
%plot(x(:, 1),fft_ft2_r_bin_test(:, 1));
plot(x(:, 1),ifft_ft2_r_bin_test(:, 1));
% plot(x(:, 1),fft_ft2_r_test(:, 1));
%plot(x(:, 1), f(:, 1));

% abs_r(:) = abs(ifft_ft2_r_test(:, 1, 1) - f(:, 1));
% abs_bin(:) = abs(ifft_ft2_r_bin_test(:, 1) - f(:, 1));
% plot(x(:, 1), abs_r(:));
% plot(x(:, 1), abs_bin(:));

err_bin = rms(ifft_ft2_r_bin_test(:, 1) - f(:, 1));
display(err_bin);


hold off





