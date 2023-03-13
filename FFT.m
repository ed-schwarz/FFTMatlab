classdef FFT
    methods (Static)

        %% Calculated the fft of a complex vector, using i
        function Y = ft2(vx)
            vx = vx';
            [m, L, d] = size(vx);

            l=log2(L);
            p=ceil(l);
            Y=vx;
            N = 2^p;
            N2=N/2;
            seq = 0 : N2-1;
            W=(exp(-pi*sqrt(-1)/N2)).^seq;
            for L = 1 : p-1
                u=Y(:,1:N2, :);
                v=Y(:,N2+1:N, :);
                t=u+v;
                S=W.*(u-v);
                Y=[t ; S];
                U=W(:,1:2:N2, :);
                W=[U ;U];
                N=N2;
                N2=N2/2;
            end
            u=Y(:,1, :);
            v=Y(:,2, :);
            Y=[u+v;u-v];
        end

        %% Calculated the fft of a complex vector, separating complex and real
        %(:, :, 1) real, (:, :, 2) complex
        function Y = ft2_r(vx)
            CNV = VHDL_Conversion;
            [m, L, d] = size(vx);
            vs(1, :, 1) = vx(:, 1);
            vs(1, :, 2) = vx(:, 2);
            vt(:, :, 1) = vs(:, :, 1)';
            vt(:, :, 2) = vs(:, :, 2)';
            l=log2(m);
            p=ceil(l);
            Y=vs;
            [my, Ly, dy] = size(Y);
            N = 2^p;
            N2=N/2;

            for k = 1:N2
                W(my, k, 1) = cos(-pi/N2*(k-1));
                W(my, k, 2) = sin(-pi/N2*(k-1));
            end

            for L = 1 : p-1
                u=Y(:,1:N2, :);
                v=Y(:,N2+1:N, :);
                t=u+v;
                S=CNV.complexmultiplication_in(W, (u-v));
                Y=[t ; S];
                U=W(:,1:2:N2, :);
                W=[U ;U];
                N=N2;
                N2=N2/2;
            end
            u=Y(:,1, :);
            v=Y(:,2, :);
            Y=[u+v;u-v];
        end

        %% Calculated the fft of a complex binary vector, separating complex and real
        %(:, :, 1) real, (:, :, 2) complex
        function Y = ft2_r_bin(vx)
            CNV = VHDL_Conversion;
            [m, L, d] = size(vx);
            vs(1, :, :, 1) = vx(:, :, 1);
            vs(1, :, :, 2) = vx(:, :, 2);
            vt(1, :, :, 1) = vs(:, :, 1)';
            vt(1, :, :, 2) = vs(:, :, 2)';
            l=log2(m);
            p=ceil(l);
            Y=vs;
            [my, ky, Ly, dy] = size(Y);
            N = 2^p;
            N2=N/2;
            for k = 1:N2
                W(1, k, :, 1) = CNV.Double2FxP(cos(-pi/N2*(k-1)), Ly/2, Ly/2);
                W(1, k, :, 2) = CNV.Double2FxP(sin(-pi/N2*(k-1)), Ly/2, Ly/2);
            end

            for L = 1 : p-1
                u=Y(:,1:N2, :, :);
                v=Y(:,N2+1:N, :, :);
                t = FFT.sum_fft(u,v);
                t_min = FFT.minus_fft(u, v);
                S = FFT.complex_binary_multiplication_fft(W, t_min);
                Y=[t ; S];
                U=W(:,1:2:N2, :, :);
                W=[U ;U];
                N=N2;
                N2=N2/2;
            end
            u=Y(:,1, :, :);
            v=Y(:,2, :, :);
            t_sum = FFT.sum_fft(u,v);

            t_minus = FFT.minus_fft(u, v);

            Y=[t_sum; t_minus];
        end

        %% Calculates the inverse fft from ft2
        function Y = ift2(vx)
            x = vx.*exp(-1i*pi);
            Y = FFT.ft2(x);
        end

        %% Calculates the inverse fft from ft2_r
        function Y = ift2_r(vx)
            vx(:, 2) = -vx(:, 2);
            Y = FFT.ft2_r(vx);
        end

        %% Calculate the inverse fft from ft2_r_bin
        function Y = ift2_r_bin(vx)
            CNV = VHDL_Conversion;
            [m, n, L, d] = size(vx);
            x(:, :, :) = vx(:, 1, :, :);
            x(:, :, 2) = CNV.invert(x(:, :, 2));
            Y = FFT.ft2_r_bin(x);
        end

        %% Calculates the binary sum of complex binary vectors for the fft
        function vy = sum_fft(vx1, vx2)
            CNV = VHDL_Conversion;
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            m = max(m1,m2);
            n = max(n1,n2);
            L = max(L1,L2);
            d = max(d1,d2);
            vy = zeros(m , n, L, d);
            vxk1 = zeros(n, L, d);
            vxk2 = zeros(n, L, d);
            for k = 1:m
                vxk1(:, :, :) = vx1(k, :, :, :);
                vxk2(:, :, :) = vx2(k, :, :, :);
                vyk = CNV.sum_signedx(vxk1, vxk2);
                vy(k, :, 1:end, :) = vyk(:, 1:end, :);

            end
        end

        %% Calculates the binary difference of complex binary vectors for the fft
        function vy = minus_fft(vx1, vx2)
            CNV = VHDL_Conversion;
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            m = max(m1,m2);
            n = max(n1,n2);
            L = max(L1,L2);
            d = max(d1,d2);
            vy = zeros(m , n, L, d);
            vxk1 = zeros(n, L, d);
            vxk2 = zeros(n, L, d);
            for k = 1:m
                vxk1(:, :, :) = vx1(k, :, :, :);
                vxk2(:, :, :) = vx2(k, :, :, :);
                vyk = CNV.minusx(vxk1, vxk2);
                vy(k, :, 1:end, :) = vyk(:, 1:end, :);

            end
        end

        %% Calculates the complex binary multiplication of complex binary vectors for the fft
        function vy = complex_binary_multiplication_fft(vx1, vx2)
            CNV = VHDL_Conversion;
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            m = max(m1,m2);
            n = max(n1,n2);
            L = max(L1,L2);
            d = max(d1,d2);
            vy = zeros(m , n, L, d);
            vxk1 = zeros(n, L, d);
            vxk2 = zeros(n, L, d);
            vxk = zeros(n, L, d);
            for k = 1:m
                vxk1(:, :, :) = vx1(k, :, :, :);
                vxk2(:, :, :) = vx2(k, :, :, :);
                vyk = CNV.complex_binary_multiplication_FxP(vxk1, vxk2);
                vxk(:, :, 1) = CNV.reduce_size_x(vyk(:, :, 1), L);
                vxk(:, :, 2) = CNV.reduce_size_x(vyk(:, :, 2), L);
                vy(k, :, :, :) = vxk(:, :, :);

            end
        end

    end
end