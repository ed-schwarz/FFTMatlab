classdef VHDL_Conversion
    methods (Static)

        %%all binary numbers use the index 1 as the least significant bit

        %% Convert the binary to unsigned numner
        function x = bin2unsigned(vx, L)
            % vx is the binary number as vector, L is the size of it
            if nargin < 2
                [m, L, d] = size(vx);
            end

            l2 = pow2(0:L-1)';
            x = vx*l2;
        end

        %%Convert unsigned number to binary vector of length L
        function [vx, clipped] = unsigned2bin(x,L)
            %x is the unsigned number

            x = reshape(x, [], 1);
            M = size(x,1);
            vx = zeros(M,L) > 0;
            clipped_neg = x < 0;
            vx(~clipped_neg, : ) = dec2bin(x (~clipped_neg),L)> '0';

            clipped_pos =  x > (pow2(L)-1);
            vx(clipped_pos, : ) =  true;
            clipped = or(clipped_neg, clipped_pos);
            vx = flip(vx, 2);

        end


        %% Convert binary to signed
        function x = bin2signed_signal(vx)
            %with the most significant bit beeing the signal
            %vx is the binary number as vector
            [m, L, d] = size(vx);
            vy(:, 1:L-1) = vx(:, 1:L-1);
            x = VHDL_Conversion.bin2unsigned(vy);
            vs = ones(1, size(x, 1))';
            vs = vs - 2*vx(:, L);
            x = vs.*x;
        end

        %% Convert binary to signed
        function x = bin2signed_2s(vx)
            %with the binary as twos complement representation
            %vx is the binary number as vector
            [m, L, d] = size(vx);
            vy(:, 1:L-1) = vx(:, 1:L-1);
            x = VHDL_Conversion.bin2unsigned(vy);
            vs = 2^(L-1)*vx(:, L);
            x = x - vs;
        end

        %% Convert signed number to binary
        function [vx, clipped] = signed2bin_signal(x,L)
            %with the most significant bit beeing the signal
            %with x beeing the number and L th desired size of the binary
            y = x;
            y(x < 0) = 2^(L - 1) - x(x < 0);
            [vx, clipped] = VHDL_Conversion.unsigned2bin(y, L);

        end

        %% Convert signed number to binary
        function [vx, clipped] = signed2bin_2s(x,L)      
            %with the binary as twos complement representation
            %with x beeing the number and L th desired size of the binary
            n = 2^(L);
            y = x;
            y(x < 0) = n + x(x < 0);
            [vx, clipped] = VHDL_Conversion.unsigned2bin(y, L);
        end


        %% Converts binary vector to binary string
        function s = bits2str(vx)
            %with vx as the binary vector
            [m, L, d] = size(vx);
            s = char(1,L);
            vx = vx(end:-1:1);
            for i=1:L
                if vx(i)
                    s(i)= '1';
                else
                    s(i)= '0';
                end
            end
            s= s';
        end


        %% Converts a double to a fixed point binary
        function [x_fxp, overflow] = Double2FxP(x, I, F)
            %with I beeing the length of the Integer part
            %with F beeing the length of the Farctional part
            %with x as the double
            x_i = fix(x*pow2(F+1));
            [vx, overflow] = VHDL_Conversion.signed2bin_signal(x_i,I+F+1);
            x_fxp = vx(:, 2:end);
        end

        %% Converts a fixed point vetor to a double
        function x = FxP2Double(vx, I, F)  
            %with I beeing the length of the Integer part
            %with F beeing the length of the Farctional part
            %with vx as the binary vector
            [m, n, d] = size(vx);
            if n ~= (I+F)
                error('wrong length');
            end
            x = VHDL_Conversion.bin2signed_signal(vx)/pow2(F);
        end

        %% Perform left shift by m bits
        function vy = lshift(vx, m)
            %with vx as the binary vector
            vy = zeros(size(vx))> 0;
            vy (:,(m+1):end) = vx(:,1:(end-m));
        end

        %% Perform right  shift by m bits
        function vy = rshift(vx, m)
            %with vx as the binary vector
            vy = zeros(size(vx))> 0;
            vy (:,1:(end-1)) = vx(:,(1+m):end);
        end

        %% Execute a sum of two unsigned binary vectors (vx1 and vx2)
        function vy = sum_unsigned(vx1, vx2)   
            [m1, L1] = size(vx1);
            [m2, L2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);

            carry = zeros(m,L+1)>0;
            vytemp = carry;
            carry(:, 2:L+1) = and(vx1(:, 1:L), vx2(:, 1:L));

            vytemp(:, 1:L) = xor(vx1(:, 1:end),vx2(:, 1:end));
            carrytemp(:, 2:L+1) = and(vytemp(:, 1:L), carry(:, 1:L));
            vy = carry;
            vy(:, 1:end) = xor(xor(vytemp(:, 1:end), carry(:, 1:end)), carrytemp(:, 1:end));

        end


        function vy = sum_signed(vx1, vx2)
            %execute a sum of two signed binary vectors (vx1 and vx2)
            [m1, L1, d1] = size(vx1);
            [m2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            d = max(d1,d2);
            vs1 = zeros(m, L - L1, d);

            vx1 = [vx1, vs1];
            vx1 = logical(vx1);
            vs2 = zeros(m, L - L2, d);

            vx2 = [vx2, vs2];
            vx2 = logical(vx2);

            for i = 1:d

                vsub = xor(vx1(:, L, i), vx2(:, L, i));

                carry = zeros(m,L+1, d)>0;
                vytemp = carry;
                carrytemp = carry;
                carry(:, 2:L+1, i) = and(vx1(:, 1:L, i), vx2(:, 1:L, i));

                vytemp(:, 1:L-1, i) = xor(vx1(:, 1:end -1, i),vx2(:, 1:end -1, i));
                carrytemp(:, 2:L, i) = and(vytemp(:, 1:L-1, i), carry(:, 1:L-1, i));
                vy(:, :, i) = carry(:, :, i);
                vy(:, 1:end-1, i) = xor(xor(vytemp(:, 1:end-1, i), carry(:, 1:end-1, i)), carrytemp(:, 1:end - 1, i));


                vx1sub = vx1;
                vx1sub(vx1(:, L, i), 1:end-1, i) = ~vx1(vx1(:, L, i), 1:end-1, i);
                vx2sub = vx2;
                vx2sub(vx2(:, L, i), 1:end-1, i) = ~vx2(vx2(:, L, i), 1:end-1, i);

                carrysub = zeros(m,L+1, d)>0;
                vytempsub = carrysub;
                carrytempsub = carrysub;
                carrysub(:, 2:L+1, i) = and(vx1sub(:, 1:L, i), vx2sub(:, 1:L, i));

                vytempsub(:, 1:L, i) = xor(vx1sub(:, 1:end, i),vx2sub(:, 1:end, i));
                carrytempsub(:, 2:L, i) = and(vytempsub(:, 1:L-1, i), carrysub(:, 1:L-1, i));
                vysub(:, :, i) = carrysub(:, :, i);
                vysub(:, 1:end-1, i) = xor(xor(vytempsub(:, 1:end-1, i), carrysub(:, 1:end-1, i)), carrytempsub(:, 1:end-1, i));
                vysub(:, L+1, i) = and(vx1sub(:, L, i), vx2sub(:, L, i));
                vysub(~carrysub(:, L, i), 1:end, i) = ~vysub(~carrysub(:, L, i), 1:end, i);

                vy(vsub(:), :, i) = vysub(vsub(:), :, i);
            end

        end

        function vy = sum_unsigned_4d(vx1, vx2)
            %execute a sum of two unsigned binary vectors (vx1 and vx2)
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            n = max(n1,n2);
            d = max(d1,d2);

            carry = zeros(m, n,L+1, d)>0;
            vytemp = carry;
            carry(:, :, 2:L+1, :) = and(vx1(:, :, 1:L, :), vx2(:, :, 1:L, :));

            vytemp(:, :, 1:L, :) = xor(vx1(:, :, 1:end, :),vx2(:, :, 1:end, :));
            carrytemp(:, :, 2:L+1, :) = and(vytemp(:, :, 1:L, :), carry(:, :, 1:L, :));
            vy(:, :, 1:L, :) = carry(:, :, 1:L, :);
            vy(:, :, 1:end, :) = xor(xor(vytemp(:, :, 1:end-1, :), carry(:, :, 1:end-1, :)), carrytemp(:, :, 1:end-1, :));

        end

        %% Execute a sum of two signed binary vectors (vx1 and vx2)
        function [vy, clipped] = sum_signedx(vx1, vx2)
            [m1, L1, d1] = size(vx1);
            [m2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            d = max(d1,d2);
            vs1 = zeros(m, L - L1, d);

            vx1 = [vx1, vs1];
            vx1 = logical(vx1);
            vs2 = zeros(m, L - L2, d);

            vx2 = [vx2, vs2];
            vx2 = logical(vx2);
            carry = zeros(m,L +1, d)>0;

            carrysub = zeros(m,L+1, d)>0;

            for i = 1:d

                vsub = xor(vx1(:, L, i), vx2(:, L, i));


                vy(:, 1, i) = xor(vx1(:, 1, i), vx2(:, 1, i));
                carry(:, 2, i) = and(vx1(:, 1, i), vx2(:, 1, i));
                for l =2:L
                    vy(:, l, i) = xor(vx1(:, l, i), vx2(:, l, i));
                    vy(:, l, i) = xor(vy(:, l, i), carry(:, l, i));
                    carry(:, l+1, i) = or(and(vx1(:, l, i), vx2(:, l, i)), and(vx1(:, l, i), carry(:, l, i)));
                    carry(:, l+1, i) = or(and(vx2(:, l, i), carry(:, l, i)), carry(:, l+1, i));
                end

                vx1sub = vx1;
                vx1sub(vx1(:, L, i), 1:end-1, i) = ~vx1(vx1(:, L, i), 1:end-1, i);
                vx2sub = vx2;
                vx2sub(vx2(:, L, i), 1:end-1, i) = ~vx2(vx2(:, L, i), 1:end-1, i);

                vysub(:, 1, i) = xor(vx1sub(:, 1, i), vx2sub(:, 1, i));
                carrysub(:, 2, i) = and(vx1sub(:, 1, i), vx2sub(:, 1, i));
                for l =2:L
                    vysub(:, l, i) = xor(vx1sub(:, l, i), vx2sub(:, l, i));
                    vysub(:, l, i) = xor(vysub(:, l, i), carrysub(:, l, i));
                    carrysub(:, l+1, i) = or(and(vx1sub(:, l, i), vx2sub(:, l, i)), and(vx1sub(:, l, i), carrysub(:, l, i)));
                    carrysub(:, l+1, i) = or(and(vx2sub(:, l, i), carrysub(:, l, i)), carrysub(:, l+1, i));
                end

                vysub(~carrysub(:, L, i), 1:end-1, i) = ~vysub(~carrysub(:, L, i), 1:end-1, i);
                vy(carrysub(:, L, i), L, i) = or(vx1(carrysub(:, L, i), L, i), vx2(carrysub(:, L, i), L, i));
                vy(vsub(:), :, i) = vysub(vsub(:), :, i);
            end
        end


        function vy = sum_signed_4d(vx1, vx2)
            %execute a sum of two signed binary vectors (vx1 and vx2)
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            n = max(n1,n2);
            d = max(d1,d2);
            vs1 = zeros(m, n, L - L1, d);

            vx1 = [vx1; vs1];
            vx1 = logical(vx1);
            vs2 = zeros(m, n, L - L2, d);

            vx2 = [vx2; vs2];
            vx2 = logical(vx2);

            for i = 1:d
                for k =1:m
                    vsub = xor(vx1(k, :, L, i), vx2(k, :, L, i));

                    carry = zeros(m, n,L+1, d)>0;
                    vytemp = carry;
                    carrytemp = carry;
                    carry(k, :, 2:L+1, i) = and(vx1(k, :, 1:L, i), vx2(k, :, 1:L, i));

                    vytemp(k, :, 1:L-1, i) = xor(vx1(k, :, 1:end -1, i),vx2(k, :, 1:end -1, i));
                    carrytemp(k, :, 2:L, i) = and(vytemp(k, :, 1:L-1, i), carry(k, :, 1:L-1, i));
                    vy(k, :, 1:L, i) = carry(k, :, 1:end-1, i);
                    vy(k, :, 1:end-1, i) = xor(xor(vytemp(k, :, 1:end-2, i), carry(k, :, 1:end-2, i)), carrytemp(k, :, 1:end - 2, i));


                    vx1sub = vx1;
                    vx1sub(k, vx1(k, :, L, i), 1:end-1, i) = ~vx1(k, vx1(k, :, L, i), 1:end-1, i);
                    vx2sub = vx2;
                    vx2sub(k, vx2(k, :, L, i), 1:end-1, i) = ~vx2(k, vx2(k, :, L, i), 1:end-1, i);

                    carrysub = zeros(m, n, L+1, d)>0;
                    vytempsub = carrysub;
                    carrytempsub = carrysub;
                    carrysub(k, :, 2:L+1, i) = and(vx1sub(k, :, 1:L, i), vx2sub(k, :, 1:L, i));

                    vytempsub(k, :, 1:L, i) = xor(vx1sub(k, :, 1:end, i),vx2sub(k, :, 1:end, i));
                    carrytempsub(k, :, 2:L, i) = and(vytempsub(k, :, 1:L-1, i), carrysub(k, :, 1:L-1, i));
                    vysub(k, :, 1:L, i) = carrysub(k, :, 1:end-1, i);
                    vysub(k, :, 1:end-1, i) = xor(xor(vytempsub(k, :, 1:end-2, i), carrysub(k, :, 1:end-2, i)), carrytempsub(k, :, 1:end-2, i));
                    vysub(k, :, L, i) = and(vx1sub(k, :, L, i), vx2sub(k, :, L, i));
                    vysub(k, ~carrysub(k, :, L, i), 1:end, i) = ~vysub(k, ~carrysub(k, :, L, i), 1:end, i);

                    vy(k, vsub(:), :, i) = vysub(k, vsub(:), :, i);
                end
            end


        end

        function vy = minus(vx1, vx2)
            %execute a subtraction of two binary vectors (vx1 and vx2)
            vy = VHDL_Conversion.sum_signed(vx1, VHDL_Conversion.invert(vx2));
        end

        %% Execute a subtraction of two binary vectors (vx1 and vx2)
        function vy = minusx(vx1, vx2)
            vy = VHDL_Conversion.sum_signedx(vx1, VHDL_Conversion.invert(vx2));
        end

        function vy = minus_4d(vx1, vx2)
            %execute a subtraction of two binary vectors (vx1 and vx2)
            vy = VHDL_Conversion.sum_signed_4d(vx1, VHDL_Conversion.invert_4d(vx2));
        end

        %% Invert the signt of a binary vector
        function [vy, clipped] = invert(vx)
            %with vx as the binary vector
            [m, L, d] = size(vx);
            vx(:, L, :) = ~vx(:, L, :);
            vy = vx;
        end

        function [vy, clipped] = invert_4d(vx)
            %invert the signt of a binary vector
            %with vx as the binary vector
            [m, n, L, d] = size(vx);
            vx(:, :, L, :) = ~vx(:, :, L, :);
            vy = vx;
        end

        %% Builds the twos complement of a binary vector
        function vy = twoscomplement(vx)
            %with vx as the binary vector
            vy = not(vx);
            carry = true;
            [m, L, d] = size(vx);
            for i = 1+(0:L-1)
                vy(i) = xor(vy(i), carry);
                carry = not(vy(i));
                if (~carry)
                    return;
                end
            end
        end

        %% Multiply two numbers with real part at row 1 and complex at row 2
        function vy = complexmultiplication(vx1, vx2)
            % z = a + j*b
            % z*w = a*c - b*d + j*(a*d + b*c)

            vy(1:end, 1) = vx1(1:end, 1)*vx2(1:end, 1) - vx1(1:end, 2)*vx2(1:end, 2);
            vy(1:end, 2) = vx1(1:end, 2)*vx2(1:end, 1) + vx1(1:end, 1)*vx2(1:end, 2);
        end

        %% Multiply two numbers with real part at row 1 and complex at row 2, with the vector with inverted axis
        function vy = complexmultiplication_in(vx1, vx2)
            % z = a + j*b
            % z*w = a*c - b*d + j*(a*d + b*c)
            [m, L, d] = size(vx1);
            vy = zeros(m, L, 2);
            vy(:, 1:end, 1) = vx1(:, 1:end, 1).*vx2(:, 1:end, 1) - vx1(:, 1:end, 2).*vx2(:, 1:end, 2);
            vy(:, 1:end, 2) = vx1(:, 1:end, 2).*vx2(:, 1:end, 1) + vx1(:, 1:end, 1).*vx2(:, 1:end, 2);
        end

        %% Multiply two binary numbers
        function vy = binarymultiplication(vx1, vx2)
            [m1, L1, d1] = size(vx1);
            [m2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            vy = zeros(m, L1 + L2 - 1);
            vx1temp = vy;
            vx1temp(:, 1:L1-1) = vx1(:, 1:L1-1);
            vx2temp = vy;
            vx2temp(:, 1:L2-1) = vx2(:, 1:L2-1);

            for i = 1:L2-1
                vytemp = zeros(m, L1 + L2 - 1);
                vytemp(:, i:L1 - 1 + i) = and(vx1temp(:, 1:L1), vx2temp(:, i));
                vy_s = VHDL_Conversion.sum_signedx(vy, vytemp);
                vy(:, 1:L1 + L2 - 1) = vy_s(:, 1:L1 + L2 - 1);
            end
            vy(:, L1 + L2 - 1) = xor(vx1(:, L1), vx2(:, L2));
        end

        function vy = binarymultiplication_4d(vx1, vx2)
            %multiply two binary numbers
            [m1, n1, L1] = size(vx1);
            [m2, n2, L2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            n = max(n1,n2);
            vy = zeros(m, n, L1 + L2 - 1);
            vx1temp = vy;
            vx1temp(:, :, 1:L1-1) = vx1(:, :, 1:L1-1);
            vx2temp = vy;
            vx2temp(:, :, 1:L2-1) = vx2(:, :, 1:L2-1);

            for i = 1:L2-1
                vytemp = zeros(m, L1 + L2 - 1);
                vytemp(:, :, i:L1 - 1 + i) = and(vx1temp(:, :, 1:L1), vx2temp(:, :, i));
                vy_s = VHDL_Conversion.sum_unsigned(vy, vytemp);
                vy(:, :, 1:L1 + L2 - 1) = vy_s(:, :, 1:L1 + L2 - 1);
            end
            vy(:, :, L1 + L2 - 1) = xor(vx1(:, :, L1), vx2(:, :, L2));
        end

        %% Multiply two fixed point numbers
        function vy = FxP_multiplication(vx1, vx2)
            [m1, L1] = size(vx1);
            [m2, L2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            vy = zeros(m, L1 + L2 - 1);
            vx1temp = vy;
            vx1temp(:, 1:L1-1) = vx1(:, 1:L1-1);
            vx2temp = vy;
            vx2temp(:, 1:L2-1) = vx2(:, 1:L2-1);

            for i = 1:L2-1
                vytemp = zeros(m, L1 + L2 - 1);
                vytemp(:, i:L1 - 1 + i) = and(vx1temp(:, 1:L1), vx2temp(:, i));
                vy_s = VHDL_Conversion.sum_signedx(vy, vytemp);
                vy(:, 1:L1 + L2 - 1) = vy_s(:, 1:L1 + L2 - 1);
            end
            vy(:, L1 + L2 - 1) = xor(vx1(:, L1), vx2(:, L2));
        end

        %% Multiply two fixed point numbers with real part at row 1 and complex at row 2
        function vy = complex_binary_multiplication_FxP(vx1, vx2) 
            % z = a + j*b
            % z*w = a*c - b*d + j*(a*d + b*c)
            [m1, L1, d1] = size(vx1);
            [m2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            vy = zeros(m, L1 + L2 -1, 2);
            vy(1:end, :, 1) = VHDL_Conversion.minusx(VHDL_Conversion.FxP_multiplication(vx1(1:end, :, 1),vx2(1:end, :, 1)), VHDL_Conversion.FxP_multiplication(vx1(1:end, :, 2),vx2(1:end, :, 2)));
            vy(1:end, :, 2) = VHDL_Conversion.sum_signedx(VHDL_Conversion.FxP_multiplication(vx1(1:end, :, 2),vx2(1:end, :, 1)), VHDL_Conversion.FxP_multiplication(vx1(1:end, :, 1),vx2(1:end, :, 2)));
        end

        %% Multiply two binary numbers with real part at row 1 and complex at row 2
        function vy = complex_binary_multiplication(vx1, vx2)
            % z = a + j*b
            % z*w = a*c - b*d + j*(a*d + b*c)
            [m1, L1, d1] = size(vx1);
            [m2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            vy = zeros(m, L1 + L2 -1, 2);
            vy(1:end, :, 1) = VHDL_Conversion.minusx(VHDL_Conversion.binarymultiplication(vx1(1:end, :, 1),vx2(1:end, :, 1)), VHDL_Conversion.binarymultiplication(vx1(1:end, :, 2),vx2(1:end, :, 2)));
            vy(1:end, :, 2) = VHDL_Conversion.sum_signedx(VHDL_Conversion.binarymultiplication(vx1(1:end, :, 2),vx2(1:end, :, 1)), VHDL_Conversion.binarymultiplication(vx1(1:end, :, 1),vx2(1:end, :, 2)));
        end
    
        
        function vy = complex_binary_multiplication_4d(vx1, vx2)
            % z = a + j*b
            % z*w = a*c - b*d + j*(a*d + b*c)
            [m1, n1, L1, d1] = size(vx1);
            [m2, n2, L2, d2] = size(vx2);
            L = max(L1,L2);
            m = max(m1,m2);
            n = max(n1,n2);
            vy = zeros(m, n, L1 + L2, 2);
            for k =1:m
                vx1k (:, 1:L1, :) = vx1(k, :, 1:L1, :);
                vx2k (:, 1:L2, :) = vx2(k, :, 1:L2, :);
                vy(k, 1:end, :, 1) = VHDL_Conversion.minus(VHDL_Conversion.binarymultiplication(vx1k(1:end, :, 1),vx2k(1:end, :, 1)), VHDL_Conversion.binarymultiplication(vx1k(1:end, :, 2),vx2k(1:end, :, 2)));
                %a1 = VHDL_Conversion.binarymultiplication(vx1(1:end, :, 2),vx2(1:end, :, 1));
                %a2 = VHDL_Conversion.binarymultiplication(vx1(1:end, :, 1),vx2(1:end, :, 2));
                vy(k, 1:end, :, 2) = VHDL_Conversion.sum_signed(VHDL_Conversion.binarymultiplication(vx1k(1:end, :, 2),vx2k(1:end, :, 1)), VHDL_Conversion.binarymultiplication(vx1k(1:end, :, 1),vx2k(1:end, :, 2)));
            end
            vy = VHDL_Conversion.reduce_size_4d(vy, L);
        end

        %% Reduce the size of a binary number to a specific size cutting it
        %cuts symetrical to the middle
        function vy = reduce_size(vx, s)
            [m, L, d] = size(vx);
            L2 = L/2;
            s2= s/2;
            vy = zeros(m, s, d);
            vy(:, 1:s -1, :) = vx(:, ceil(L2-s2):floor(L2+s2) - 1, :);
            vy(:, s, :) = vx(:, L, :);
        end
        
        %% Reduce the size of a binary number to a specific size converting it
        function vy = reduce_size_x(vx, s)
            [m, L, d] = size(vx);
            x = VHDL_Conversion.FxP2Double(vx, floor(L/2), ceil(L/2));
            vy = VHDL_Conversion.Double2FxP(x, floor(s/2), ceil(s/2));
        end

        function vy = reduce_size_4d(vx, s)
            [m, n, L, d] = size(vx);
            L2 = ceil(L/2);
            s2= ceil(s/2);
            vy = zeros(m, n, s, d);
            vy(:, :, 1:s -1, :) = vx(:, :, L2-s2:L2+s2 - 2, :);
            vy(:, :, s, :) = vx(:, :, L, :);
        end

        %% Reduce the size of a binary vector with the axis inverted
        function vy = reduce_size_in(vx, s)
            [L, m, d] = size(vx);
            L2 = ceil(L/2);
            s2= ceil(s/2);
            vy = zeros(m, s, d);
            vy(1:s -1, :, :) = vx(L2-s2:L2+s2 - 2, :, :);
            vy(s, :, :) = vx(L, :, :);
        end




        %convert vectors
        function vy = bin2unsigned_vector(vx)
            %Convert a vector of binary to a vector of unsigned numbers
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i) = VHDL_Conversion.bin2unsigned(vx(i, :));
            end
        end

        function vy = bin2signed_signal_vector(vx)
            %Convert a vector of binary to signed
            %with the most significant bit beeing the signal
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i) = VHDL_Conversion.bin2signed_signal(vx(i, :));
            end
        end

        function vy = bin2signed_2s_vector(vx)
            %Convert a vector of binary to signed
            %with the binary as twos complement representation
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i) = VHDL_Conversion.bin2signed_2s(vx(i, :));
            end
        end

        function [vy, clipped] = unsigned2bin_vector(vx, L)
            %Convert a vector of unsigned numbers to a vector of binary
            %vx is the vector of numbers
            %L is the length of the binary vectors
            [m, n] = size(vx);
            vy = zeros(m, L);
            for i = 1:m
                [vy(i, :), clipped] = VHDL_Conversion.unsigned2bin(vx(i), L);
            end
        end

        function [vy, clipped] = signed2bin_signal_vector(vx, L)
            %Convert a vector of signed numbers to a vector of binary
            %with the most significant bit beeing the signal
            %vx is the vector of numbers
            %L is the length of the binary vectors
            [m, n] = size(vx);
            vy = zeros(m, L);
            for i = 1:m
                [vy(i, :), clipped] = VHDL_Conversion.signed2bin_signal(vx(i), L);
            end
        end

        function [vy, clipped] = signed2bin_2s_vector(vx, L)
            %Convert a vector of signed numbers to a vector of binary
            %with the binary as twos complement representation
            %vx is the vector of numbers
            %L is the length of the binary vectors
            [m, n] = size(vx);
            vy = zeros(m, L);
            for i = 1:m
                [vy(i, :), clipped] = VHDL_Conversion.signed2bin_2s(vx(i), L);
            end
        end

        function [vy, overflow] = Double2FxP_vector(vx, I, F)
            %Convert a vector of doubles to a vector of binary fixed point
            %with I beeing the length of the Integer part
            %with F beeing the length of the Farctional part
            %vx is the vector of numbers
            [m, n] = size(vx);
            vy = zeros(m, I+F);
            for i = 1:m
                [vy(i, :), overflow] = VHDL_Conversion.Double2FxP(vx(i), I, F);
            end
        end

        function vy = FxP2Double_vector(vx, I, F)
            % Convert a vector of binary fixed point to a vector of doubles
            %with I beeing the length of the Integer part
            %with F beeing the length of the Farctional part
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i) = VHDL_Conversion.FxP2Double(vx(i, :), I, F);
            end
        end

        function vy = lshift_vector(vx, s)
            % Perform left shift by s bits to every binary vector
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.lshift(vx(i, :), s);
            end
        end

        function vy = rshift_vector(vx, s)
            % Perform right shift by s bits to every binary vector
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.rshift(vx(i, :), s);
            end
        end

        function vy = invert_vector(vx)
            %invert the signt of all binary vectors
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.invert(vx(i, :));
            end
        end

        function vy = sum_vector(vx1, vx2)
            %execute a sum of two binary vectors in the same line
            %(vx1 and vx2) beeing the vector of binary numbers as vectors
            [m, n] = size(vx1);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.sum(vx1(i, :), vx2(i, :));
            end
        end

        function vy = minus_vector(vx1, vx2)
            %execute a subtraction of two binary vectors in the same line
            %(vx1 and vx2) beeing the vector of binary numbers as vectors
            [m, n] = size(vx1);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.minus(vx1(i, :), vx2(i, :));
            end
        end

        function vy = twoscomplement_vector(vx)
            %builds the twos complement of all binary vectors
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vy(i, :) = VHDL_Conversion.twoscomplement(vx(i, :));
            end
        end

        function vs = bits2str_vector(vx)
            %converts all binary vectors to binary string
            %vx is the vector of binary numbers as vectors
            [m, n] = size(vx);
            for i = 1:m
                vs(i, :) = VHDL_Conversion.bits2str(vx(i, :));
            end
        end

    end
end