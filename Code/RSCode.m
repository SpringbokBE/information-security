% File name:  RSCode.m
% Authors:    Samuel Rimbaut,
              Arne Vandenberghe,
              Johannes Van Wonterghem
% Date:       30/04/2020
% Emails:     samuel.rimbaut@ugent.be,
              arne.vandenberghe@ugent.be,
              johannes.vanwonterghem@ugent.be
% Brief:      E003600B, Information Theory, Project 
% About:      Implementation of Reed-Solomon coding.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef RSCode

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        m; % GF(2^m) field
        
        n; % Code length
        k; % Information length
        t; % Error correction capability
        l; % Shortened information length (-> shortened code length = l+n-k)
        
        m0; % m0 of the Reed-Solomon code, determines first root of generator
        
        g; % rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=RSCode(m,t,l,m0)
			% Constructor of the RSCode class
			% INPUT:
            % -m: Defines the GF(2^m) field
            % -t: Error correction capability
			% -l: Shortened information length (-> shortened code length = l+n-k)
			% -m0: m0 of the Reed-Solomon code, determines first root of generator
            % OUTPUT:
            % -obj: the RSCode object
            obj.m = m;
          
            obj.t = t;
            obj.l = l;
            
            obj.n = 2^m-1;
            obj.k = obj.n-2*t;
           
            obj.m0 = m0;
            
            obj.g = RSCode.makeGenerator(m,t,m0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function code = encode(obj,msg)
            % Systematically encode information words using the Reed-Solomon code
            % INPUT:
            % -obj: the RSCode object, defines the properties of the code
            % -mesig: every row contains a length obj.l information word consisting of GF(2^m) elements
            % OUTPUT:
            % -code: every row contains a GF(2^m) codeword corresponding to systematic Reed-Solomon coding of the corresponding information word            
            
            assert(size(msg,2) == obj.l);
            msg = [msg zeros(size(msg,1),obj.t*2)];
            enc = gf(zeros(size(msg)),obj.m);
            for i = 1:size(msg,1)
                [~,enc(i,1:end)] = deconv(msg(i,1:end),obj.g);
            end
            
            code = msg + enc;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [decoded,nERR] = decode(obj,code)
            % Decode Reed-Solomon codes
            % INPUT:
            % -obj: the RSCode object, defines the properties of the code
            % -code:  every row contains a GF(2^m) codeword of length obj.l+2*obj.t
            % OUTPUT:
            % -decoded: every row contains a GF(2^m) information word corresponding to decoding of the corresponding Reed-Solomon codeword
            % -nERR: column vector containing the number of corrected symbols for every codeword, -1 if error correction failed
            
            %frequency domain decoding
            msg = code;
            assert(size(msg,2) == obj.l+2*obj.t);
            msg = [zeros(size(msg,1), (obj.k-obj.l)) msg];%add zeros for complete FT
            decoded = gf(zeros(size(msg,1), obj.l),obj.m);
            nERR = zeros(size(msg,1), 1);
            %matrix of primitives for fourier
            alfa = gf(2,obj.m);
            
            for i = 1:size(msg,1)
                Sz = gf(zeros(1,2*obj.t),obj.m);
                %fourier transform
                for j = 1:(2*obj.t)%j goes from m0 to 2t+m0, only use it in calc
                    for fi = 1:size(msg,2)%use all the values of r for every Sj
                        Sz(j) = Sz(j) + (msg(i,(size(msg,2)+1-fi))*alfa^(j+obj.m0-1)^(fi-1));
                    end
                end
                Sz = RSCode.flip(Sz, obj.m);
                %Euclids algorithm
                %initialization
                delta0 = 0;
                delta1 = 1;
                delta2 = 1;
                omega0 = gf(zeros(1,(2*obj.t+1)), obj.m);
                omega0(1) = 1;%z^2t
                omega1 = Sz;
                [omega1, deg] = RSCode.Shorten(omega1, obj.m);
                while (deg >= obj.t)
                    [q,~] = deconv(omega0,omega1);
                    [omega0, qomega] = RSCode.SameSize(omega0, conv(q, omega1), obj.m);
                    [omega2, deg] = RSCode.Shorten(omega0 - qomega, obj.m);
                    [delta0, qdelta] = RSCode.SameSize(delta0, conv(q, delta1), obj.m);
                    [delta2, ~] = RSCode.Shorten(delta0 - qdelta, obj.m);
                    delta0 = delta1;%pass on for next rotation
                    delta1 = delta2;
                    omega0 = omega1;
                    omega1 = omega2;
                end
                delta2 = RSCode.flip(delta2, obj.m);
                Sz = RSCode.flip(Sz, obj.m);
                %recursive extension
                E = gf(zeros(1,size(msg,2)), obj.m);
                for j = 1:(2*obj.t)%assign already known values
                    E(j+obj.m0) = Sz(j);
                end
                %calculate other values with key equation
                %set delta(0) to 1
                if (delta2(1) ~= 0)
                    delta2 = delta2/delta2(1);
                end
                for j = 2*obj.t+1+obj.m0:size(msg,2)
                    for kk = 2:size(delta2,2)
                        E(j) = E(j) + delta2(kk)*E(j-kk+1);%+because of mod
                    end
                end
                if (delta2(size(delta2,2)) ~= 0)
                    delta2 = delta2/delta2(size(delta2,2));
                end
                for j = obj.m0:-1:1
                    for kk = size(delta2,2)-1:-1:1
                        E(j) = E(j) + delta2(kk)*E(j+size(delta2,2)-kk);%+because of mod
                    end
                end
                
                %reverse fourier
                e = gf(zeros(1,size(msg,2)), obj.m);
                nerrors = 0;
                for fi = 1:size(msg,2)
                    for fj = 1:size(msg,2)
                        e(fi) = e(fi) + E(fj)*alfa^(-(fi-1)*(fj-1));
                    end
                    if (e(fi) ~= 0)
                        nerrors = nerrors+1;
                    end
                end
                e = RSCode.flip(e, obj.m);
                %receive word
                nERR(i) = (nerrors <=obj.t)*nerrors+(nerrors>obj.t)*-1; %too many errors cannot be corrected
                c = msg(i, 1:end) - e;
                decoded(i,1:end) = c((obj.k-obj.l+1):obj.k);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function generator = makeGenerator(m, t, m0)
            % Generate the Reed-Solomon generator polynomial with error correcting capability t over GF(2^m)
            % INPUT:
            % -m: positive integer
            % -t: error correction capability of the Reed-Solomon code, positive integer > 1
            % -m0: determines first root of generator polynomial, positve integer >= 0
            % OUTPUT:
            % -generator: rowvector containing the GF(2^m) generator polynomial coefficients in order of descending powers
            alfa = gf(2,m);
            generator = gf([1,alfa^m0],m);
            for i = 1:2*t-1
                generator = conv(generator,gf([1,alfa^(m0+i)],m));
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [vector, degree] = Shorten(g, m)
            % conv wont work with zeros in array
            degree = size(g,2);
            for i = 1:size(g,2)
                if g(i) == 0
                    degree = degree - 1;
                else
                    break;
                end
            end
            vector = gf(zeros(1,(degree)), m);
            offset = size(g,2) - (degree);
            for i = 1:(degree)
                vector(i) = g(offset + i);
            end
            degree = degree -1; %5 elements -> z^4 -> degree =4
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function f = flip(g, m)
            % in matlab polynomials are in reverse order than gf
            f = gf(zeros(1, size(g,2)), m);
            for i = 1:size(g,2)
                f(i) = g(size(g,2)+1-i);
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [f2, g2] = SameSize(f, g, m)
            % + and - wont work when not the same size....
            offset = size(f,2) - size(g,2);
            if offset > 0
                f2 = f;
                g2 = gf(zeros(1, size(f,2)), m);
                for i = 1:(size(f,2))
                    if (i <= offset)
                        g2(i) = 0;
                    else
                        g2(i) = g(i-offset);
                    end
                end
            elseif offset < 0
                g2 = g;
                f2 = gf(zeros(1, size(g,2)), m);
                for i = 1:(size(g,2))
                    if (i <= -offset)
                        f2(i) = 0;
                    else
                        f2(i) = f(i+offset);
                    end
                end
            else
                g2 = g;
                f2 = f;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function test()
            % Test the Matlab code of this class
            
            m0 = 0; % Also test with other values of m0!
            %alfa = gf(2,4);
            %rs = RSCode(4,3,6,2);
            %msg = gf([alfa,alfa^5,0,alfa,alfa^4,1],4);
            %code = rs.encode(msg);
            %code(1, [7,12]) = [alfa^7,alfa]
            %[decoded,nERR] = rs.decode(code)
            %assert(1==0)
            rs = RSCode(8,5,10,m0); % Construct the RSCode object
            
            msg = gf(randi([0,2^8-1],5,10),8) % Generate a random message of 5 information words
            code = rs.encode(msg); % Encode this message
            % Introduce errors
            code(2,[3 18]) = code(2,[5 18])+1;
            code(3,8) = 0;
            code(4,[4 2 19 20 6]) = gf(randi([0,2^8-1],1,5),8);
            code(5,[4 2 19 20 6 13]) = gf(randi([0,2^8-1],1,6),8);
            
            
            [decoded,nERR] = rs.decode(code); % Decode
            
            nERR
            decoded
            assert(all(all(decoded(1:4,:) == msg(1:4,:))))
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%