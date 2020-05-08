% File name:  AudioCD.m
% Authors:    Gerbrand De Laender,
%             Jan Mariën,
%             Samuel Rimbaut,
%             Arne Vandenberghe,
%             Johannes Van Wonterghem
% Date:       30/04/2020
% Emails:     gerbrand.delaender@ugent.be,
%             jan.marien@ugent.be,
%             samuel.rimbaut@ugent.be,
%             arne.vandenberghe@ugent.be
%             johannes.vanwonterghem@ugent.be
% Brief:      E003600B, Information Theory, Project 
% About:      Implementation of Cross-Interleaved Reed-Solomon coding on an
              Audio CD.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef AudioCD_opt

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties
        Fs; % Sample rate of the audio
        
        configuration; % % 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
        max_interpolation; % The maximum number of interpolated audio samples
        
        enc2; % Outer RS code of the CIRC
        dec2;
        gpoly2;
        
        enc1; % Inner RS code of the CIRC
        dec1;
        gpoly1;
        
        gpoly_8_parity; % Single RS code for configuration 3
        enc_8_parity;
        dec_8_parity;
                
        cd_bits; % Bits written to disk (before EFM)
        
        scaled_quantized_padded_original; % Reference to compare the output of readCD to
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=AudioCD_opt(Fs,configuration,max_interpolation)  
            % Constructor of the AudioCD_opt class
			% INPUT:
            % -Fs: The sample rate of the audio
            % -configuration: 0: no CIRC; 1: CIRC as described in standard; 2: Concatenated RS, no interleaving; 3: Single 32,24 RS
			% -max_interpolation: The maximum number of interpolated audio samples
            % OUTPUT:
            % -obj: the AudioCD_opt object

            obj.Fs = Fs;
            obj.max_interpolation = max_interpolation;
            
            % Initialize the RS encoders and decoders
            primpoly = [1,0,0,0,1,1,1,0,1]; % Primitive polynomial from the standard
            if (configuration == 1) || (configuration == 2)
                obj.gpoly2 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc2 = comm.RSEncoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec2 = comm.RSDecoder(255,251,obj.gpoly2,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);

                obj.gpoly1 = rsgenpoly(255,251,bi2de(fliplr(primpoly)),0);
                obj.enc1 = comm.RSEncoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec1 = comm.RSDecoder(255,251,obj.gpoly1,28,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly,'ErasuresInputPort',true);
            elseif configuration == 3
                obj.gpoly_8_parity = rsgenpoly(255,247,bi2de(fliplr(primpoly)),0);
                obj.enc_8_parity = comm.RSEncoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);
                obj.dec_8_parity = comm.RSDecoder(255,247,obj.gpoly_8_parity,24,'PrimitivePolynomialSource','Property','PrimitivePolynomial',primpoly);           
            end
            obj.configuration = configuration;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=writeCd(obj,audiofile)
            % Write an audiofile to the CD
            % INPUT:
            % -obj: the AudioCD_opt object
            % -audiofile: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype
            % OUTPUT:
            % -obj: the updated AudioCD_opt object            

            assert(size(audiofile,2) == 2);
            
            xscaled = audiofile / max(max(abs(audiofile))); % normalize to -1:1
            x = uencode(xscaled,16); % convert to 16 bit signed values

            xlr16 = reshape(x',[],1); % serialize left and right audio channel
            xlr8 = typecast(xlr16,'uint8'); % split into 8 bit words

            xlr8_padded = [xlr8 ; zeros(24-(rem(numel(xlr8)-1,24)+1),1)]; % pad with zeros to fill an integer number of frames
            n_frames = numel(xlr8_padded)/24; % every frame contains 24 8 bit words
            
            ylr16 = typecast(uint8(xlr8_padded),'uint16');
            y = reshape(ylr16,2,[])';
            obj.scaled_quantized_padded_original = udecode(y,16); % Reference to compare the output of readCD to
            
            switch obj.configuration
                case 0 % no CIRC
                    xlrb = de2bi(xlr8,8);
                case 1 % CIRC as described in standard
                    [delay_interleaved,n_frames] = obj.CIRC_enc_delay_interleave_opt(xlr8_padded,n_frames);
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2_opt(delay_interleaved,n_frames);
                    [delay_unequal,n_frames] = obj.CIRC_enc_delay_unequal_opt(C2_encoded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1_opt(delay_unequal,n_frames);
                    [delay_inv,n_frames] = obj.CIRC_enc_delay_inv_opt(C1_encoded,n_frames);

                    xlrb = de2bi(delay_inv,8);
                case 2 % Concatenated RS, no interleaving
                    [C2_encoded,n_frames] = obj.CIRC_enc_C2(xlr8_padded,n_frames);
                    [C1_encoded,n_frames] = obj.CIRC_enc_C1(C2_encoded,n_frames);
                    
                    xlrb = de2bi(C1_encoded,8);
                case 3 % Single 32,24 RS
                    [encoded,n_frames] = obj.C3_enc_8_parity(xlr8_padded,n_frames);
                    
                    xlrb = de2bi(encoded,8);
                otherwise
                    error('Invalid configuration selected');
            end
            
            xlrbserial = reshape(xlrb',[],1);
            
            obj.cd_bits = logical(xlrbserial);            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj=bitErrorsCd(obj,p)
            % Add uniform bit errors to cd
            % INPUT:
            % -obj: the AudioCD_opt object
            % -p: the bit error probability, i.e., an obj.cd_bits bit is flipped with probability p
            % OUTPUT:
            % -obj: the updated AudioCD_opt object            

            noise = rand(size(obj.cd_bits))<p;
            
            obj.cd_bits = xor(obj.cd_bits,noise);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function obj=scratchCd(obj,length_scratch,location_scratch)
            % Add a scratch to the cd
            % INPUT:
            % -obj: the AudioCD_opt object
            % -length_scratch: the length of the scratch (in number of bit)
            % -location_scratch: the location of the scratch (in bit offset from start of obj.cd_bits)
            % OUTPUT:
            % -obj: the updated AudioCD_opt object            
            
            obj.cd_bits(location_scratch:min(location_scratch+length_scratch-1,numel(obj.cd_bits))) = 0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [audio_out,interpolation_flags]=readCd(obj)
            % Read an audiofile from the CD
            % INPUT:
            % -obj: the AudioCD_opt object, containing the data in obj.cd_bits
            % OUTPUT:
            % -audio_out: [Nsamples x 2] matrix containing the left and right audio track as samples of the double datatype  
            % -interpolation_flags: [Nsamples x 2] matrix containing a 0 where no erasure was flagged, a 1 where an erasure was interpolated and a -1
            % where interpolation failed

            
            ylrb = reshape(obj.cd_bits,8,[])';
            ylr8 = bi2de(ylrb);
            
            switch obj.configuration
                case 0 % no CIRC
                    ylr16 = typecast(uint8(ylr8),'uint16');
                    y = reshape(ylr16,2,[])';
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                case 1
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [delay_inv,n_frames] = obj.CIRC_dec_delay_inv_opt(ylr8,n_frames);
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1_opt(delay_inv,n_frames);
                    [delay_unequal,erasure_flags,n_frames] = obj.CIRC_dec_delay_unequal_opt(C1_decoded,erasure_flags,n_frames);
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2_opt(delay_unequal,erasure_flags,n_frames);
                    [deinterleave_delay,erasure_flags,n_frames] = obj.CIRC_dec_deinterleave_delay_opt(C2_decoded,erasure_flags,n_frames);

                    ylr16 = typecast(uint8(deinterleave_delay),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                    
                case 2 % Concatenated RS, no interleaving
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [C1_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C1(ylr8,n_frames);
                    erasure_flags_t = erasure_flags;
                    [C2_decoded,erasure_flags,n_frames] = obj.CIRC_dec_C2(C1_decoded,erasure_flags,n_frames);
                    
                    if(numel(erasure_flags)  ~= numel(C2_decoded))
                        disp('Something wrong!');
                    end
                    
                    ylr16 = typecast(uint8(C2_decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                case 3 % Single 32,24 RS
                    n_frames = numel(ylr8)/32;
                    assert(n_frames*32 == numel(ylr8));
                    
                    [decoded,erasure_flags,n_frames] = obj.C3_dec_8_parity(ylr8,n_frames);
                    
                    ylr16 = typecast(uint8(decoded),'uint16');
                    y = reshape(ylr16,2,[])';
                    
                    erasure_flags = reshape(erasure_flags,2,[]);
                    erasure_flags = or(erasure_flags(1,:),erasure_flags(2,:))';
                    erasure_flags = reshape(erasure_flags,2,[])';
                    
                    % Linear Interpolation
                    interpolation_failed = zeros(size(erasure_flags),'logical');
                    [y(:,1),interpolation_failed(:,1)] = obj.interpolator(y(:,1),erasure_flags(:,1)); % Left
                    [y(:,2),interpolation_failed(:,2)] = obj.interpolator(y(:,2),erasure_flags(:,2)); % Right
                    
                    audio_out = udecode(y,16);
                    interpolation_flags = zeros(size(audio_out));
                    interpolation_flags(erasure_flags) = 1;
                    interpolation_flags(interpolation_failed) = -1;
                otherwise
                    error('Invalid configuration selected');
            end
            
        end
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_enc_delay_interleave(obj,input,n_frames)
            % CIRC Encoder: Delay of 2 frames + interleaving sequence
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            buf1 = zeros(24,1,'uint8');
            buf2 = zeros(24,1,'uint8');
            output = zeros((n_frames+2)*24,1,'uint8');
            
            for i = 1:n_frames
                frame = input((i-1)*24+1:i*24);
                output((i-1)*24+1:i*24) = [buf1(1:2); buf1(9:10); buf1(17:18); buf1(3:4); buf1(11:12); ...
                    buf1(19:20); frame(5:6); frame(13:14); frame(21:22); frame(7:8); frame(15:16); frame(23:24)];
                buf1 = buf2;
                buf2 = frame;
            end
            for i = n_frames+1:n_frames+2
                frame = zeros(24,1,'uint8');
                                output((i-1)*24+1:i*24) = [buf1(1:2); buf1(9:10); buf1(17:18); buf1(3:4); buf1(11:12); ...
                buf1(19:20); frame(5:6); frame(13:14); frame(21:22); frame(7:8); frame(15:16); frame(23:24)];
                buf1 = buf2;
                buf2 = frame;
            end
            n_frames = n_frames+2;        
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~70x faster
        function [output,n_frames] = CIRC_enc_delay_interleave_opt(obj,input,n_frames)
            % CIRC Encoder: Delay of 2 frames + interleaving sequence
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            input = reshape(input, 24, n_frames);
            output = zeros(24, n_frames + 2,'uint8');

            output([1:2 7:8 3:4 9:10 5:6 11:12],3:end) = input([1:4 9:12 17:20],:);
            output([13:14 19:20 15:16 21:22 17:18 23:24],1:end-2) = input([5:8 13:16 21:24],:);

            output = reshape(output,[],1);
            n_frames = n_frames + 2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_enc_C2(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C2)
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*28, 1, 'uint8');
            for i = 1:n_frames
                frame = input((i-1)*24+1:i*24);
                temp_frame = step(obj.enc2, frame);
                output((i-1)*28+1:i*28) = [temp_frame(1:12); temp_frame(25:28); temp_frame(13:24)]; % Put parity bits in the middle
            end 

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~1.05x faster
        function [output,n_frames] = CIRC_enc_C2_opt(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C2)
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames         
            %             input = reshape(input, n_frames, 24);
            input = reshape(input, 24, n_frames);
            output = zeros(28, n_frames, 'uint8');
            
            for i = 1:n_frames
                output(:,i) = obj.enc2(input(:,i));
            end
            
            output = output([1:12 25:28 13:24],:);
            output = reshape(output,[],1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_enc_delay_unequal(obj,input,n_frames)
            % CIRC Encoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
           
           D = 4; 
           buffer = zeros(27*D + 1, 28, 'uint8');
           output = zeros((n_frames + 108) * 28, 1, 'uint8');
           for i = 1:n_frames+108
                buffer = circshift(buffer, [1,0]); % Shift rows
                if (i > n_frames)
                    frame = zeros(28, 1, 'uint8');
                else
                    frame = input((i-1)*28+1:i*28);
                end
                buffer(1,:) = frame;
                new_frame = zeros(28, 1, 'uint8');
                for k = 0:27
                    new_frame(k+1) = buffer((k*D + 1), k+1);
                end
                output((i-1)*28+1:i*28) = new_frame;
           end
           n_frames = n_frames + 108;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~40x faster
        function [output,n_frames] = CIRC_enc_delay_unequal_opt(obj,input,n_frames)
            % CIRC Encoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            input = reshape(input, 28, n_frames);
            output = zeros(28, n_frames + 108,'uint8');
            
            for d = 0:27
                output(d+1,4*d+1:end-4*(27-d)) = input(d+1,:);
            end

            output = reshape(output,[],1);
            n_frames = n_frames + 108;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_enc_C1(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C1)
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames
            output = zeros(n_frames*32, 1, 'uint8');
            for i = 1:n_frames
                frame = input((i-1)*28+1:i*28);
                temp_frame = step(obj.enc1, frame);
                output((i-1)*32+1:i*32) = temp_frame;
            end 

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~1.05x faster
        function [output,n_frames] = CIRC_enc_C1_opt(obj,input,n_frames)
            % CIRC Encoder: Generation of 4 parity symbols (C1)
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames
            input = reshape(input, 28, n_frames);
            output = zeros(32, n_frames,'uint8');
            
            for i = 1:n_frames
                output(:,i) = obj.enc1(input(:,i));
            end
            
            output = reshape(output,[],1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_enc_delay_inv(obj,input,n_frames)
            % CIRC Encoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            output = zeros(32*(n_frames+1),1,'uint8');
            for i = 1:n_frames
                for(j = 1:32)
                    if(mod(j,2) == 0)
                        output((i-1)*32+j) = input((i-1)*32+j);
                        if ((j>12 && j <17) || (j>28))
                            output((i-1)*32+j) = bitcmp(output((i-1)*32+j));
                        end
                    else
                        output(i*32+j) = input((i-1)*32+j);
                        if ((j>12 && j <17) || (j>28))
                            output(i*32+j) = bitcmp(output(i*32+j));
                        end
                    end
                end
            end
            n_frames = n_frames+1;
            

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~5x faster
        function [output,n_frames] = CIRC_enc_delay_inv_opt(obj,input,n_frames)
            % CIRC Encoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC encoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC encoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            input = reshape(input, 32, n_frames);
            output = zeros(32, n_frames + 1,'uint8');
            
            output(1:2:end,2:end) = input(1:2:end,:);
            output(2:2:end,1:end-1) = input(2:2:end,:);
            output([13 15 29 31],:) = bitcmp(output([13 15 29 31],:),'uint8');
            output([14 16 30 32],:) = bitcmp(output([14 16 30 32],:),'uint8');
            
            output = reshape(output,[],1);
            n_frames = n_frames + 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = CIRC_dec_delay_inv(obj,input,n_frames)
            % CIRC Decoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_inv(obj.CIRC_enc_delay_inv(x)) == x!
            output = zeros(32*(n_frames),1,'uint8');
            for i = 1:n_frames
                for j = 1:32
                    if(mod(j,2) == 1)
                        output((i-1)*32+j) = input((i-1)*32+j);
                        if ((j>12 && j <17) || (j>28))
                            output((i-1)*32+j) = bitcmp(output((i-1)*32+j));
                        end
                    else
                        output(i*32+j) = input((i-1)*32+j);
                        if ((j>12 && j <17) || (j>28))
                            output(i*32+j) = bitcmp(output(i*32+j));
                        end
                    end
                    
                        
                end
                
            end
            output = output(33:end-32); %remove first frame
            n_frames = n_frames-1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~2.5x faster
        function [output,n_frames] = CIRC_dec_delay_inv_opt(obj,input,n_frames)
            % CIRC Decoder: Delay of 1 frame + inversions
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -n_frames: the length of the output expressed in frames (changed from input because of delay!)
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_inv(obj.CIRC_enc_delay_inv(x)) == x!
            input = reshape(input, 32, n_frames);
            output = zeros(32, n_frames - 1,'uint8');

            output(1:2:end,:) = input(1:2:end,2:end);
            output(2:2:end,:) = input(2:2:end,1:end-1);
            output([13 15 29 31],:) = bitcmp(output([13 15 29 31],:),'uint8');
            output([14 16 30 32],:) = bitcmp(output([14 16 30 32],:),'uint8');

            output = reshape(output,[],1);
            n_frames = n_frames - 1;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,erasure_flags_out,n_frames] = CIRC_dec_C1(obj,input,n_frames)
            % CIRC Decoder: C1 decoder
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            output = zeros(n_frames*28,1,'uint8');
            erasure_flags_out = zeros(n_frames*28,1,'logical');
            for i = 1:n_frames
                frame = input((i-1)*32+1:i*32);
                erasure_in = zeros(32,1,'logical');
                [output_dec,ERR] = step(obj.dec1,frame,erasure_in);
                if (ERR == 0 || ERR == 1) 
                    output((i-1)*28+1:i*28) = output_dec;
                else
                    output((i-1)*28+1:i*28) = output_dec;%not really needed
                    erasure_flags_out((i-1)*28+1:i*28) = 1;
                end
                
            end

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~same execution time
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C1_opt(obj,input,n_frames)
            % CIRC Decoder: C1 decoder
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            input = reshape(input, 32, n_frames);
            output = zeros(28, n_frames, 'uint8');
            erasure_flags_in = zeros(32, 1, 'logical');
            erasure_flags_out = zeros(28, n_frames, 'logical');
            
            for i = 1:n_frames
                [output(:,i), err] = obj.dec1(input(:,i),erasure_flags_in);
                if (err ~= 0 && err ~= 1)
                    erasure_flags_out(:,i) = 1;
                end
            end
                
            output = reshape(output,[],1);
            erasure_flags_out = reshape(erasure_flags_out,[],1);
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,erasure_flags_out,n_frames] = CIRC_dec_delay_unequal(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (eraure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_unequal(obj.CIRC_enc_delay_unequal(x)) == x!
           
           D = 4; 
           buffer = zeros(27*D + 1, 28, 'uint8');
           erasure_buffer = zeros(27*D + 1, 28, 'logical');
           output = zeros((n_frames) * 28, 1, 'uint8');
           erasure_flags_out = zeros(n_frames*28,1,'logical');
           for i = 1:n_frames
                frame = input((i-1)*28+1:i*28);
                err = erasure_flags_in((i-1)*28+1:i*28);
                buffer = circshift(buffer, [1,0]); % Shift rows
                erasure_buffer = circshift(erasure_buffer, [1,0]);
                buffer(1,:) = frame;
                erasure_buffer(1, :) = err;
                
                new_frame = zeros(28, 1, 'uint8');
                new_err = zeros(28, 1, 'logical');
                for k = 0:27
                    new_frame(k+1) = buffer((27-k)*D + 1, k+1);
                    new_err(k+1) = erasure_buffer((27-k)*D + 1, k+1);
                end
                output((i-1)*28+1:i*28) = new_frame;
                erasure_flags_out((i-1)*28+1:i*28) = new_err;
           end
           n_frames = n_frames - 108;
           output = output(108 * 28 + 1 : end);
           erasure_flags_out=  erasure_flags_out(108 * 28 + 1 : end);

        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

         % ~45x faster
        function [output,erasure_flags_out,n_frames] = CIRC_dec_delay_unequal_opt(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: Delay lines of unequal length
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (eraure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_delay_unequal(obj.CIRC_enc_delay_unequal(x)) == x!
            input = reshape(input, 28, n_frames);
            erasure_flags_in = reshape(erasure_flags_in, 28, n_frames);
            output = zeros(28, n_frames - 108,'uint8');
            erasure_flags_out = zeros(28, n_frames - 108,'logical');
                        
            for d = 0:27
                output(d+1,:) = input(d+1,4*d+1:end-4*(27-d));
                erasure_flags_out(d+1,:) = erasure_flags_in(d+1,4*d+1:end-4*(27-d));
            end

            output = reshape(output,[],1);
            erasure_flags_out = reshape(erasure_flags_out,[],1);
            n_frames = n_frames - 108;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,erasure_flags_out,n_frames] = CIRC_dec_C2(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: C2 decoder
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder            
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            for i = 1:n_frames    
                %step(decoder,input,input_erasures_from_c1_deinterleaved)
                %put parity bits at the end first
                flags_temp = [erasure_flags_in((i-1)*28+1:(i-1)*28+12); ...
                    erasure_flags_in((i-1)*28+17:i*28); erasure_flags_in((i-1)*28+13:(i-1)*28+16)];
                f = sum(flags_temp);
                encode_temp = [input((i-1)*28+1:(i-1)*28+12); ...
                    input((i-1)*28+17:i*28); input((i-1)*28+13:(i-1)*28+16)];
                if(f>=5) %step throws an error if f>4
                    ERR = -1;
                    output_dec = encode_temp(1:24);
                    erasure_flags_out((i-1)*24+1:i*24) = flags_temp(1:24);
                else
                    [output_dec,ERR] = step(obj.dec2,encode_temp,flags_temp);
                end
                
                if (ERR == 0) || (ERR == 1) %0 or 1 errors -> decode
                    output((i-1)*24+1:i*24) = output_dec;
                else
                    if f>2 %copy erasure flags and input to the output
                        output((i-1)*24+1:i*24) = encode_temp(1:24);
                        erasure_flags_out((i-1)*24+1:i*24) = flags_temp(1:24);
                    elseif (f==2) && (ERR ~= -1) %2 erasures and decoding succesful
                        output((i-1)*24+1:i*24) = output_dec;
                    else
                        output((i-1)*24+1:i*24) = output_dec(1:24); 
                        erasure_flags_out((i-1)*24+1:i*24) = 1;
                    end
                end
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~1.3x faster
        function [output,erasure_flags_out,n_frames] = CIRC_dec_C2_opt(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: C2 decoder
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder            
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block, follow the decoding algorithm from the assignment
            % -n_frames: the length of the output expressed in frames
            input = reshape(input, 28, n_frames);
            input = input([1:12 17:28 13:16],:); % Put parity bits at the end.
            erasure_flags_in = reshape(erasure_flags_in, 28, n_frames);
            erasure_flags_in = erasure_flags_in([1:12 17:28 13:16],:); % Put parity bits at the end.
            output = zeros(24, n_frames, 'uint8');
            erasure_flags_out = zeros(24, n_frames, 'logical');
            
            f = sum(erasure_flags_in);
            
            for i = 1:n_frames
                if (f(i) > 4)
                    output(:,i) = input(1:end-4,i);
                    erasure_flags_out(:,i) = erasure_flags_in(1:end-4,i);
                else
                    [output(:,i), err] = obj.dec2(input(:,i),erasure_flags_in(:,i));
                    if (err ~= 0 && err ~= 1)
                        if (f(i) < 2)
                            erasure_flags_out(:,i) = 1;
                        elseif (f(i) > 2)
                            output(:,i) = input(1:end-4,i);
                            erasure_flags_out(:,i) = erasure_flags_in(1:end-4,i);
                        end
                    end
                end    
            end
            
            output = reshape(output,[],1);
            erasure_flags_out = reshape(erasure_flags_out,[],1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,erasure_flags_out,n_frames] = CIRC_dec_deinterleave_delay(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: De-interleaving sequence + delay of 2 frames
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (erasure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_deinterleave_delay(obj.CIRC_enc_delay_interleave(x)) == x!
            buf1 = zeros(24,1,'uint8');
            buf2 = zeros(24,1,'uint8');
            fbuf1 = zeros(24,1,'logical');
            fbuf2 = zeros(24,1,'logical');
            output = zeros((n_frames-2)*24,1,'uint8');
            for i = 1:n_frames
                frame = input((i-1)*24+1:i*24);
                er_frame = erasure_flags_in((i-1)*24+1:i*24);
                output((i-1)*24+1:i*24) = [frame(1:2); frame(7:8); buf2(13:14); buf2(19:20); frame(3:4); ...
                    frame(9:10); buf2(15:16); buf2(21:22); frame(5:6); frame(11:12); buf2(17:18); buf2(23:24)];
                erasure_flags_out((i-1)*24+1:i*24) = [er_frame(1:2); er_frame(7:8); fbuf2(13:14); ...
                    fbuf2(19:20); er_frame(3:4); er_frame(9:10); fbuf2(15:16); fbuf2(21:22); ...
                    er_frame(5:6); er_frame(11:12); fbuf2(17:18); fbuf2(23:24)];
                buf2 = buf1;
                fbuf2 = fbuf1;
                buf1 = frame;
                fbuf1 = er_frame;
            end
            %first 2 frames should be empty because everything is delayed
            %by 2 frames
            output = output(2*24+1:end);
            n_frames = n_frames - 2;
            erasure_flags_out = erasure_flags_out(2*24+1:end);
            
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % ~120x faster
        function [output,erasure_flags_out,n_frames] = CIRC_dec_deinterleave_delay_opt(obj,input,erasure_flags_in,n_frames)
            % CIRC Decoder: De-interleaving sequence + delay of 2 frames
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block of the CIRC decoder
            % -erasure_flags_in: the erasure flags at the input of this block of the CIRC decoder
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block of the CIRC decoder  
            % -erasure_flags_out: the erasure flags at the output of this block (erasure flags follow same interleaving as the data)
            % -n_frames: the length of the output expressed in frames
            
            % Note: remove empty frames such that obj.CIRC_dec_deinterleave_delay(obj.CIRC_enc_delay_interleave(x)) == x!
            input = reshape(input, 24, n_frames);
            erasure_flags_in = reshape(erasure_flags_in, 24, n_frames);
            output = zeros(24, n_frames - 2,'uint8');
            erasure_flags_out = zeros(24, n_frames - 2,'logical');
            
            output([1:2 9:10 17:18 3:4 11:12 19:20],:) = input(1:12,3:end);
            erasure_flags_out([1:2 9:10 17:18 3:4 11:12 19:20],:) = erasure_flags_in(1:12,3:end);
            output([5:6 13:14 21:22 7:8 15:16 23:24],:) = input(13:24,1:end-2);
            erasure_flags_out([5:6 13:14 21:22 7:8 15:16 23:24],:) = erasure_flags_in(13:24,1:end-2);
            
            output = reshape(output,[],1);
            erasure_flags_out = reshape(erasure_flags_out,1,[]); % TODO : Why isn't this a column vector, same dimension as output?
            
            n_frames = n_frames - 2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,n_frames] = C3_enc_8_parity(obj,input,n_frames)
            % Configuration 3: Generation of 8 parity symbols
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -n_frames: the length of the output expressed in frames
            
            output = zeros(n_frames*32,1,'uint8');
            for i = 1:n_frames
                output((i-1)*32+1:i*32) = step(obj.enc_8_parity,input((i-1)*24+1:i*24));
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,erasure_flags_out,n_frames] = C3_dec_8_parity(obj,input,n_frames)
            % Configuration 3: Decoder
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block
            % -n_frames: the length of the input expressed in frames
            % OUTPUT:
            % -output: the output of this block
            % -erasure_flags_out: the erasure flags at the output of this block
            % -n_frames: the length of the output expressed in frames

            output = zeros(n_frames*24,1,'uint8');
            erasure_flags_out = zeros(n_frames*24,1,'logical');
            for i = 1:n_frames
                [output_dec,ERR] = step(obj.dec_8_parity,input((i-1)*32+1:i*32));
                if ERR == -1
                    output((i-1)*24+1:i*24) = output_dec;
                    erasure_flags_out((i-1)*24+1:i*24) = 1;
                else
                    output((i-1)*24+1:i*24) = output_dec;
                end
                
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [output,interpolation_failed] = interpolator(obj,input,erasure_flags_in)
            % Interpolation: Linear interpolation
            % INPUT:
            % -obj: the AudioCD_opt object
            % -input: the input to this block
            % -erasure_flags_in: the erasure flags at the input of this block
            % OUTPUT:
            % -output: linear interpolation of the input where there are no more than obj.max_interpolation consecutive erasures
            % -interpolation_failed: equal to one at the samples where interpolation failed

			% Make sure first and last sample are defined. If erasure: set
			% to 0 (corresponds to the unsigned integer value 2^15).
            if erasure_flags_in(1) ~= 0
                erasure_flags_in(1) = 0;
                input(1) = int16(2^15);
            end
            if erasure_flags_in(end) ~= 0
                erasure_flags_in(end) = 0;
                input(end) = int16(2^15);
            end
            
            % Find the number of consecutive erasures.
            n_consecutive_erasures = zeros(size(erasure_flags_in));
            x = [0,erasure_flags_in',0];
            y = strfind(x, [0 1]);   
            n_consecutive_erasures(y) = strfind(x, [1 0]) - y;
            
            % Apply linear interpolation as long as the number of 
            % consecutive erasures is smaller than the max.
            output = input;
            interpolation_failed = erasure_flags_in;
            for i = find( (n_consecutive_erasures>0) & (n_consecutive_erasures<= obj.max_interpolation))'
                output(i:i+n_consecutive_erasures(i)-1) = round(double(output(i-1))+(1:n_consecutive_erasures(i))*(double(output(i+n_consecutive_erasures(i)))-double(output(i-1)))/(n_consecutive_erasures(i)+1));
                interpolation_failed(i:i+n_consecutive_erasures(i)-1) = 0;
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods(Static)
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function plot(configuration)
            % Test the Matlab code of this class
            
            % audio_file = audioread('Hallelujah.wav');
            % Fs = 44.1e3;
            
            audio_file = audioread('Hallelujah_22050.wav');
            Fs = 22.05e3;
            x = logspace(-1-log10(2),-3,10);
            cd = AudioCD_opt(Fs,configuration,8);
            prob_flag = zeros(size(x));
            prob_failed = zeros(size(x));
            
            cd = cd.writeCd(audio_file);
            ogbits = cd.cd_bits;
            for i=1:size(x,2)
                cd.cd_bits = ogbits;
                cd = cd.bitErrorsCd(x(i));
                [out,interpolation_flags] = cd.readCd();
                sum1 = sum(sum(interpolation_flags~=0));
                sum2 = sum(sum(interpolation_flags==-1));
                sizeflags = size(interpolation_flags,1)*size(interpolation_flags,2);
                prob_flag(i) = sum1 /sizeflags;
                prob_failed(i) = sum2 /sizeflags;      
            end
            semilogx((x), (prob_flag));
            hold on
            semilogx(x, prob_failed);
            hold off
            xlabel('bit error probability')
            ylabel('probability')
            legend('probability of erasure flag', 'probability of failed interpolation');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function test()
            % Test the Matlab code of this class
            
            % audio_file = audioread('Hallelujah.wav');
            % Fs = 44.1e3;
            
            audio_file = audioread('Hallelujah_22050.wav');
            Fs = 22.05e3;


            cd = AudioCD_opt(Fs,1,8);
            tic
            cd = cd.writeCd(audio_file);
            toc
            
            T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
            for i = 1:floor(numel(cd.cd_bits)/T_scratch)
               cd = cd.scratchCd(4000,1+(i-1)*T_scratch);
            end
            % cd = cd.scratchCd(5000,1);
            % cd = cd.bitErrorsCd(0.0001);
            tic
            [out,interpolation_flags] = cd.readCd();
            toc

            fprintf('Number samples with erasure flags: %d\n',sum(sum(interpolation_flags~=0)));
            fprintf('Number samples with failed interpolations: %d\n',sum(sum(interpolation_flags==-1)));
            fprintf('Number undetected errors: %d\n',sum(sum(out(interpolation_flags==0) ~= cd.scaled_quantized_padded_original(interpolation_flags==0))));

            % sound(out,Fs);
            audiowrite("output.wav",out,Fs);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function testStep1()
            % Simple test to see if decoding encoded input after step 1
            % equals the input.
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            frames = 10;
            input = round(rand(24*frames, 1)*128);
            flags = round(rand(24*frames, 1)*128);
            [encoded, n_frames] = CIRC_enc_delay_interleave(cd, input, frames);
            [erasure_flags_in, n_frames] = CIRC_enc_delay_interleave(cd, flags, frames);
            
            [decoded,erasure_flags_out,n_frames2] = CIRC_dec_deinterleave_delay(cd,encoded,erasure_flags_in,n_frames);
            
            correct = 1;
            for i=1:(n_frames2*24)
                if decoded(i) ~= input(i) || erasure_flags_out(i) ~= flags(i)
                    correct = 0;
                    break
                end
            end
            
            if correct == 0
                disp('Decoded does not equal input')
            else
                disp('Decoding after encoding succeeded!')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testStep2()
            % Simple test to see if decoding encoded input after step 2
            % equals the input.
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            n_frames = 10;
            input = round(rand(24*n_frames, 1)*128);
            err_frames = zeros(28*n_frames, 1, 'logical');
            [encoded, n_frames] = CIRC_enc_C2(cd, input, n_frames);
            [decoded, err_out, n_frames] = CIRC_dec_C2(cd, encoded, err_frames, n_frames);
            
            correct = 1;
            for i=1:(n_frames*24)
                if decoded(i) ~= input(i)
                    correct = 0;
                    break
                end
            end
            if correct == 0
                disp('Decoded does not equal input')
            else
                disp('Decoding after encoding succeeded!')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testStep3()
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            n_frames = 10;
            input = round(rand(28*n_frames, 1)*128);
            err = zeros(28*(n_frames+108), 1, 'logical');
            [delayed, n_frames_2] = CIRC_enc_delay_unequal(cd,input,n_frames);
            [undelayed, n_frames3] = CIRC_dec_delay_unequal(cd, delayed, err, n_frames_2);
            correct = 1;
            for i=1:(n_frames*28)
                if undelayed(i) ~= input(i)
                    correct = 0;
                    break
                end
            end
            if correct == 0
                disp('Decoded does not equal input')
            else
                disp('Decoding after encoding succeeded!')
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testC1()
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            n_frames = 1;
            input = round(rand(28*n_frames,1)*128);
            [enc,n_frames2] = CIRC_enc_C1(cd,input,n_frames);
            %enc(3) = 0; %1 error
            enc(3:4) = 0; %errors
            [output,erasure_flags_out,n_frames] = CIRC_dec_C1(cd,enc,n_frames2);
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [input,output,same] = testDecoderC2()
            %tests functionality of decoder with a few predefined errors
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            n_frames = 1;
            input = round(rand(24*n_frames, 1)*128);
            encoded = step(cd.enc2,input); %parity bits currently at the end
            encoded_pCenter = [encoded(1:12); encoded(25:28); encoded(13:24)];
            flags = zeros(28,1,'logical');
            flags(4:5) = 1; %set a few flags to 1
            encoded_pCenter(4:6) = 0;
            %encoded_pCenter(8) = 0;
            [output,flags_out,~] = CIRC_dec_C2(cd,encoded_pCenter,flags,n_frames);
            correct = 1;
            for i=1:n_frames*24
                if output(i) ~= input(i);
                    correct = 0;
                    break
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function testFinalInterleave()
            Fs = 22.05e3;
            cd = AudioCD_opt(Fs,1,8);
            n_frames = 4;
            input = uint8(round(rand(32*n_frames, 1)*128));
            [output_interleave,int_frames] = CIRC_enc_delay_inv(cd,input,n_frames);
            [output,n_frames] = CIRC_dec_delay_inv(cd,output_interleave,int_frames);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function evaluateScratchesAndRandomErrors()
            % Evaluate the results on a 20x100 log-log grid, varying
            % both the scratch width and the percentage of random errors.
            
            audio_file = audioread('Hallelujah.wav');
            Fs = 44.1e3;
            
            X = logspace(0, 5.60206, 100);
            Y = logspace(-3, -1, 20);
            Z11 = zeros(100,20);
            Z12 = zeros(100,20);
            Z13 = zeros(100,20);
            
            % Make use of a parallel pool to speed up processing.
            parfor i = 1:100
                for j = 1:20
                        cd = AudioCD_opt(Fs,1,8);
                        cd = cd.writeCd(audio_file);

                        T_scratch = 600000; % Scratch at a diameter of approx. 66 mm
                        for ii = 1:floor(numel(cd.cd_bits)/T_scratch)
                           cd = cd.scratchCd(round(X(i)),1+(ii-1)*T_scratch);
                        end

                        cd = cd.bitErrorsCd(Y(j));

                        [out,interpolation_flags] = cd.readCd();

                        nef = sum(sum(interpolation_flags~=0));
                        nfi = sum(sum(interpolation_flags==-1));
                        nue = sum(sum(out(interpolation_flags==0) ~= cd.scaled_quantized_padded_original(interpolation_flags==0)));

                        fprintf('i = %d, (%d, %d, %d)\n', i, nef, nfi, nue);

                        Z11(i,j) = nef;
                        Z12(i,j) = nfi;
                        Z13(i,j) = nue;
                    end
            end
            
            save('Z11_Final.mat','Z11');
            save('Z12_Final.mat','Z12');
            save('Z13_Final.mat','Z13');
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%