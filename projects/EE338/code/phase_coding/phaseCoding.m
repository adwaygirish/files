close all; clear; clc;

%% encoding
% write path to base song file
[signal, FS] = audioread('in.wav');

file = 'message.txt'; % path to file with the message that is to be hidden
fid  = fopen(file, 'r');
text = fread(fid,'*char')';
fclose(fid);

success = false;

base = signal(:,1);
L_msg = length(text);
msg = getBits(text);
I = length(base);
m = length(msg);    % Length of the bit sequence to be hidden
if I <= 2*m
    disp(['Text message too long - find a longer base signal']);
    out = 0;
    return;
end

L = min(2^nextpow2(2*m), I); % length of each segment
if L == 2*m
    L = L + 2;
    if L >= I
        disp(['Text message too long - find a longer base signal']);
        out = 0;
        return;
    end
end
halfL = floor(L/2);
N = floor(I/L);        % Number of points per segment

% Segment the signal
s = reshape(base(1:N*L,1), L, N);

% obtain L-point DFT of each signal, and the phase and amplitudes
w = fft(s);
Phi = angle(w);        % Phases - N columns, each of L rows
A = abs(w);             % Amplitudes

% Calculating phase differences between adjacent segments
DeltaPhi = zeros(L,N);
for i=2:N
    DeltaPhi(:,i)=Phi(:,i)-Phi(:,i-1);
end

% Convert the data bits into phase {'0' : pi/2, '1' : -pi/2}
PhiData = zeros(1, m);
phi0 = pi/2;
for i=1:m
    if msg(i) == '0'
        PhiData(i) =phi0;
    else
        PhiData(i) = -phi0;
    end
end

Phi_new = zeros(size(DeltaPhi));

% Save message in column 1 of the phase matrix
% remember to maintain hermitian symmetry
Phi_new(:,1) = Phi(:,1);
Phi_new(halfL-m+1:halfL,1) = PhiData;             % Hermitian symmetry
Phi_new(halfL+1+1:halfL+1+m,1) = -flip(PhiData);  % Hermitian symmetry

% Constructing phase matrix using phase differences DeltaPhi
for i=2:N
    Phi_new(:,i) = Phi_new(:,i-1) + DeltaPhi(:,i);
end

% Reconstructing the sound signal by taking IFFT
z = real(ifft(A .* exp(1i*Phi_new)));
out = reshape(z, N*L, 1);               % concatenating
if N*L + 1 <= I
    % Adding rest of signal - I/L may not be an integer
    out  = [out; base(N*L+1:I)];
end
success = true;

if ~success
    disp(['Failed. Try with another song/message']);
    return;
end
audiowrite(['out.wav'], out, FS);
disp(['Stego signal is saved in ', 'out_stego.wav']);

%% decoding

% % Uncomment the  lines 94-97 if you want to test effect of noise,
% % else comment
% disp(['********Adding noise********']);
% out = awgn(out, 60);
% audiowrite('out_noisy.wav', out, FS);
% disp(['Noisy stego signal is saved in ', 'out_noisy_stego.wav']);

modifiedSignal = out;

m_decode   = 8*L_msg;             % Length of bit sequence (assuming 8bits per character)
I_decode = length(modifiedSignal);
L_decode = min(2^nextpow2(2*m_decode), I_decode); % length of each segment
if L_decode == 2*m_decode
    L_decode = L_decode + 2;
end
halfL_decode = floor(L_decode/2);
x   = modifiedSignal(1:L_decode,1);       % First segment
Phi_decode = angle(fft(x));       % Phase angles of first segment

% Retrieving data back from phases stored in first segment
data_decode = char(zeros(1,m_decode));
for k=1:m_decode
	if Phi_decode(halfL_decode - m_decode + k)>=0
    	data_decode(k)='0';
    else
        data_decode(k)='1';
	end
end

bin_decode = reshape(data_decode(1:m_decode), 8, m_decode/8)';
out_decode = char(bin2dec(bin_decode))';

%% Check error rate
ratio = BER(out_decode, text); % Bit error rate
disp(['******** RESULTS ********']);
fprintf('Original Text: %s\n', text); 
fprintf('Decoded Text: %s\n', out_decode); 
fprintf('BER : %d\n', ratio);

function out = BER(a, b)
    %BER Bit Error Rate
    a_bits = getBits(a); 
    b_bits = getBits(b);
    len_a = length(a_bits);
    len_b = length(b_bits);
    len = min(len_a, len_b);
    ber = 0;
    for i=1:len
        ber = ber +  (a_bits(i) ~= b_bits(i));
    end
    ber = ber + abs(len_b - len_a);
    out = 100*(ber/len);
end

function bit_seq = getBits(text)
    % we are assuming that we store txt - char data type
    matrix  = dec2bin(uint8(text),8);
    bit_seq = reshape(matrix', 1, 8*length(text));
end