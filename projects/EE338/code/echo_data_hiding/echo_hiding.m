close all;
clearvars;
clc;

% Fixed parameters

d0 = 150;     %Delay rate for bit0
d1 = 200;     %Delay rate for bit1
amp = 0.5;  %Echo amplitude
L = 8*1024;  %Length of frames


%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  ENCODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[signal,Fs] = audioread('sounds/cello-double.wav');
% sound(signal,Fs)

[I,ch] = size(signal);

file = 'data.txt';
fid  = fopen(file, 'r');
text = fread(fid,'*char')';
matrix  = dec2bin(uint8(text),8);
msg = reshape(matrix', 1, 8*length(text));
fprintf('Encoded Text: %s\n', text);
nframe = floor(I/L);
N = nframe - mod(nframe,8);

ker0 = [zeros(d0, 1); 1]*amp;        %Echo kernel for bit0
ker1 = [zeros(d1, 1); 1]*amp;        %Echo kernel for bit1

echo_zero = filter(ker0, 1, signal);    %Echo signal for bit0
echo_one = filter(ker1, 1, signal);    %Echo signal for bit1

if (length(msg) > N)
    fprintf('Message is too long, being cropped!\n');
    bits = msg(1:N);
else
    bits = [msg, num2str(zeros(N-length(msg), 1))'];
end
K=256;
K = K - mod(K, 4);                       %Divisibility by 4
encbit = str2num(reshape(bits, N, 1))';      %char -> double
m_sig  = reshape(ones(L,1)*encbit, N*L, 1);  %Mixer signal
c      = conv(m_sig, hanning(K));            %Hann windowing
w_sig  = c(K/2+1:end-K/2+1) / max(abs(c));   %Normalization
figure(3),subplot(2,1,1);
plot(w_sig);
ylim([-1 2]);
figure(3),subplot(2,1,2);
plot(1-w_sig);
ylim([-1 2]);
mix = w_sig;  %Mixer signal
encoded = signal(1:N*L, :) + echo_zero(1:N*L, :) .* abs(mix-1) + echo_one(1:N*L, :) .* mix;
encoded = [encoded; signal(N*L+1:I, :)];   %Rest of the signal
figure(2),subplot(3,1,1)
u=5000;
l=8000;
plot(u:l, signal(u:l, :))
subplot(3,1,2)
plot([u:l],echo_zero(u:l))
subplot(3,1,3)
plot([u:l],echo_one(u:l))
audiowrite('sounds/outputs/decoded.wav',encoded,Fs);
% sound(encoded,Fs)
%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECODING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% len_msg = 0;

N = floor(length(encoded)/L);             %Number of frames
segments = reshape(encoded(1:N*L,1), L, N);   %Dividing signal into frames
data = char.empty(N, 0);

for k=1:N
	rceps = ifft(log(abs(fft(segments(:,k)))));  %Real cepstrum
	if (rceps(d0+1) >= rceps(d1+1))
        data(k) = '0';
	else
        data(k) = '1';
	end
end

m   = floor(N/8);
bin = reshape(data(1:8*m), 8, m)';   %Retrieved in binary
msg = char(bin2dec(bin))';           %bin=>char

% if (len_msg~=0)
% 	msg = msg(1:len_msg);
% end

fprintf('Decoded Text: %s\n', msg);



