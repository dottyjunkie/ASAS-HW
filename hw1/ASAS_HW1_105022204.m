%%1
y=zeros(16000,1);
x=zeros(16000,1);
f0=4400;
T=1/16000;
fs=16000;
L=16000
M=8000;

for i = 1:16000
    
	y(i)=0.5*g(i-8000,M)*cos(2*pi*f0*i*T);
    x(i)=i;
end
sound(y, fs);
audiowrite("sound_a_long_wavepacket_high.wav",y,fs)
plot(x,y)


%% 2(a) test part of matlab fft
fft_y=fft(y);
fft_y = 20*log10(fft_y + 1e-6);
P2 = abs(fft_y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;

plot(f,P1) 
%% my version of dft 
my_dft_freq=my_dft(y);
my_dft_freq=my_dft_freq.';
%% draw it
my_dft_freq = 20*log10(my_dft_freq + 1e-6);
P2 = abs(my_dft_freq/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f,P1) 
%% 3
[y_sample,fs_sample] = audioread('sample.mp3');

wlen = 1024;                        % window length (recomended to be power of 2)
hop = wlen/4;                       % hop size (recomended to be power of 2)
nfft = 4096;                        % number of fft points (recomended to be power of 2)
win = hann(wlen, 'periodic');
%% my version of stft
[my_s,my_f,my_t]=my_stft(y_sample, win, hop, nfft, fs_sample);
C = sum(win)/wlen;
my_s = abs(my_s)/wlen/C;

if rem(nfft, 2)                    
    my_s(2:end, :) = my_s(2:end, :).*2;
else                                
    my_s(2:end-1, :) = my_s(2:end-1, :).*2;
end

my_s = 20*log10(my_s + 1e-6);
%%
figure(1)
surf(my_t, my_f, my_s)
shading interp
axis tight
view(0, 90)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Amplitude spectrogram of the signal')
hcol = colorbar;
set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(hcol, 'Magnitude, dB')
%% use matlab build in spectrogram
[s_buildin, f_buildin, t_buildin] = spectrogram(y_sample(:,1), win, hop, nfft, fs_sample);
C = sum(win)/wlen;
s_buildin = abs(s_buildin)/wlen/C;

if rem(nfft, 2)                    
    s_buildin(2:end, :) = s_buildin(2:end, :).*2;
else                                
    s_buildin(2:end-1, :) = s_buildin(2:end-1, :).*2;
end

s_buildin = 20*log10(s_buildin + 1e-6);
%%
figure(1)
surf(t_buildin, f_buildin,s_buildin)
shading interp
axis tight
view(0, 90)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
xlabel('Time, s')
ylabel('Frequency, Hz')
title('Amplitude spectrogram of the signal from Matlab')
hcol = colorbar;
set(hcol, 'FontName', 'Times New Roman', 'FontSize', 14)
ylabel(hcol, 'Magnitude, dB')
%


%% fuction part
function output = g(i,M)
    if i>=(-M) && i<=(M)
        output=cos(pi/2*i/M)^2;
        %output=1
    else
        output=0;
    end
end
function [STFT, f, t] = my_stft(x, win, hop, nfft, fs)
% Input:
% x - signal in the time domain
% win - analysis window function
% hop - hop size
% nfft - number of FFT points
% fs - sampling frequency, Hz
%
% Output:
% STFT - STFT-matrix (only unique points, time 
%        across columns, frequency across rows)
% f - frequency vector, Hz
% t - time vector, s
x = x(:); 
xlen = length(x);
wlen = length(win);
% stft matrix size estimation and preallocation
NUP = ceil((1+nfft)/2);     % calculate the number of unique fft points
L = 1+fix((xlen-wlen)/hop); % calculate the number of signal frames
STFT = zeros(NUP, L);       % preallocate the stft matrix
% STFT (via time-localized FFT)
for l = 0:L-1
    % windowing
    xw = x(1+l*hop : wlen+l*hop).*win;
    % FFT
    X = fft(xw, nfft);
    % update of the stft matrix
    STFT(:, 1+l) = X(1:NUP);
end
% calculation of the time and frequency vectors
t = (wlen/2:hop:wlen/2+(L-1)*hop)/fs;
f = (0:NUP-1)*fs/nfft;
end

function xk = my_dft(xn)
    ln=length(xn);
    i=sqrt(-1);
    xk=zeros(1,ln);
    
    for k=0:ln-1
        for n=0:ln-1
            xk(k+1)=xk(k+1)+(xn(n+1)*exp((-i)*2*pi*k*n/ln));
        end
    end
end