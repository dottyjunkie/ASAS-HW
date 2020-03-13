%%
%%1(a)
[x,fs] = audioread("sample.wav");

sound(x,fs);

%% 1(b)
temp=0;
p_size=100;
for i=p_size:264600
    for p=1:p_size
        temp=temp+x(i-p+1);
    end
    y(i)=temp/p_size;
    temp=0;
end
%%
sound(y,fs);
audiowrite('p=100.wav',y,fs)
%% 1(c)
y_conv = conv(x, 1/p_size*ones(p_size,1));
sound(y_conv,fs);
audiowrite('y_conv.wav',y,fs)
%% 1(d)
h = 1/p_size*ones(p_size,1);
freqz(h);

%% 1(e)
x_con_demo=[1,1,2,3,5]
x_con_demo_after=conv(x_con_demo, [1 -1]);

x_con_after=conv(x, [1 -1]);
%% 1(f)
noise=0.1*randn(44100,1);
soundsc(noise,fs);
audiowrite('noise.wav',noise,fs)
plot(psd(spectrum.periodogram,noise,'Fs',fs,'NFFT',length(noise)));
%%
noise_filter= conv(noise, 1/p_size*ones(p_size,1));
noise_filter=noise_filter(1:44100,1);
soundsc(noise_filter,fs);
audiowrite('noise_filter.wav',noise_filter,fs)
plot(psd(spectrum.periodogram,noise_filter,'Fs',fs,'NFFT',length(noise_filter)));
%% 2(a)
a=0.999;
y_iir(1)=1;
for i=2:264600
    y_iir(i)=a*y_iir(i-1)+x(i);
end
soundsc(y_iir,fs);
audiowrite('y_iir.wav',noise_filter,fs)
%% filter check
y_build_filter=filter([1],[1 -a],x);
soundsc(y_build_filter,fs);
audiowrite('y_build_filter.wav',y_build_filter,fs)
%% 2(b)
noise=0.1*randn(44100,1);
soundsc(noise,fs);
%% filter the noise by matlab buildin function

noise_buildin_filter=filter([1],[1 -a],noise);
soundsc(noise_buildin_filter,fs);
audiowrite('y_build_filter.wav',noise_buildin_filter,fs)
%% 2(c)
x_pulse = zeros(44100,1);
count=0
for i=1:44100
    count=count+1;
    if count==80
       x_pulse(i)=0.5;
       count=0;
    end
end
sound(x_pulse,fs);
audiowrite('x_pulse.wav',x_pulse,fs)
plot(psd(spectrum.periodogram,x_pulse,'Fs',fs,'NFFT',length(x_pulse)));
%% filter
x_pulse_filter=filter([1],[1 -a],x_pulse);
soundsc(x_pulse_filter,fs);
audiowrite('x_pulse_filter.wav',x_pulse_filter,fs)
plot(psd(spectrum.periodogram,x_pulse_filter,'Fs',fs,'NFFT',length(x_pulse_filter)));
%% 2(d)
freqz(filter([1],[1 -0],x_pulse));
%%
freqz(filter([1],[1 -0.2],x_pulse));
%%
freqz(filter([1],[1 -0.4],x_pulse));
%%
freqz(filter([1],[1 -0.6],x_pulse));
%%
freqz(filter([1],[1 -0.8],x_pulse));
%%
freqz(filter([1],[1 -0.99],x_pulse));
%%
%%
freqz(filter([1],[1 -1],x_pulse));