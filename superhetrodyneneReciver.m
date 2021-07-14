clc;
close all;
clear all;
%Transmitter Tx
%reading the 4 signals determine Fs
[M1,Fs1] = audioread('Short_BBCArabic2.wav');%M1 represents The first message
%(first channel)(modulating signal 1)
[M2,Fs2] = audioread('Short_FM9090.wav');%M2 represents the second message(second
%channel)(modulating signal 2)and so on
[M3,Fs3] = audioread('Short_QuranPalestine.wav');
[M4,Fs4] = audioread('Short_SkyNewsArabia.wav');
%converting Them to 1_d vector in order to use interp function
M1=M1'; %Transpose The Data Matrix
M2=M2';
M3=M3';
M4=M4';
M1(2,:) = [];% to eliminate the 2nd CHANNEL vector
M2(2,:) = [];
M3(2,:) = [];
M4(2,:) = [];

% making the 4 signals have the same length (padding by zeroes)
n1=length(M1);
n2=length(M2);
n3=length(M3);
n4=length(M4);
for i=n2:n1
M2(i)=0;
end
for i=n3:n1
M3(i)=0;
end
for i=n4:n1
M4(i)=0;
end
%PLOT THE 4 SIGNALS TO CHECK NYQUIST CRITERIA and to get THE BW of each %message
N=length(M1);
f=(-N/2:N/2-1)*Fs1/N; %Frequency vector (horizontal axis)
MW1=fft(M1,N); %MW1 represents the fourier transform of the first channel
MW1=fftshift(MW1); % To make the spectrum centered at zero
figure(1)
subplot(4,1,1);
plot(f,abs(MW1))
xlabel('F(Hz)')
ylabel('|M1(W)|(Volt)')
title('The Spectrum of M1')
MW2=fft(M2,N);%MW2 represents the fourier transform of the second channel
MW2=fftshift(MW2);
subplot(4,1,2);
plot(f,abs(MW2))
xlabel('F(Hz)')
ylabel('|M2(W)|(Volt)')
title('The Spectrum of M2')
MW3=fft(M3,N);%MW3 represents the fourier transform of the third channel
MW3=fftshift(MW3);
subplot(4,1,3);
plot(f,abs(MW3))
xlabel('F(Hz)')
ylabel('|M3(W)|(Volt)')
title('The Spectrum of M3')
MW4=fft(M4,N);%MW3 represents the fourier transform of the fourth channel
MW4=fftshift(MW4);
subplot(4,1,4);
plot(f,abs(MW4))
xlabel('F(Hz)')
ylabel('|M4(W)|(Volt)')
title('The Spectrum of M4')
%Multiply FS by 15 to guarantee normalized freq of RX_filters will be within[0,1]
M1=interp(M1,15);%Multiply FS by 15 also to guarantee nyquist criteria
M2=interp(M2,15);
M3=interp(M3,15);
M4=interp(M4,15);
Fs=15*Fs1;
n1=length(M1);
%AM MODULATOR DSB-Sc
for i=1:n1
ModulatedCarrier1(i)=M1(i)*cos(2*pi*(100000/Fs)*i);
ModulatedCarrier2(i)=M2(i)*cos(2*pi*(150000/Fs)*i);

ModulatedCarrier3(i)=M3(i)*cos(2*pi*(200000/Fs)*i);
ModulatedCarrier4(i)=M4(i)*cos(2*pi*(250000/Fs)*i);
end
%channel
TX_output=ModulatedCarrier1+ModulatedCarrier2+ModulatedCarrier3+ModulatedCarrier4;
%TX_output represents the sum of the modulated signals(FDM)(output of Tx)
%PLOT The Spectrum of The Output of The transmitter
N=length(TX_output);
f=(-N/2:N/2-1)*Fs/N;
TX_output_W=fft(TX_output,N);%TX_output_W represents the Fourier Transform of
%Tx_out
TX_output_W=fftshift(TX_output_W);
figure(2)
plot(f,abs(TX_output_W))
xlabel('F(Hz)')
ylabel('|TXoutput(W)|(Volt)')
title('The Spectrum of The Output of The transmitter')
%Receiver Rx
% RF stage (rf bandpass filter + local oscillator) Tunable
n=input("enter the number of the channel you want to hear [0,1,2,3] ");
f=100000+n*50000;
% RF band Pass Filter(To reject image signal) (Tunable)
[b,a]=butter(5,[((f-12500)*2)/(Fs),((f+12500)*2)/(Fs)]);
%(12500<<2Wif) it wont be Sharp so minimize BW to reject The image signal
%(All Frequencies are normalized to FS/2)
RFFilter_OutPut=filter(b,a,TX_output);
%PLOT The Spectrum of The Output of The Rf filter
N=length(RFFilter_OutPut);
f1=(-N/2:N/2-1)*Fs/N;
RFFilter_OutPut_w=fft(RFFilter_OutPut,N);
RFFilter_OutPut_w=fftshift(RFFilter_OutPut_w);
figure(3)
plot(f1,abs(RFFilter_OutPut_w))
xlabel('F(Hz)')
ylabel('|RFFilterOutPut(W)|(Volt)')
title('The Spectrum of The Output of The RF Filter')
%local oscillator Wlo=Wrf-Wif (Tunable)(TO convert the RF_signal TO if_signal)
for i=1:n1
RFMixer_output(i)= RFFilter_OutPut(i)*10*cos(2*pi*((f/Fs)-(25000/Fs))*i);
end
%PLOT The Spectrum of The Output of Mixer
N=length(RFMixer_output);
f=(-N/2:N/2-1)*Fs/N;
RFMixer_output_W=fft(RFMixer_output,N);
RFMixer_output_W=fftshift(RFMixer_output_W);
figure(4)
plot(f,abs(RFMixer_output_W))
title('The Spectrum of The Output of The RF Mixer')
xlabel('F(Hz)')
ylabel('|RFMixeroutput(W)|(Volt)')
%if stage (Fixed)(To get The IF_signal only)
[b,a]=butter(6,[((25000-9000)*2)/(Fs),((25000+9000)*2)/(Fs)]);
IFFilter_Output=filter(b,a,RFMixer_output);
%PLOT The Spectrum of The Output of the if stage
N=length(IFFilter_Output);
f=(-N/2:N/2-1)*Fs/N;
IFFilter_Output_W=fft(IFFilter_Output,N);

IFFilter_Output_W=fftshift(IFFilter_Output_W);
figure(5)
plot(f,abs(IFFilter_Output_W))
title('The Spectrum of The Output of the if stage')
xlabel('F(Hz)')
ylabel('|IFFilterOutput(W)|(Volt)')
%if Mixer (Fixed) (To convert The if_signal To baseband_signal )
for i=1:n1
IFMixer_output(i)=IFFilter_Output(i)*10*cos(2*pi*(25000/Fs)*i);
end
%PLOT The Spectrum of The Output of the if Mixer
N=length(IFMixer_output);
f=(-N/2:N/2-1)*Fs/N;
IFMixer_output_w=fft(IFMixer_output,N);
IFMixer_output_w=fftshift(IFMixer_output_w);
figure(6)
plot(f,abs(IFMixer_output_w))
title('The Spectrum of The Output of the if Mixer')
xlabel('F(Hz)')
ylabel('|IFMixeroutput(W)|(Volt)')
%LPf stage (Fixed)(To get the Baseband_signal only)
[b,a]=butter(6,2*9000/Fs);
LPf_out=filter(b,a,IFMixer_output);
%PLOT The Spectrum of The Output of the LPf stage
N=length(LPf_out);
f=(-N/2:N/2-1)*Fs/N;
LPf_out_w=fft(LPf_out,N);
LPf_out_w=fftshift(LPf_out_w);
figure(7)
plot(f,abs(LPf_out_w))
title('The Spectrum of The Output of the LPf stage')
xlabel('F(Hz)')
ylabel('|LPfout(W)|(Volt)')
%normalizing the data between[0,1] to avoid clipping
LPf_out=LPf_out/max(abs(LPf_out));
%divide FS by 15 in order to use sound function
RX_Out=decimate(LPf_out,15);
sound(RX_Out,44100)