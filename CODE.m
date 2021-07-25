%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PART 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fo=(1e4)/(2*pi);    %freq of the signal
fs=8e3;             %sampling rate
ts=1/fs;            %sampling interval
T=50e-6;            %Pulse width interval
RFACT = 25e-4;      %Range factor to control the range

t = 0:1/(100*fs):(2*RFACT-(1/(100*fs)));

mt=sin(2*pi*fo*t);                  %input signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              PAM SIGNAL IN TIME DOMAIN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberofpulses = 2*RFACT*fs;        %number of pulses
duty = 40;                          %duty cycle of the train of the pulses

v = length(t);   % vector of time to manipulate the error of the matrix dimensions

samples = v/numberofpulses;               %number of samples
%number samples = 5e-4/ts;
pulseindex = [1:samples:v];
PWM_HIGH =  ceil(samples * duty/100);    % ceil the values of the falt topp rectangular pulse
PAMsignal = zeros(1,v);
for i = 1 : length(pulseindex)           %for loop to draw the PAM signal 
    PAMsignal(pulseindex(i):pulseindex(i) + PWM_HIGH) = mt(pulseindex(i));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            PLOTTING PAM SIGNAL IN TIME DOMAIN 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(t,PAMsignal,'r')
hold on;
plot(t,mt,'b');
hold on;
ylabel('Amplitude')
xlabel('Time')
%title('PAM signal in Time Domain')
%legend('s(t)','m(t)')
grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            PLOTTING PAM SIGNAL IN FREQEUNCY DOMAIN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%fourier transform and frequency domain
N=length(PAMsignal);  %define Vector of the same Length of N
%To plot the frequency specterum the frequency and time should be in the
%same length
%f=0:1:N-1 -(N-1)/2:1/N:(N-1)/2  -fs*(N-1)/2:fs/N:fs*(N-1)/2  
f=(-fs/2):(fs/N):((fs/2)-(fs/N)); %Define frequency axis to draw the spectrum on it shifting it 
%-fs/N to remove one iteration which resemples the zero value to get the
%same length of the array of the fourier transform
SF=fftshift(fft(PAMsignal));
%plot frequency spectrum
figure(2)
plot(f,((10/N*7))*abs(SF),'b') % ((10/N*7)) to normalize the signal  
%title('Frequency Spectrum');
xlabel('f(Hz)');
ylabel('PAM(F)');
%legend('Spectrum of PAM(t)');
grid

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%PART 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Decleration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs=50;
Ts=1/fs;
step=0.2;
t=0:Ts:2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                 I/P Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
inputSignal=2*sin(2*pi*t);
%plotting the input signal
figure(3)
plot(inputSignal,'b');
% title('Input Signal')
ylabel('Amplitude')
xlabel('Time')
xlim([0 100]); %2 cycles of the sin
grid
hold on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            Delta Modulated Signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bitstream = zeros(1,101);
%intialization of bit stream to store the value of DM PCM code
Stair_Case_approx= zeros(1,101);
%intialization of the strais array to hold the value to be compared with
%value of the incoming input signal

length_ip_signal=length(inputSignal);
%vector array to hold the length the input signal depending on the time
%range stated at the Decleration top then this length is looped on through
%for loop that loops on all the input signal value and hold value of
%bitstream and stair case depending on the relative postion between each
%two signals at the iteration of the for loop
for i=1:length_ip_signal
    %if the signal value is higher than the stair case hold 1
    if inputSignal(i)>Stair_Case_approx(i)
        Bitstream(i)=1;
        Stair_Case_approx(i+1)=Stair_Case_approx(i)+step;
    else
    %if the signal value is higher than the stair case hold 0
        Bitstream(i)=0;
        Stair_Case_approx(i+1)=Stair_Case_approx(i)-step;
    end
end

%plotting the value of Stair_Case_approx
stairs(Stair_Case_approx);
grid
hold on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             Recovered signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%for loop that loops on all the bitstream values and compare it to the
%value of the stair case to recover the signal 
for i=1:Bitstream
    if Bitstream(i)>Stair_Case_approx(i)
        Bitstream(i)=0;
        Stair_Case_approx(i+1)=Stair_Case_approx(i)-step;
    else
        Bitstream(i)=1;
        Stair_Case_approx(i+1)=Stair_Case_approx(i)+step;
    end
end
%plotting the recovered signal
 plot(Stair_Case_approx,'r');
legend ('Input Signal m(t)','Delta Modulated Signal s(t)','Recovered Signal');
title('Time domain signals')
