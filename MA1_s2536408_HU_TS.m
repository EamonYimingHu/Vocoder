%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Part One
%
% A basic time stretcher
% For Q=0.75;suitable frame_time=0.057（N=2514）O=0.75(HA=629)
% For Q=1.25;suitable frame_time=0.07  (N=2558) O=0.8(HA=617),pre-set this already
%
% Yiming HU 11/11/23 
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                            % clear the command window
clear;                          % clear workspace
close all;                      % close all plots

% Read in the input WAV file
[x, Fs] = audioread('Cath_cut.wav');
% Average left and right channels to mono in case the input audio is stereo.
if size(x, 2) == 2
    x = mean(x, 2);
end


% Define the anaLysis frame time length in seconds
frame_time = 0.07;
O = 0.8;                        %set overlap factor
Q = 1.25;

% Convert frame time length to samples
N = round(frame_time * Fs);     %def N as the length of analysis frame
HA = round((1 - O) * N);        %def anaLysis hop size from N and O
HS = round(Q * HA);             %Calculate hop size HS from Q

if Q==1                         %This is for Q10 reconstruction
    N = 4 * HA;                 %We sightly adjust N to make N/HA a integer for O=0.75
end

% Generate a Hann window of length N manualLy
win = 0.5 * (1 - cos(2 * pi * (0:N - 1).' / N)); 
win(1) = 0;                     % Make sure the first value is zero to make it periodic
% Note: The window is already periodic as the last value will not be zero


% firstLy Zero-padding at the beginning---------------
x = [zeros(N, 1); x];
L=length(x);
% calculate exactly how many samples you need to pad the end, and zero-padding at the end----------------------
end_padding = N -mod(L-N,HA);
if end_padding == N             %For the special case that mod=0, which means there is no need for endpadding
    end_padding = 0;
end
x = [x; zeros(end_padding, 1)];
% Now x has been zero-padded at both the beginning and the end, The length of x is a perfect fit for samples 
L=length(x);
NF = floor((L - N) / HA) + 1; 


% Q7 Create and initialize output vector y--------------
y = zeros((NF - 1) * HS + N, 1);

% Q8 read in mth analysis frame from x------------------
for m = 1:NF-1  
    y(m * HS + 1:m * HS + N) = y(m * HS + 1:m * HS + N) +x(m * HA + 1:m * HA + N) .* win .* win;
end
% mth analysis frame start point is m * HS + 1，end at m*HS+N which is exact N sample


% Q10Calculate the gain introduced by the analysis and synthesis windows---------------------------
gain = sum(win .^ 2)/ HA;

% Compensate for the gain in the output
y = y / gain;

% Plot the input and output (overlaid), and the error
if Q == 1
    t = (0:length(x) - 1) * (1/Fs); % Time vector
    % Create a new figure
    figure;

    % First subplot: Input and Output overlaid
    subplot(2, 1, 1);
    plot(t, x,'r', t, y,'b');
    title('Input and Output Signals');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Input', 'Output');
    grid on;

    % Second subplot: Error
    subplot(2, 1, 2);
    plot(t, y-x, 'm');
    title('Error Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    legend('Error');
    grid on;
end

% Play the compensated output
soundsc(y, Fs);

% Call plot spectrogram function
MA1_s2536408_HU_myspec(x, Fs, N, O)

