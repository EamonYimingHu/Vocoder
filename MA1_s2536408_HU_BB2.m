%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Beyond basics 
%
% Phase-vocoder with pitch shifting
%
% I've noticed that the pitch shift combined with the previous phase vocoder can make it a voice changer
% To use as a voicce changer, ensure pitch_shift_factor*Q=1,(e.g I've pre-set this as 3 tone lower,for female voice changing to a male voice)
% if set pitch_shift_factor =1,play with Q, it's phase vocoder,changing speed without changing pitch
% if set Q=1, play with pitch_shift_factor, it will change both speed and pitch
%
% Yiming HU 20/11/23
%-------------------------------------------------


% Clear the command window, workspace and close all plots -----------------
clc;                            % clear the command window
clear;                          % clear workspace
close all;                      % close all plots

% Read in the input WAV file
[x, Fs] = audioread('Cath_cut.wav');
%Average left and right channels to mono in case the input audio is stereo.
if size(x,2) == 2
    x=mean(x,2);
end


% Define the anaLysis frame time length in seconds
frame_time = 0.03;
O = 0.75;                       %set overlap factor
Q = 1/2 ^ (6/12);
pitch_shift_factor = 2^(6/12);  %Here I introduced this to normalise Ra so that the sound can be reconstructed when factor&Q==1
%My pre-set is 2^(6/12) stands for 3 tone lower, changing originial female voice to male voice
%Set pitch shift factor<1, for higher pitch; factor>1, for lower pitch

% Convert frame time length to samples
N = round(frame_time * Fs);     %def N as the length of analysis frame
HA = round((1 - O) * N);        %def anaLysis hop size from N and O
HS = round(Q * HA);             %Calculate hop size HS from Q
RA = pitch_shift_factor * HS;

% Generate a Hann window of length N manualLy
win = 0.5 * (1 - cos(2 * pi * (0:N - 1).' / N)); 
win(1) = 0;                     % Make sure the first value is zero to make it periodic

% Note: The window is already periodic as the last value will not be zero

% firstLy Zero-padding at the beginning---------------
x = [zeros(N, 1); x];
L = length(x);


% calculate exactly how many samples you need to pad the end, and zero-padding at the end----------------------
end_padding = N -mod(L - N, HA);
%For the special case that mod=0, which means there is no need for endpadding
if end_padding == N
    end_padding = 0;
end
x = [x; zeros(end_padding, 1)];
% Now x has been zero-padded at both the beginning and the end, The length of x is a perfect fit for samples
L = length(x);
NF = floor((L - N) / HA) + 1;   %number of frames


% Q3 Calculate NFFT-----------------------------------
NFFT = 2 ^ nextpow2(N);

% Q5 Initialize variables for phase vocoder-----------
phi_m = zeros(NFFT / 2 + 1, 1);
phi_m_plus_1 = zeros(NFFT / 2 + 1, 1);
theta_m = zeros(NFFT / 2 + 1, 1);
theta_m_plus_1 = zeros(NFFT / 2 + 1, 1);
omega_hat_k = (2 * pi * (0:NFFT / 2)') / NFFT;

% Create and initialize output vector y--------------
y = zeros(floor((NF - 1) * RA) + floor(N * RA / HS), 1);

% Extract the m-th frame
for m = 1:NF - 1
    xm = win .* x(m * HA + 1:m * HA + N); 

    % Q7Apply FFT on every frame with length NFFT--------
    X = fft(xm, NFFT);
    %For simplicity, use fft’s in-built functionality, i.e if length parameter is greater than the vector fft will automatically pad the vector up with zero

    % Q8Keep only the first NFFT/2 + 1 bins and Separate X into magnitude and phase--------------
    Xmag = abs(X(1:NFFT / 2 + 1));
    Xang = angle(X(1:NFFT / 2 + 1));

    % Calculate the phase differences and wrap to pi
    phi_m_plus_1 = angle(X(1:NFFT / 2 + 1));
    deltaphi = ppa(phi_m_plus_1 - phi_m - omega_hat_k * HA);
    phi_m = phi_m_plus_1;

    % Q9Estimate the instantaneous frequencies of the current frame-----------------------------------
    IF_m_plus_1 = omega_hat_k + deltaphi / HA;                        %IF(instantaneous frequencies)

    % Q10Calculate the modified phases
    theta_m_plus_1 = ppa(theta_m + IF_m_plus_1 * HS);
    theta_m = theta_m_plus_1;

    % Q11Create and Calculate new modified DFT frame Ym+1[k] of full-length NFFT
    Ym_plus_1_half = Xmag .* exp(1j * theta_m_plus_1);
    Ym_plus_1 = zeros(NFFT, 1);                                       %initialise Ym+1[k]
    Ym_plus_1(1:NFFT / 2 + 1) = Ym_plus_1_half;
    Ym_plus_1(NFFT / 2 + 2:end) = conj(Ym_plus_1_half(end - 1:-1:2)); %fill in the remaining using Hermitian symmetry

    %ensure the DC bin and Nyquist bin are real,now we have Hermitian symmetric, in next step the ifft result are then pure real
    Ym_plus_1(NFFT / 2 + 1) = real(Ym_plus_1(NFFT / 2 + 1));
    Ym_plus_1(1) = real(Ym_plus_1(1));

    % Q12Take the inverse DFT
    y_modified = ifft(Ym_plus_1, NFFT);

    % Q13Truncate the frame from NFFT to N samples and apply synthesis window-------
    y_modified = y_modified(1:N);
    y_modified = y_modified .* win;
    % Apply resampling pitch shifting and generate resampling window
    t_original = (0:1 / Fs:(N - 1) / Fs);
    t_resampled = (0:(1 / (Fs * RA / HS)):(N*RA/HS - 1) / (Fs * RA / HS));
    %t_resampled = (0 : 1/Fs: (N*RA/HS-1)/Fs);%Use this if you want to keep Fs all the time, same feature but lower quality
    resampledSignal = interp1(t_original, y_modified, t_resampled, 'linear', 'extrap');
    %after trying spline,pchip,linear. the difference was not significant, so we choose linear to improve efficiency

    y((1:floor(N * RA / HS)) + floor(m * RA)) = y((1:floor(N * RA / HS)) + floor(m * RA)) + resampledSignal';
    % Add the windowed synthesis frame into the output vector y
end

sound(y, Fs);

%function ppa: wrap any phase angle to phase over[-pi,pi]
function wrapped_phase = ppa(phase)
    wrapped_phase = mod(phase + pi, 2 * pi) - pi;
end

