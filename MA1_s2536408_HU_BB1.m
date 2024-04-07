%-------------------------------------------------
% MAFTDSP Matlab Assignment 1 - Beyond Basics
%
% Phase-vocoder with time variant Q(t)
%
% I made Q as a oscillatory function of time using cos function,and adjust to appropriate range [0.8,1.2] to fit audio
% I made Q as an exponential function of time as well
% Q is at line 60 for cos, and at line 61 for exponential
%
% Yiming HU 11/11/23
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
frame_time = 0.02;
O = 0.75;                       %set overlap factor


% Convert frame time length to samples
N = round(frame_time * Fs);     %def N as the length of analysis frame
HA = round((1 - O) * N);        %def anaLysis hop size from N and O


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
NF = floor((L-N)/HA)+1;         %number of frames


% Define Q as a oscillatory function of time, adjust to appropriate range [0.8,1.2] to fit audio
Q = 1 + 0.2 * cos( pi * (0:NF - 1) / NF);
%Q = exp( -1 * (0:NF - 1) / NF);
HS = round(Q * HA);             %Calculate hop size HS from Q

% Q3 Calculate NFFT-----------------------------------
NFFT = 2 ^ nextpow2(N);

% Q5 Initialize variables for phase vocoder-----------
phi_m = zeros(NFFT / 2 + 1, 1);
phi_m_plus_1 = zeros(NFFT / 2 + 1, 1);
theta_m = zeros(NFFT / 2 + 1, 1);
theta_m_plus_1 = zeros(NFFT / 2 + 1, 1);
omega_hat_k = (2 * pi * (0:NFFT / 2)') / NFFT;


y = zeros(400512, 1); 
%400512 is the y-length for Q=1.5, Pre-allocate an maxium length

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
    IF_m_plus_1 = omega_hat_k + deltaphi / HA;                       %IF(instantaneous frequencies)

    % Q10Calculate the modified phases
    theta_m_plus_1 = ppa(theta_m + IF_m_plus_1 * HS(m));
    theta_m = theta_m_plus_1;

    % Q11Create and Calculate new modified DFT frame Ym+1[k] of full-length NFFT
    Ym_plus_1_half = Xmag .* exp(1j * theta_m_plus_1);
    Ym_plus_1 = zeros(NFFT, 1);                                       %initialise Ym+1[k]
    Ym_plus_1(1:NFFT / 2 + 1) = Ym_plus_1_half;
    Ym_plus_1(NFFT / 2 + 2:end) = conj(Ym_plus_1_half(end - 1:-1:2)); % fill in the remaining using Hermitian symmetry

    %ensure the DC bin and Nyquist bin are real,now we have Hermitian symmetric, in next step the ifft result are then pure real
    Ym_plus_1(NFFT / 2 + 1) = real(Ym_plus_1(NFFT / 2 + 1));
    Ym_plus_1(1) = real(Ym_plus_1(1));

    % Q12Take the inverse DFT
    y_modified = ifft(Ym_plus_1, NFFT);

    % Q13Truncate the frame from NFFT to N samples and apply synthesis window-------
    y_modified = y_modified(1:N);
    y_modified = y_modified .* win;

    % Add the windowed synthesis frame into the output vector y
    y(m * HS(m) + 1:m * HS(m) + N) = y(m * HS(m) + 1:m * HS(m) + N) + y_modified;
end

sound(y,Fs);

%function ppa: wrap any phase angle to phase over[-pi,pi]
function wrapped_phase = ppa(phase)
    wrapped_phase = mod(phase + pi, 2 * pi) - pi;
end
