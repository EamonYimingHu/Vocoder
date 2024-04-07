%-------------------------------------------------
% Audio Signal Processing with Matlab
%
% Code Spectrogram for TS and PV
%
% Yiming HU 11/11/23
%-------------------------------------------------


function ShowSpectrogram(x, Fs, N, O)
    L = length(x);
    HA = floor((1 - O) * N);
    
    % Generate a Hann window of length N manualLy
    win = 0.5 * (1 - cos(2 * pi * (0:N - 1).' / N));    
    win(1) = 0;                                         % Make sure the first value is zero to make it periodic
    NFFT = 2 ^ nextpow2(N);                             % Work out NFFT for power-of-two FFT
    NF = floor((L - N) / HA) + 1;

    % Pre-allocate the STFT matrix for efficiency
    stft_matrix = zeros(NFFT / 2 + 1, NF);              % Initilise matrix for STFT results


    for m = 1:NF - 1
        xm = win .* x(m * HA + 1:m * HA + N);           %extract mth analysis frame rom vector x and apply Hann win
        X = fft(xm, NFFT);                              % FFT on every frame with length NFFT
        %For simplicity, use fft's in-built functionality, i.e if length parameter is greater than the vector fft will automatically pad the vector up with zero
        stft_matrix(:, m) = X(1:NFFT / 2 + 1);          %store the DFT of each overlapping frame into one column
    end

    % Convert to magnitude in decibels
    stft_matrix_db = 20 * log10(abs(stft_matrix));
    
    % Set axis
    f = (0:NFFT / 2) * Fs / NFFT;                       % frequency vector from zero to the Nyquist frequency
    t = (0:length(x) - 1) * (1 / Fs);                   % Time vector 

    % Plot the spectrogram
    figure;
    imagesc(t, f, stft_matrix_db);                      % plot as as a two-dimensional image
    axis xy;                                            % frequency ascend going up, time ascend going right
    colormap('jet');                                    % color map
    colorbar;                                           % show color scale
    max_db = max(stft_matrix_db(:));
    caxis([max_db - 60,max_db]);                        % set color axis as 60 dB range 
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title_str = sprintf('Spectrogram - Frame Length: %.3f seconds | Overlap: %.1f%%', N/Fs, O * 100);
    title(title_str);
end
