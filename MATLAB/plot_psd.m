function plot_psd(data, fs)
% plot_psd. 
% Plots the power spectral density of the given data.
%
% Inputs: 
%   data: complex vector of IQ samples
%   fs:   sample rate of the data in Hz
%
% Example:
% plot_psd(ofdm_data, 30.72e6);

Nfft = 1024;
Window = kaiser(1000, 9);

X = (-1:2/Nfft:1-2/Nfft)*((fs) / (2e6));
Signal_PSD = 10 * log10(fftshift(pwelch(data, Window))/sqrt(Nfft));
density = fs / Nfft; % In kHz
hold on;
plot(X, Signal_PSD, 'LineWidth', 0.5);
xlabel('Frequency (MHz)');
ylabel(sprintf('PSD (dBm/%d kHz)', density/1e3));
end