classdef Signal < handle
    %Signal. Main entity that is passed between all modules. Assumes OFDM.
    %
    % Example
    %   my_signal = Signal(data, n_streams, domain, f_s)
    %
    %   If data is 'time' domain, then the dimensions should be
    %      (stream, sample)
    %   If data is 'freq' domain, then the dimensions should be
    %      (stream, bin, symbol)
    
    properties
        data
        n_streams % Number of parallel streams in this data.
        domain    % Domain of the data. 'time' or 'freq'
        fs        % Sample rate of the data in Hz.
        modulator % Holds the modulator to go back/forth from td to fd
        name
        figure_style
        rms_power
        papr
    end
    
    methods
        function obj = Signal(data, n_streams, domain, fs, modulator, name)
            %Signal Construct an instance of this class.
            if ~(strcmp(domain, 'freq') || strcmp(domain,'time'))
                error('This isnt a real domain. Choose freq or time');
            end
            obj.domain = domain;
            obj.fs = fs;
            obj.data = data;
            obj.n_streams = n_streams;
            obj.modulator = modulator;
            if nargin == 5
                name = '';
            end
            obj.name = name;
        end
        
        function match_this(obj, domain, fs)
            if nargin == 2
                fs = obj.fs;  % Assume current fs is good.
            end
            
            % Make sure the domain is correct
            obj.change_domain(domain);
            
            % Make sure the sample rate is correct
            % This is probably only for time domain data.
            if obj.fs ~= fs
                obj.change_fs(fs)
            end
        end
        
        function ber = calculate_bit_errors(obj, ref_bits)
            if nargin == 1
                ref_bits = obj.modulator.user_bits;
            end
            
            rx_bits = obj.modulator.demodulate(obj.data);
            n_errors = sum(ref_bits~=rx_bits);
            total_bits = length(rx_bits);
            ber = n_errors/total_bits;
        end
        
        function evm = calculate_evm(obj, ref_symbols)
            if nargin == 1
                ref_symbols = obj.modulator.original_fd;
            end
            
            error = ref_symbols - obj.data;
            mag = abs(ref_symbols).^2;
            evm = sum(error(:)) / sum(mag(:));
            
        end
        
        function freq_shift(obj,offset)
            % TODO.
        end
        
        function powers = calculate_current_rms_dbm(obj)
            % TODO.
        end
        
        function normalize_to_this_rms(obj, this_rms)
            % TODO.
        end
        
        function measure_channels(obj, channels)
            n_channels = numel(channels);
            all_rows = zeros(obj.n_streams, n_channels); 
            for i = 1:obj.n_streams
                all_channels = zeros(1, n_channels);
                for i_channel = -2:2
                    % Get the absolute powers for each channel in the stream
                    lower_limit = i_channel*carrier_bw - carrier_bw/2;
                    upper_limit = i_channel*carrier_bw + carrier_bw/2;
                    % add 3 to index this array from 1.
                    all_channels(i_channel+3) = obj.signal_array(i).measure_power([lower_limit, upper_limit]);
                end
                this_main = all_channels(3);
                
                % Convert to relative dBc
                this_row = all_channels - this_main;
                s
                % Put absolute main back in
                this_row(3) = this_main;
                all_rows(i, :) = this_row;
            end
        end
        
        function gain(obj, gain_amount)
            %gain. Apply a fixed gain or attenuation to each stream.
            % gain_amount: dB. Amplify if > 0. Attenuate if < 0
            
            current_power = obj.calculate_current_rms_dbm;
            for i_stream = 1:obj.n_streams
                obj.signal_array(i_stream).normalize_to_this_rms(current_power(i_stream) + gain_amount);
            end
        end
        
        function change_domain(obj, desired_domain)
            %change_domain. Method to parse the desired domain and call the
            %correct conversion method.
            if strcmp(desired_domain, 'bypass')
                return
            elseif strcmp(obj.domain, 'freq')  && strcmp(desired_domain, 'time')
                obj.data = obj.modulator.fd_to_td(obj.data);
                obj.domain = 'time';
            elseif strcmp(obj.domain, 'time')  && strcmp(desired_domain, 'freq')
                obj.data = obj.modulator.td_to_fd(obj.data);
                obj.domain = 'freq';
            end
        end
        
        function change_fs(obj, desired_fs)
            if strcmp(desired_fs, 'bypass')
                return
            elseif obj.fs < desired_fs
                warning('Upsampling from %d to %d', obj.fs, desired_fs);
                obj.upsample(desired_fs);
            elseif obj.fs > desired_fs
                warning('Downsampling from %d to %d', obj.fs, desired_fs);
                obj.downsample(desired_fs);
            else
                warning('Unexpected case where there shouldnt be up/down sampling');
            end
        end
        
        function plot_spectrogram(obj)
            if strcmp(obj.domain, 'freq')
                s_copy = obj.copy;
                s_copy.match_this('time');
                s_copy.plot_spectrogram;
                return;
            end
            figure()
            spectrogram(obj.data(1,:), 100,80,100,obj.fs, 'centered','yaxis')
        end
        
        function plot_psd(obj, fig_id)
            % PLOT_PSD() Plot all signals to input fig_id (98 if default).
            % Signals will be plotted in grid subplots.
            
            % Only works if Signal is time domain. If not, make copy,
            % convert, and plot that.
            if strcmp(obj.domain, 'freq')
                s_copy = obj.copy;
                s_copy.match_this('time');
                s_copy.plot_psd;
                return;
            end
            
            if nargin == 1
                fig_id = 98;
            end
            
            figure(fig_id)
            
            tile_set = [];
            for i_channel = 1:obj.n_streams
                [x, y, d] = obj.get_psd(i_channel);
                tile(i_channel) = subplot(floor(sqrt(obj.n_streams)),ceil(obj.n_streams/floor(sqrt(obj.n_streams))),i_channel);
                try
                    plot(x,y,obj.figure_style, 'DisplayName', obj.name);
                catch
                    try
                        plot(x,y);
                    catch
                        warning('Cant Plot PSD.');
                    end
                end
                hold on; grid on;
                title(sprintf('Rx Stream %d', i_channel));
                xlabel('Frequency (MHz)');
                ylabel(sprintf('PSD (dBm/%d kHz)', d/1e3));
                
                % Only show 1 legend for all plots
                if i_channel == obj.n_streams
                    legend show;
                    legend('Position',[0.91, 0.5, 0.1, 0.1])
                end
                
                tile_set = [tile_set tile(i_channel)];
            end
            
            linkaxes(tile_set,'xy')
            ylim([-120 0]);
        end
        
        function plot_constellation(obj, fig_id)
            if strcmp(obj.domain, 'time')
                s_copy = obj.copy;
                s_copy.match_this('freq');
                s_copy.plot_constellation;
                return;
            end
            if nargin == 1
                fig_id = 333;
            end
            figure(fig_id)
            
            for i_channel = 1:obj.n_streams
                this_fd_data = obj.data(i_channel,:,:);
                plot(this_fd_data(:), 'o');
            end
            hold on; grid on;
        end
        
        function plot_iq(obj, fig_id)
            if strcmp(obj.domain, 'freq')
                s_copy = obj.copy;
                s_copy.match_this('time');
                s_copy.plot_iq;
                return;
            end
            if nargin == 1
                fig_id = 101;
            end
            figure(fig_id)
            tile_set = [];
            
            N = length(obj.data);
            period = 1 / obj.fs;
            t = period * (1:N);
            
            for i_channel = 1:obj.n_streams
                y = obj.data(i_channel,: );
                tile(i_channel) = subplot(2,ceil(obj.n_streams/2),i_channel);
                plot(t,real(y)); hold on
                plot(t, imag(y));grid on;
                title(sprintf('Stream %d', i_channel));
                xlabel('Time (s)');
                ylabel('Amplitude');
                tile_set = [tile_set tile(i_channel)];
            end
            
            linkaxes(tile_set,'xy')
        end
        
        function S_copy = copy(obj)
            % Construct new mSignal based on obj
            S_copy = Signal(obj.data, obj.n_streams, obj.domain, ...
                obj.fs, obj.modulator, obj.name);
            
            % Copy all other optional/non constructor input properties
            S_copy.figure_style = obj.figure_style;
        end
    end
    
    methods (Access=private)
        function [X, Signal_PSD, density] = get_psd(obj, index)
            Nfft = 1024;
            Window = kaiser(1000, 9);
            
            X = (-1:2/Nfft:1-2/Nfft)*((obj.fs) / (2e6));
            Signal_PSD = 10 * log10(fftshift(pwelch(obj.data(index,:), Window))/sqrt(Nfft));
            
            density = obj.fs / Nfft; % In kHz
        end
        
        function upsample(obj, desired_fs)
            for i=1:obj.n_streams
                obj.signal_array(i).upsample(desired_fs);
            end
            obj.fs = desired_fs;
        end
        
        function downsample(obj, desired_fs)
            for i=1:obj.n_streams
                obj.signal_array(i).downsample(desired_fs);
            end
            obj.fs = desired_fs;
        end
    end
    
    methods (Static)
        function obj = make_ofdm(n_users, ofdm_settings)
            my_ofdm = OFDM(ofdm_settings, 'n_users', n_users);
            fd_data = my_ofdm.modulate();
            obj = Signal(fd_data, n_users, 'freq', my_ofdm.sampling_rate, my_ofdm);
        end
    end
end
