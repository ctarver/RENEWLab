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
        debug = 0
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
        function argos_sync(obj, reference)
           % Course frame estimation.
           power_threshold = 0.005;
           indexes_above = obj.data > power_threshold;
           pattern = zeros(1,100);
           pattern(end) = 1;  % Find a rising edge
           pos = strfind(indexes_above, pattern);
           extract_index = pos(1);
           extract_length = 3288;
           extracted_frame = obj.data(extract_index: extract_index+extract_length-1);
           
           figure;
           plot(indexes_above);
           hold on;
           plot(abs(obj.data));
           plot(diff(indexes_above));
            
           %angles = obj.data./reference.data;
           %cf0 = f_cfr_fitz(angles)
           %cyclosync()
        end
        
        function v = f_cfr_fitz(obj,z)
            %f_cfr_fitz Method for coarse CFO correction and estimation.
            
            L = length(z);
            N = floor(L/2);
            Fs = obj.fs;
            
            R=zeros(1,N);
            
            for m=1:N
                summa=0;
                for k=m+1:L
                    summa=summa+z(k)*conj(z(k-m));
                end
                R(m)=1/(L-m)*summa;
            end
            
            v=1/(pi*N*(N+1)/Fs)*sum(atan2(imag(R),real(R)));
        end        
        
        function ber = calculate_bit_errors(obj, ref_bits)
            obj.match_this('freq');
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
        
        function freq_shift(obj, offset)
            shift = exp(2*pi*1i*offset/obj.fs*(1:length(obj.data)));
            if(size(shift) == size(obj.data))
                obj.data = obj.data .* shift;
            else
                obj.data = obj.data .* shift.';
            end
        end
        
        function normalize_to_this_rms(obj, this_rms)
            scale_factor = obj.calculate_scale_factor(this_rms);
            obj.data = obj.data * scale_factor;
            obj.calculate_current_rms_dbm();
            if abs(this_rms-obj.rms_power) > 0.01
                error('RMS is wrong.');
            end
        end
        
        function calculate_current_rms_dbm(obj, index)
            if nargin == 1
                index = 1;
            end
            obj.rms_power = 10 * log10(norm(obj.data(index, :))^2/50/length(obj.data(index, :))) + 30;
        end
        
        function scale_factor = calculate_scale_factor(obj, desired_dbm_power, index)
            if nargin == 2
                index = 1;
            end 
            scale_factor = sqrt(50*length(obj.data(index,:))*10^((desired_dbm_power - 30) / 10)) / norm(obj.data(index,:));
        end
        
        function location = align_to(obj, S_reference_in, align_type, offset)
            %Aligns signal data to S_reference signal. Supports integer
            %sample alignment and subsample alignment with a few alignment
            %types:
            %   'all' - aligns all rx copies to usable correlation
            %       peak, returns all aligned copies array
            %   'best' - aligns rx data with highest correlation peak
            %   'average' - aligns all rx copies to usable correlation
            %       peaks and averages them to return one vector
            %   'none' - returns rx signal with no alignment or cropping
            %
            % offset: int. Optional, default = 0. alignment offset applied
            % to the crosscorrelation return value.
            
            S_reference = S_reference_in.copy();
            S_reference.match_this('time');
            
            raw = obj.data;
            data_ref = S_reference.data(1, :);  % Need only 1 user!
            
            if nargin == 2
                align_type = 'first';
                offset = 1;
            elseif nargin == 3
                offset = 0;
            end
            
            assert(obj.n_streams==1, 'This align function only supports 1 stream.')
            
            location = 0; % Just a default value.
            switch align_type
                case 'all'
                    % Do  normal sync stuff to get subsamp and phase.
                    out_aligned_and_averaged = obj.align_sample(obj.data, data_ref, false, offset);
                    
                    [~, phase_shift] = obj.align_phase(out_aligned_and_averaged, data_ref);
                    % Now do it but for multiple copies. Use calculations
                    % from above.
                    [r, lags] = xcorr(obj.data, data_ref);
                    peak_threshold = 0.9* max(abs(r));
                    [~, locs] = findpeaks(abs(r), lags, 'MinPeakHeight', peak_threshold);
                    aligned = out_up_raw(locs(2):end-1);  % Assumes offset = 0
                    X = [aligned, [0; aligned(1:end-1)]];
                    aligned = X*subsamp;
                    aligned = aligned*phase_shift/norm(phase_shift);
                    
                    obj.data = aligned;
                case 'first'
                    [out_aligned_and_averaged, location] = obj.align_sample(obj.data, data_ref, false, offset);
                    out_aligned_and_averaged = obj.align_phase(out_aligned_and_averaged, data_ref);
                    obj.data = out_aligned_and_averaged.';
                case 'average'
                    [out_aligned_and_averaged, location]  = obj.align_sample(obj.data, data_ref, true, offset);
                    out_synced = obj.align_phase(out_aligned_and_averaged, data_ref);
                    obj.data = out_synced;
                case 'none'
                    obj.data = raw;
                otherwise
                    error('Unknown align_type is Signal align.');
            end
        end
        
        function [out, coeffs] = align_phase(obj, align_this, to_this)
            %align_phase. Uses a LS fit to rotate x to fit y.
            % Inputs:
            %    - align_this. Vector of data that we will rotate
            %    - to_this.    Vector of data that we will align to.
            
            % Make sure both are column vectors
            align_this = reshape(align_this, [], 1);
            to_this = reshape(to_this, [], 1);
            coeffs = obj.perform_ls_estimation(align_this, to_this);
            out = align_this * coeffs / norm(coeffs);
        end
        
        function beta = perform_ls_estimation(obj, X, y)
            lambda = 0.001;
            beta = (X' * X + lambda * eye(size((X' * X)))) \ (X' * y);
        end
        
        
        
        function normalize_to_this_amp(obj, desired_amp)
            current_real_max = max(max(real(obj.data)));
            current_imag_max = max(max(imag(obj.data)));
            current_max = max(current_real_max, current_imag_max);
            
            obj.data = desired_amp * obj.data / current_max;
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
            
            %linkaxes(tile_set,'xy')
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
        
        function [out, location] = align_sample(obj, align_this, to_this, ...
                do_averaging, offset)
            % Uses crosscorrelation to align signals in time domain.
            % out = align_sample(align_this, to_this, do_averaging)
            % Inputs:
            %    - align_this.   Vector of data that we will align to the next input
            %    - to_this.      Vector of data that we will align to.
            %    - do_averaging. (Optional; default is false) Bool indicating that we
            %                    should averagemultiple transmissions.
            %    - offset        (Optional, default = 0) Int.
            %                    Offset to use in returning the aligned
            %                    version of the signal
            % Outputs:
            %    - out:          Output vector that is a modified version
            %                    of align_this
            % Averaging:
            %    Averaging is helpful in DPD and SIC learning. The
            %    do_averaging flag will cause it to average across each
            %    copy of the signal that was received. It automatically
            %    throws away the first and last version since they may be
            %    partial transmissions.
            
            % make sure that align_this is a column vector.
            align_this = reshape(align_this, [], 1);
            
            if nargin == 3
                offset = 1;
            end
            
            %normalize align_this and to_this to 1, such that absolute
            %thresholding can be used for finding correlation peaks
            if nnz(to_this) == 0
                norm_to_this = to_this;
                norm_align_this = align_this;
                transmit_data_ratio = 1;
            else
                transmit_data_ratio = nnz(to_this)/length(to_this);
                norm_align_this = align_this/max(align_this);
                norm_to_this = to_this/max(to_this);
            end
            
            location = nan; % default case.
            [r, lags] = xcorr(norm_align_this, norm_to_this);
            peak_threshold = max([0.9*max(abs(r)), 0]); % 1000*transmit_data_ratio]);
            [~, locs] = findpeaks(abs(r), lags, 'MinPeakHeight', peak_threshold);
            x = locs >= 0;
            locs = locs(x);
            if length(locs) == 1
                do_averaging = 0;
            elseif isempty(locs)
                disp('No alignment peak found. Maybe you have great isolation. Maybe there is a problem...');
                do_averaging = 0;
            end
            
            if do_averaging
                n_averages = max(length(locs)-2, 1); % Plan on throwing away the first and last peak.
                n_averages = min(n_averages, 100);
                out = zeros(length(to_this), 1);
                for i = 1:n_averages
                    out = out + align_this(locs(i+1)+offset:locs(i+1)+length(to_this)-1+offset);
                end
                out = out / n_averages;
                location = locs(2) + offset;
            else
                n_averages = 1;
                try
                    out = align_this(locs(2)+offset:locs(2)+length(to_this)-1+offset);
                    location = locs(2) + offset;
                catch
                    try
                        out = align_this(locs(1)+offset:locs(1)+length(to_this)-1+offset);
                        location = locs(1) + offset;
                    catch
                        out = zeros(length(to_this), 1);
                    end
                end
            end
            if obj.debug
                fprintf('\tNumber of averages: %d\n', n_averages);
                fprintf('\tAlignment Peaks: %d\n', locs);
                figure
                plot(lags, abs(r));
                hold on;
                plot(lags, (abs(r) > peak_threshold)*peak_threshold);
                xlabel('Sample Index')
                ylabel('Cross Correlation')
                title('Sample Alignment')
            end
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
