classdef RealChannel < Module
    %REALCHANNEL. This channel subclass is meant to be used with the FPGA
    %platform for learning and loading channels.
    %
    %
    % Chance Tarver
    % February 2021
    
    properties
        H           % Channel matrix. Freq domain. (user_rx, tx_ant, subcarrier)
        n_ants      % int. number of TX antennas in the array.
        n_users     % int. Number of UE RXs
        theta       % int. array of angles to users. In learn mode, there should only be 1 angle.
        distance    % int. array of distances to the users. In learn mode, there should only be 1 angle.
        n_scs       % int. Number of subcarriers to learn over.
        fft_size    % int. Size of the fft to use. Should be power of 3 bigger than n_scs.
        sc_spacing  % int. spacing between subcarriers in the learn signal.
        f_center    % int. Center frequency in Hz for test. Used for setting correct string when loading and saving channel.
        pilots      % complex matrix. saved when pilots are generated. compared against demod data to learn channel.
        p_sig       % msig. All the pilots wrapped in a msig. % TODO. could remove above.
        one_shot    % bool. Flag to make pilots all in 1 msig.
    end
    
    methods
        function obj = RealChannel(varargin)
            % Parse the inputs.
            vars = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validBool = @(x) islogical(x);
            addParameter(vars, 'name', 'RealChannel', @(x) any(validatestring(x,{'RealChannel'})));
            addParameter(vars, 'required_domain', 'freq', @(x) any(validatestring(x,{'freq'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'n_ants', 8, validScalarPosNum);
            addParameter(vars, 'n_users', 1, validScalarPosNum);
            addParameter(vars, 'sc_spacing', 15e3, validScalarPosNum);
            addParameter(vars, 'f_center', 3.5e9, validScalarPosNum);
            addParameter(vars, 'fft_size', 4096, validScalarPosNum);
            addParameter(vars, 'n_scs', 1024, validScalarPosNum);
            addParameter(vars, 'theta', 90, validScalarPosNum);
            addParameter(vars, 'distance', 200, validScalarPosNum);
            addParameter(vars, 'one_shot', true, validBool);
            parse(vars, varargin{:});
            
            % Save inputs to obj
            fields = fieldnames(vars.Results);
            for i = 1:numel(fields)
                obj.(fields{i}) = vars.Results.(fields{i});
            end
            
            %if strcmp(p.channel.mode, 'load')
            %    obj.load();
            %end
        end
        
        
        function normalize(obj)
            %the_max = max(max(abs(obj.H)));
            %obj.H = obj.H/the_max;
        end
        
        
        
        function learn(obj, pilots_in, ue_msig_in, i)
            %learn. Use the array class passed in to learn the channel.
            % Create a signal.
            pilots = pilots_in.copy();
            pilots.match_this('freq');
            ue_msig = ue_msig_in.copy();
            %ue_msig.match_this('time', obj.required_fs);
            
            if obj.one_shot
                obj.learn_one_shot(pilots, ue_msig);
            else
                obj.learn_this_tx_channel(pilots, ue_msig, i);
            end
            
            % Save the channel in the archive.
            obj.save;
        end
        
        function S = create_pilots(obj, ofdm_settings, i_tx)
            if obj.one_shot
                S = obj.make_one_shot_pilot(ofdm_settings);
            else
                S = obj.make_single_pilot_symbol(ofdm_settings, i_tx);
            end
        end
        
        function save(obj)
            channel = obj.H;
            %str = obj.format_file_name(1);
            %path_to_lib = obj.path_of_library;
            %name = strcat(path_to_lib, str);
            name = 'matlab_channel';
            save(name, 'channel') % TODO, get it in the correct path.
        end
        
        function load(obj)
            % load. Load channel vectors from memory and create 1 channel matrix.
            n_ues = length(obj.theta);
            path_to_lib = obj.path_of_library();
            obj.H = zeros(n_ues, obj.n_tx, obj.n_scs);
            for i_ue = 1:n_ues
                str = obj.format_file_name(i_ue);
                name = strcat(path_to_lib, str);
                this_h = load(name);
                obj.H(i_ue, :, :) = this_h.channel;
            end
        end
        
        function report(obj)
            
        end
    end
    
    methods (Access = protected)
        function learn_this_tx_channel(obj, pilots, ue_msig, i)
            % Learn the channel in the case that we do
            % Maybe need to do an align and extract to the data from user 1.
            % extract just the pilots for tx 1.
            %ue_msig = mSignal(ue_msig_full.signal_array(obj.i_ue).data, 1, pilots.domain, pilots.fs, pilots.mod_settings);
            %ue_msig.match_this('time');
            
            %s1 = mSignal(pilots.signal_array(i).data, 1, pilots.domain, pilots.fs, pilots.mod_settings);
            %s1.match_this('time', obj.required_fs);
            %ue_msig.align_to_group(s1);
            
            %% Demod
            ue_msig.match_this('freq');
            Y = ue_msig.data;
            
            % Rearange into grid.
            Y = squeeze(Y(1,:,:));
            %Y = Y.'; % Subcarriers should be columns
            
            [n_subcarriers, n_symbols] = size(Y);
            
            if isempty(obj.H)
                obj.H = zeros(1, obj.n_ants, obj.n_scs); % 1st index is user.
            end
            % For each symbols compare to perfect pilot that we
            % transmitted.
            
            % For each subcarrier, compute
            this_x = squeeze(pilots.data(1, :, :)); % Grab pilot for this TX
            data_scs = abs(this_x(1,:))>0;
            %this_x = this_x.';
            this_h = Y(:,data_scs)./this_x(:,data_scs);
            figure(22);
            plot(angle(fftshift(this_h.')));  grid on;
            change_in_phase = mean(angle(this_h(6,:)) - angle(this_h(1,:)));
            change_in_time = (6*548)/7.68e6;
            radians_per_sec = change_in_phase / change_in_time;
            cfo = radians_per_sec / (2*pi)
            obj.H(1, i, data_scs) = mean(this_h);
        end
        
        function S = make_single_pilot_symbol(obj, ofdm_settings, i_tx)
            %Will create one pilot symbol. This is intended to be used to
            %do one-at-a-time learning across all the TXs.
            
            % We will create 1 symbol.
            %Will create a pilot sequence to fit within the number of
            %samples specified.
            
            % We will create 1 symbol per TX. They will not overlap in
            % time.
            ofdm_settings.n_users = 1;
            ofdm_settings.n_symbols = 6;
            ofdm_settings.n_scs = ofdm_settings.fft_size;
            S = Signal.make_ofdm(1, ofdm_settings);
            % Delete data for other antennas..
            %complement = setxor(i_tx, 1:obj.n_ants);
            %S.data(complement, :, :)  = zeros(size(S.data(complement, :, :)));
            obj.p_sig = S;
        end
        
        function learn_one_shot(obj, pilots, ue_msig)
            % Learn the channel in the case that we send one-shot pilots.
            % Maybe need to do an align and extract to the data from user 1.
            % extract just the pilots for tx 1.
            
            %% Align
            %s1 = Signal(pilots.data(1,:,:), 1, pilots.domain, pilots.fs, pilots.modulator);
            %s1.match_this('time', obj.required_fs);
            %ue_msig.align_to_group(s1);
            
            %% Demod
            ue_msig.match_this('freq');
            Y = ue_msig.data;
            
            % Rearange into grid.
            Y = squeeze(Y(1,:,:));  % ASSUMING 1 USER. TODO
            Y = Y.'; % Subcarriers should be columns
            
            [n_subcarriers, n_symbols] = size(Y);
            
            obj.H = zeros(1, obj.n_ants, obj.n_scs); % 1st index is user.
            % For each symbols compare to perfect pilot that we
            % transmitted.
            for i = 1:n_symbols
                for i_ant = 1:obj.n_ants
                    this_x = squeeze(pilots.data(i_ant, i, :));
                    data_scs = abs(this_x)>0.00001;
                    % For each subcarrier, compute
                    this_x = squeeze(pilots.data(i_ant, i, :));
                    this_h = Y(data_scs, i_ant)./this_x(data_scs);
                    data_ind = find(data_scs==1);
                    
                    % Interpolate
                    full_h = interp1(data_ind, angle(this_h),1:obj.n_scs,'linear','extrap');
                    
                    
                    obj.H(1, i_ant, :) =  exp((full_h') * 1i);
                end
            end
            
            % Interpolate.
            
        end
        
        function S = make_one_shot_pilot(obj, ofdm_settings)
            %Will create a pilot sequence to fit within the number of
            %samples specified.
            
            % We will create 1 symbol per TX. They will not overlap in
            % time.
            ofdm_settings.n_users = obj.n_ants;
            ofdm_settings.n_scs = ofdm_settings.fft_size;
            S = Signal.make_ofdm(obj.n_ants, ofdm_settings);
            
            % Delete data so ortho in frequency.
            % Loop over the subcarriers.
            for i_sc = 1:ofdm_settings.fft_size
                active_antenna = mod(i_sc, obj.n_ants)+1;
                for i_tx = 1:obj.n_ants
                    if i_tx == active_antenna
                        % Do Nothing
                    else
                        S.data(i_tx, :, i_sc)  = zeros(size(S.data(i_tx, :, i_sc) ));
                    end
                end
            end
            
            obj.p_sig = S;
        end
        
        function str = format_file_name(obj, i_ue)
            % i_ue: int for which ue index you are trying to access. Will
            % use the object.
            this_theta = obj.theta(i_ue); this_dist = obj.distance(i_ue);
            f = obj.f_center; n_sc = obj.n_scs;
            str = sprintf('%d_degrees_%d_meters_%d_hz_%d_scs.mat', this_theta, this_dist, f, n_sc);
        end
        
    end
    methods (Static)
        function path = path_of_library()
            root_path = Helpers.get_root_path;
            path = strcat(root_path,'Modules\Channels\channel_library\');
        end
    end
end