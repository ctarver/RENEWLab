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
        
        function learn(obj, pilots, ue_msig_in, i)
            %learn. Use the array class passed in to learn the channel.
            % Create a signal.
            ue_msig = ue_msig_in.copy();
            ue_msig.match_this('time', obj.required_fs);
            
            if obj.one_shot
                obj.learn_one_shot(pilots, ue_msig);
            else
                obj.learn_this_tx_channel(pilots, ue_msig, i);
            end
            
            % Save the channel in the archive.
            obj.save;
        end
        
        function S = create_pilots(obj, i_tx)
            if obj.one_shot
                S = obj.make_one_shot_pilot();
            else
                S = obj.make_single_pilot_symbol(i_tx);
            end
        end
        
        function save(obj)
            channel = obj.H;
            str = obj.format_file_name(1);
            path_to_lib = obj.path_of_library;
            name = strcat(path_to_lib, str);
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
        function learn_this_tx_channel(obj, pilots, ue_msig_full, i)
            % Learn the channel in the case that we do
            % Maybe need to do an align and extract to the data from user 1.
            % extract just the pilots for tx 1.
            ue_msig = mSignal(ue_msig_full.signal_array(obj.i_ue).data, 1, pilots.domain, pilots.fs, pilots.mod_settings);
            ue_msig.match_this('time');
            
            s1 = mSignal(pilots.signal_array(i).data, 1, pilots.domain, pilots.fs, pilots.mod_settings);
            s1.match_this('time', obj.required_fs);
            ue_msig.align_to_group(s1);
            
            %% Demod
            ue_msig.match_this('freq');
            Y = ue_msig.extract_data;
            
            % Rearange into grid.
            Y = squeeze(Y(1,:,:));
            Y = Y.'; % Subcarriers should be columns
            
            [n_subcarriers, n_symbols] = size(Y);
            
            if isempty(obj.H)
                obj.H = zeros(1, obj.n_tx, obj.n_scs); % 1st index is user.
            end
            % For each symbols compare to perfect pilot that we
            % transmitted.
            
            % For each subcarrier, compute
            this_x = squeeze(obj.pilots(i, 1, :)); % Grab pilot for this TX
            this_x = this_x.';
            this_h = Y./this_x;
            obj.H(1, i, :) = this_h;
        end
        
        function S = make_single_pilot_symbol(obj, i)
            %Will create one pilot symbol. This is intended to be used to
            %do one-at-a-time learning across all the TXs.
            
            % We will create 1 symbol.
            n_symbols = 1;
            obj.pilots = zeros(obj.n_tx, n_symbols, obj.n_scs);
            symbol_alphabet = [-1-1i;-1+1i;1+1i;1-1i];
            rng(0);
            obj.pilots(i, 1, : ) = symbol_alphabet(randi(4,obj.n_scs,1));
            
            ofdm_settings.name = 'ofdm';
            ofdm_settings.n_symbols = n_symbols;
            ofdm_settings.n_scs = obj.n_scs;
            ofdm_settings.sc_spacing = obj.sc_spacing;
            ofdm_settings.fft_size = obj.fft_size;
            ofdm_settings.cp_length = 0; % TODO. Do we need this for right now?
            ofdm_settings.window_length = 0;
            ofdm_settings.rrc_taps = 1;
            ofdm_settings.clip_index = 0;
            fs =  obj.sc_spacing * obj.fft_size;
            assert(fs==obj.required_fs, 'Calculated fs doesnt match specified fs');
            S = mSignal(obj.pilots, obj.n_tx, obj.required_domain, obj.required_fs, ofdm_settings);
            obj.p_sig = S;
        end
        
        function learn_one_shot(obj, pilots, ue_msig)
            % Learn the channel in the case that we sed one-shot pilots.
            % Maybe need to do an align and extract to the data from user 1.
            % extract just the pilots for tx 1.
            s1 = mSignal(pilots.signal_array(1).data, 1, pilots.domain, pilots.fs, pilots.mod_settings);
            s1.match_this('time', obj.required_fs);
            ue_msig.align_to_group(s1);
            
            %% Demod
            ue_msig.match_this('freq');
            Y = ue_msig.extract_data;
            
            % Rearange into grid.
            Y = squeeze(Y(1,:,:));
            Y = Y.'; % Subcarriers should be columns
            
            [n_subcarriers, n_symbols] = size(Y);
            
            obj.H = zeros(1, obj.n_tx, obj.n_scs); % 1st index is user.
            % For each symbols compare to perfect pilot that we
            % transmitted.
            for i = 1:n_symbols
                % For each subcarrier, compute
                this_x = squeeze(obj.pilots(i,i,:));
                this_h = Y(:, i)./this_x;
                obj.H(1,i,:) = this_h;
            end
        end
        
        function S = make_one_shot_pilot(obj)
            %Will create a pilot sequence to fit within the number of
            %samples specified.
            
            % We will create 1 symbol per TX. They will not overlap in
            % time.
            n_symbols = obj.n_tx;
            obj.pilots = zeros(obj.n_tx, n_symbols, obj.n_scs);
            symbol_alphabet = [-1-1i;-1+1i;1+1i;1-1i];
            rng(0);
            for i_tx = 1:obj.n_tx
                % The pilots can be anything as long as we keep track of
                % them. We will just do a QPSK on each subcarrier
                obj.pilots(i_tx, i_tx, : ) = symbol_alphabet(randi(4,obj.n_scs,1));
            end
            
            ofdm_settings.name = 'ofdm';
            ofdm_settings.n_symbols = n_symbols;
            ofdm_settings.n_scs = obj.n_scs;
            ofdm_settings.sc_spacing = obj.sc_spacing;
            ofdm_settings.fft_size = obj.fft_size;
            ofdm_settings.cp_length = 0; % TODO. Do we need this for right now?
            ofdm_settings.window_length = 0;
            ofdm_settings.rrc_taps = 1;
            ofdm_settings.clip_index = 0;
            fs =  obj.sc_spacing * obj.fft_size;
            assert(fs==obj.required_fs,'Calculated fs doesnt match specified fs');
            S = Signal(obj.pilots, obj.n_tx, obj.required_domain, obj.required_fs, ofdm_settings);
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