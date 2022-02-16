classdef Quadriga < Channel
    %Quadriga. Wraps Quadriga
    
    properties
        layout
        scenario
        f_c
        fft_size
        sc_spacing
        theta_array
        dist_array
    end
    
    methods
        function obj = Quadriga(varargin)
            % QUADRIGA() Constructor for this class
            % Inputs:
            % p             Struct
            
            % Parse the inputs.
            vars = inputParser;
            validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x > 0);
            validStr = @(x) str(x) && isscalar(x);
            
            addParameter(vars, 'name', 'Quadriga', @(x) any(validatestring(x,{'Quadriga'})));
            addParameter(vars, 'required_domain', 'freq', @(x) any(validatestring(x,{'freq'})));
            addParameter(vars, 'required_fs', 122.88e6, validScalarPosNum);
            addParameter(vars, 'index', 1, validScalarPosNum);
            addParameter(vars, 'n_users', 1, validScalarPosNum);
            addParameter(vars, 'n_ants', 16, validScalarPosNum);
            addParameter(vars, 'f_c', 3.5e9, validScalarPosNum);
            addParameter(vars, 'fft_size', 4096, validScalarPosNum);
            addParameter(vars, 'sc_spacing', 15e3, validScalarPosNum);
            addParameter(vars, 'scenario', 'LOSonly', validStr)
            addParameter(vars, 'theta_array', 70);
            addParameter(vars, 'dist_array', 200);
            
            parse(vars, varargin{:});
            obj.save_inputs_to_obj(vars);
            
            s = obj.setup_sim_params(obj.f_c);
            a = obj.build_array();
            obj.build_layout(s, a);
            %bj.layout.visualize;
            obj.create_channel_matrix();
        end
        
        function report(obj)
            fprintf('Quadriga channel.\n')
        end
        
        function [Y_down] = subclass_use(obj, X1)
            % USE() Use channel
            
            % Go through the channel for each subcarrier
            
            [n_antennas, n_symbols, larger_fft] = size(X1);
            %[n_users, ~, ~] = size(X2);
            
            assert(n_antennas==obj.n_ants, 'Downlink data isnt correct dimensions for channel');
            %assert(n_users==obj.n_users, 'Uplink data isnt correct dimensions for channel');
            
            Y_down = zeros(obj.n_users, n_symbols, larger_fft);
            %Y_up = zeros(n_antennas, n_symbols, larger_fft);
            for i_subcarrier = 1:larger_fft
                HX1 = obj.H(:,:,i_subcarrier)*X1(:,:,i_subcarrier);
                %self_interference = obj.C(:,:,i_subcarrier) * X1(:,:,i_subcarrier);
                
                %HX2 = transpose(obj.H(:,:,i_subcarrier))*X2(:,:,i_subcarrier);
                Y_down(:,:,i_subcarrier) = HX1; %+ sqrt(N0)*obj.N(:,:,i_subcarrier);
                %Y_up(:,:,i_subcarrier) = HX2 + self_interference; %+ sqrt(N0)*obj.N(:,:,i_subcarrier);
            end
        end
        function Y_down = subclass_use_down(obj, X1)
            % USE() Use channel
            
            % Go through the channel for each subcarrier
            [~, n_symbols, larger_fft] = size(X1);
            
            Y_down = zeros(obj.n_ue, n_symbols, larger_fft);
            for i_subcarrier = 1:larger_fft
                HX1 = obj.H(:,:,i_subcarrier)*X1(:,:,i_subcarrier);
                Y_down(:,:,i_subcarrier) = HX1; %+ sqrt(N0)*obj.N(:,:,i_subcarrier);
            end
        end
        function Y_up = subclass_use_up(obj, X2)
            % USE() Use channel
            
            % Go through the channel for each subcarrier
            [~, n_symbols, larger_fft] = size(X2);
            Y_up = zeros(n_antobj.n_ant, n_symbols, larger_fft);
            for i_subcarrier = 1:larger_fft
                HX2 = transpose(obj.H(:,:,i_subcarrier))*X2(:,:,i_subcarrier);
                Y_up(:,:,i_subcarrier) = HX2; %+ sqrt(N0)*obj.N(:,:,i_subcarrier);
            end
        end
    end
    
    methods (Access = protected)
        function build_layout(obj, s, a)
            ant_height = 1.5;
            lambda = physconst('LightSpeed')/obj.f_c;
            length_of_array = obj.n_ants * 0.5 * lambda;
            obj.layout = qd_layout();
            obj.layout.simpar = s;
            obj.layout.tx_array = a;
            obj.layout.no_rx = obj.n_users;
            %  Get user locations from the users
            for i_user = 1:obj.n_users
                rx_distance = obj.dist_array(i_user);
                theta = obj.theta_array(i_user);
                x_loc = rx_distance * sin(theta*pi/180);
                y_loc = -rx_distance * cos(theta*pi/180);
                obj.layout.rx_position(:,i_user) = [x_loc, y_loc, ant_height];
                obj.layout.rx_array(i_user).center_frequency = obj.f_c;
                obj.layout.rx_track(i_user).no_snapshots = 1;
            end
            obj.layout.set_scenario(obj.scenario);
            obj.layout.tx_position = [0;length_of_array/2-(0.5 * lambda/2); ant_height];
        end
        
        function s = setup_sim_params(~, f_c)
            s = qd_simulation_parameters;
            s.center_frequency = f_c;
            s.use_random_initial_phase;
        end
        
        function a = build_array(obj)
            n_verticle_elements = 1;
            n_horizontal_elements = obj.n_ants;
            a = qd_arrayant('3gpp-3d', n_verticle_elements, n_horizontal_elements,...
                obj.f_c, 1, 0.5);
        end
        
        function create_channel_matrix(obj)
            % CREATE_CHANNEL_MATRIX() Creates channel matrix, H, between
            % all antenna elements
            % Inputs:
            % p             Struct
            channel = obj.layout.get_channels;
            f_up = obj.fft_size * obj.sc_spacing;
            
            % Get the channel for each user?
            % This is in RX Antenna, TX Antenna, Subcarrier, Time index;
            dummy = channel(1).fr(f_up,  obj.fft_size);
            [~,~,~, n_time_indexes] = size(dummy);
            H_all = zeros(obj.n_users, 1, obj.n_ants,  obj.fft_size, n_time_indexes);
            for i_user = 1:obj.n_users
                H_all(i_user, :,:,:,:) = channel(i_user).fr(f_up,  obj.fft_size);
            end
            
            % Reorganize the channel. Pick 1 time index
            % obj.n_ue_antennas, obj.n_enb_antennas
            obj.H = zeros(obj.n_users, obj.n_ants,  obj.fft_size);
            
            for i_user = 1:obj.n_users
                obj.H(i_user,:, :) = H_all(i_user, 1, :, :, 1);
                for i_sc = 1: obj.fft_size
                    obj.H(i_user,:, i_sc) = obj.H(i_user,:, i_sc)/norm(obj.H(i_user,:, i_sc)); % TODO. Fix the array so I don't need to do this
                end
            end
            obj.H = fftshift(obj.H, 3);
        end
    end
end