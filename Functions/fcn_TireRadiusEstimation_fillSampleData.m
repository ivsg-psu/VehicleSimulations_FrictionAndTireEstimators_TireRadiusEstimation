function cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData(varargin)
% fcn_TireRadiusEstimation_fillSampleData
% Produces dummy data to test TireRadiusEstimation functions. 
%
% FORMAT:
%
%       cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(fig_num))
%
% INPUTS:
%
%     (optional inputs)
%
%      specific_test_cases: an integer or vector of integers listing which
%      test cases will be loaded. For example, [1 2 3] will load test cases
%      1 through 3. Default is to load all test cases with default
%      test_options.
%
%      test_options: a cell array, one for each specific_test_cases, to set
%      optional inputs for each
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. 
%
% OUTPUTS:
%
%      cellArrayOf_sampleData: an cell array of data for testing. The
%      meaning of the data depends on the test cases. The test cases are as
%      follows:
%
%   CASE: xx - ideal cases, no noise
%      1: [distances thetas] for perfect data, r = 1 meter, angle rotating
%         from 0 to 2*pi.
%         Test options: [N] where N is repeats of data as rows (default is
%         10)
%
%      2: [velocities omegas] in units of [meters/sec], [radians/sec]
%         respectively,for perfect data, r = 1 meter, angle rotating
%         for 10 samples.
%         Test options: [N] where N is repeats of data as rows (default is
%         10)
%
%   CASE: 70xx - data collected for HSOV
%   7001: [distances thetas] in units of [meters], [radians]
%         respectively, for data collected on HSOV in Reber parking lot.
%         Radius is approximately 0.09 meters. 3 laps total.
%
%   7002: [velocities omegas] in units of [meters/sec], [radians/sec]
%         respectively, for data collected on HSOV in Reber parking lot.
%         Radius is approximately 0.09 meters. 3 laps total.
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%       See the script:
%       script_test_fcn_TireRadiusEstimation_fillSampleData for a full
%       test suite.
%
% This function was written on 2025_07_09 by S. Brennan
% Questions or comments? sbrennan@psu.edu

% Revision history:
% 2025_07_09 - S. Brennan
% -- wrote the code, using Laps repo data loading as starter

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 3; % The largest Number of argument inputs to the function
flag_max_speed = 0; % The default. This runs code with all error checking
if (nargin==MAX_NARGIN && isequal(varargin{end},-1))
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 0; % Flag to perform input checking
    flag_max_speed = 1;
else
    % Check to see if we are externally setting debug mode to be "on"
    flag_do_debug = 0; % % % % Flag to plot the results for debugging
    flag_check_inputs = 1; % Flag to perform input checking
    MATLABFLAG_TIRERADIUSESTIMATION_FLAG_CHECK_INPUTS = getenv("MATLABFLAG_TIRERADIUSESTIMATION_FLAG_CHECK_INPUTS");
    MATLABFLAG_TIRERADIUSESTIMATION_FLAG_DO_DEBUG = getenv("MATLABFLAG_TIRERADIUSESTIMATION_FLAG_DO_DEBUG");
    if ~isempty(MATLABFLAG_TIRERADIUSESTIMATION_FLAG_CHECK_INPUTS) && ~isempty(MATLABFLAG_TIRERADIUSESTIMATION_FLAG_DO_DEBUG)
        flag_do_debug = str2double(MATLABFLAG_TIRERADIUSESTIMATION_FLAG_DO_DEBUG);
        flag_check_inputs  = str2double(MATLABFLAG_TIRERADIUSESTIMATION_FLAG_CHECK_INPUTS);
    end
end

% flag_do_debug = 1;

if flag_do_debug % If debugging is on, print on entry/exit to the function
    st = dbstack; %#ok<*UNRCH>
    fprintf(1,'STARTING function: %s, in file: %s\n',st(1).name,st(1).file);
    debug_fig_num = 999978; %#ok<NASGU>
else
    debug_fig_num = []; %#ok<NASGU>
end

%% check input arguments?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____                   _
%  |_   _|                 | |
%    | |  _ __  _ __  _   _| |_ ___
%    | | | '_ \| '_ \| | | | __/ __|
%   _| |_| | | | |_) | |_| | |_\__ \
%  |_____|_| |_| .__/ \__,_|\__|___/
%              | |
%              |_|
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0==flag_max_speed
    if flag_check_inputs
        % Are there the right number of inputs?
        narginchk(0,MAX_NARGIN);

    end
end

% Does user specify specific_test_cases?
specific_test_cases = []; % Default is to use them all
if 1 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        specific_test_cases = temp; 

        if flag_check_inputs
            % Check the specific_test_cases variable
            fcn_DebugTools_checkInputsToFunctions(specific_test_cases, '1column_of_numbers');
        end
    end
end

% Does user specify test_options?
test_options = [];
if 2 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        test_options = temp; 
        if flag_check_inputs
            % Check the test_options variable
            assert(iscell(test_options));
        end
    end
end

% Does user want to show the plots?
flag_do_plots = 0; % Default is to NOT show plots
if (0==flag_max_speed) && (MAX_NARGIN == nargin) 
    temp = varargin{end};
    if ~isempty(temp) % Did the user NOT give an empty figure number?
        fig_num = temp;
        figure(fig_num);
        flag_do_plots = 1;
    end
end

%% Main code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   __  __       _
%  |  \/  |     (_)
%  | \  / | __ _ _ _ __
%  | |\/| |/ _` | | '_ \
%  | |  | | (_| | | | | |
%  |_|  |_|\__,_|_|_| |_|
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(specific_test_cases)
    NtestCases = 2;
    specific_test_cases = (1:NtestCases)';
else
    NtestCases = length(specific_test_cases);
end


cellArrayOf_sampleData = cell(NtestCases,1);
dataTypes = zeros(NtestCases,1);
plotTitles = cell(NtestCases,1);
for ith_test = 1:length(specific_test_cases)
    thisTestCase = specific_test_cases(ith_test,1);
    switch thisTestCase
        case 1
            %      1: [distances thetas] for perfect data, r = 1 meter, angle rotating
            %         from 0 to 2*pi.
            %         Test options: [N] where N is repeats of data as rows (default is
            %         10)
            if isempty(test_options)
                Npoints = 10;
            else
                Npoints = test_options{ith_test};
            end
            thetas = linspace(0,2*pi,Npoints)';
            distances = thetas*1.0;
            cellArrayOf_sampleData{ith_test,1} = [distances thetas]; 
            dataTypes(ith_test,1) = 1;
            plotTitles{ith_test,1} = sprintf('Perfect position data, N=%.0d samples',Npoints);
        case 2
            %      2: [velocities omegas] in units of [meters/sec], [radians/sec]
            %         respectively,for perfect data, r = 1 meter, angle rotating
            %         for 10 samples.
            %         Test options: [N] where N is repeats of data as rows (default is
            %         10)
            if isempty(test_options)
                Npoints = 10;
            else
                Npoints = test_options{ith_test};
            end
            omegas = ones(Npoints,1);
            velocities = omegas*1.0;
            cellArrayOf_sampleData{ith_test,1} = [velocities omegas]; 
            dataTypes(ith_test,1) = 2;
            plotTitles{ith_test,1} = sprintf('Perfect velocity data, N=%.0d samples',Npoints);

        case 7001
            % CASE: 70xx - data collected for HSOV
            %      7001: [distances thetas] in [meters, radians] for data
            %      collected on HSOV in Reber parking lot. Radius is
            %      approximately 0.09 meters. 3 laps total.
            load('wheel_radius_HSOV_data_3lapsAtReber.mat', 'all_data_trimmed_to_laps', 'reference_LLA','trimmingBoundaries');
            Nlaps = length(all_data_trimmed_to_laps);
            cellArrayOf_distances = cell(Nlaps,1);
            cellArrayOf_thetas    = cell(Nlaps,1);
            previousDistance = 0;
            previousTheta = [0 0 0 0];
            for ith_lap = 1:Nlaps
                this_distance = all_data_trimmed_to_laps(ith_lap).gpsStation - all_data_trimmed_to_laps(ith_lap).gpsStation(1) + previousDistance;
                cellArrayOf_distances{ith_lap,1} = this_distance;
                previousDistance = this_distance(end,1);

                rawThetas = all_data_trimmed_to_laps(ith_lap).encoderRadsAtGPSTimes;
                Nraw = length(rawThetas(:,1));

                % Remove initial offsets in the lap, and force the angle to
                % start at same point of previous lap
                this_theta = rawThetas + ones(Nraw,1)*(previousTheta - rawThetas(1,:));
                cellArrayOf_thetas{ith_lap,1} = this_theta;
                previousTheta = this_theta(end,:);

            end
            distances_thetas = [cell2mat(cellArrayOf_distances) cell2mat(cellArrayOf_thetas)];

            % Remove repeats
            distances_thetas_unique = unique(distances_thetas,'rows','stable');
            cellArrayOf_sampleData{ith_test,1} = distances_thetas_unique;
            dataTypes(ith_test,1) = 1;
            Npoints = length(cellArrayOf_sampleData{ith_test,1}(:,1));
            plotTitles{ith_test,1} = sprintf('HSOV data, Reber 3 laps, N=%.0d',Npoints);


        case 7002
            %   7002: [velocities omegas] in units of [meters/sec], [radians/sec]
            %         respectively, for data collected on HSOV in Reber parking lot.
            %         Radius is approximately 0.09 meters. 3 laps total.
            
            load('wheel_radius_HSOV_data_3lapsAtReber.mat', 'all_data_trimmed_to_laps', 'reference_LLA','trimmingBoundaries');
            Nlaps = length(all_data_trimmed_to_laps);
            cellArrayOf_velocities = cell(Nlaps,1);
            cellArrayOf_omegas    = cell(Nlaps,1);

            for ith_lap = 1:Nlaps
                this_velocity = all_data_trimmed_to_laps(ith_lap).GroundVelocity;
                cellArrayOf_velocities{ith_lap,1} = this_velocity;

                this_omega = all_data_trimmed_to_laps(ith_lap).encoderRadPerSecAtGPSTime;
                cellArrayOf_omegas{ith_lap,1} = this_omega;
            end
            velocities_omegas = [cell2mat(cellArrayOf_velocities) cell2mat(cellArrayOf_omegas)];

            cellArrayOf_sampleData{ith_test,1} = velocities_omegas;
            dataTypes(ith_test,1) = 2;
            Npoints = length(cellArrayOf_sampleData{ith_test,1}(:,1));
            plotTitles{ith_test,1} = sprintf('HSOV data, Reber 3 laps, N=%.0d',Npoints);

        otherwise
            error('Unknown test given: %.0d',thisTestCase);
    end % Ends switch

end % Ends for loop through cases

%% Any debugging?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   _____       _
%  |  __ \     | |
%  | |  | | ___| |__  _   _  __ _
%  | |  | |/ _ \ '_ \| | | |/ _` |
%  | |__| |  __/ |_) | |_| | (_| |
%  |_____/ \___|_.__/ \__,_|\__, |
%                            __/ |
%                           |___/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag_do_plots
    % Prep a figure
    figure(fig_num);
    clf;

    tiledlayout('flow');

    % Show results
    for ith_test = 1:length(cellArrayOf_sampleData)
        thisData = cellArrayOf_sampleData{ith_test};
        caseNumber = specific_test_cases(ith_test);
        Nencoders  = length(thisData(1,:))-1;

        nexttile

        hold on;
        grid on;
          
        for ith_encoder = 1:Nencoders
            plot(thisData(:,1),thisData(:,ith_encoder+1),'.-','LineWidth',3,'DisplayName',sprintf('Case %.0d, Encoder %.0d',caseNumber,ith_encoder));
        end
        title(plotTitles{ith_test});

        if 1==dataTypes(ith_test,1)
            ylabel('Thetas [rad]');
            xlabel('Distances [m]');
        elseif 2==dataTypes(ith_test,1)
            ylabel('Omegas [rad/sec]');
            xlabel('Velocities [m/s]');
        end
        
    end
    sgtitle('All test case data');
end

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end
end
