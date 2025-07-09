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
%      1: [distances thetas] for perfect data, r = 1 meter, angle rotating
%         from 0 to 2*pi.
%         Test options: [N] where N is repeats of data as rows (default is
%         10)
%
%      2: [velocities omegas] for perfect data, r = 1 meter, angle rotating
%         for 10 samples.
%         Test options: [N] where N is repeats of data as rows (default is
%         10)
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
if 1 <= nargin
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
for ith_test = 1:length(specific_test_cases)
    thisTestCase = specific_test_cases(ith_test,1);
    switch thisTestCase
        case 1
            %      1: [distances thetas] for perfect data, r = 1 meter.
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
        case 2
            %      2: [velocities omegas] for perfect data, r = 1 meter, angle rotating
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

        nexttile

        hold on;
        axis equal;
        grid on;
          
        plot(cellArrayOf_sampleData{ith_test}(:,1),cellArrayOf_sampleData{ith_test}(:,2),'.-','LineWidth',3,'DisplayName',sprintf('Test data %.0d',ith_test));
        title(sprintf('Test data %.0d',ith_test))

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
