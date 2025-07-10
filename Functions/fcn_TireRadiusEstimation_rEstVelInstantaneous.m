function rEff = fcn_TireRadiusEstimation_rEstVelInstantaneous(velocities_omegas, varargin)
% fcn_TireRadiusEstimation_rEstVelInstantaneous
% Produces instantaneous estimates of effective wheel radius based on
% velocity inputs.
%
% FORMAT:
%
%       rEff = fcn_TireRadiusEstimation_rEstVelInstantaneous(velocities_omegas, (plotXvalues), (plotXlabelString), (fig_num))
%
% INPUTS:
%      velocities_omegas: an NxM matrix of [distances thetas] in units of
%      [meters/sec], [radians/sec], where distances is always a Nx1 and
%      thetas is Nx(M-1) to allow multiple encoder readings at same time.
%
%      (optional inputs)
%
%      plotXvalues: an Nx1 array of values associated with each
%      velocities_omegas row. This is useful to indicate time samples for
%      each data, or station location for each data. Default, if left
%      empty, is to show indices, e.g. (1:N)'. Only used for plotting.
%
%      plotXlabelString: a string for labeling the x-axis according to the
%      plotXvalues. This can be, for example, 'Time [sec]' or 'Station
%      [m]'. The default, if left empty, is 'Indices (unitless)'. Only used
%      for plotting.
%
%      fig_num: a figure number to plot results. If set to -1, skips any
%      input checking or debugging, no figures will be generated, and sets
%      up code to maximize speed. 
%
% OUTPUTS:
%
%      rEff: an Nx(M-1) vector of estimated effective wheel radii
%
% DEPENDENCIES:
%
%      fcn_DebugTools_checkInputsToFunctions
%
% EXAMPLES:
%
%      See the script:
%      script_test_fcn_TireRadiusEstimation_rEstVelInstantaneous for a
%      full test suite.
%
% This function was written on 2025_04_16 by G. Gayoso for HSOV library,
% with modifications for the TireRadiusEstimation library by S. Brennan.
% Questions or comments? sbrennan@psu.edu

% Revision history:    
% 2025_04_16  - G. Gayoso
% -- wrote the code originally for HSOV library
% 2025_07_10  - S. Brennan
% -- modified the code in prep for export to TireRadiusEstimation library

% TO-DO
% (none)

%% Debugging and Input checks

% Check if flag_max_speed set. This occurs if the fig_num variable input
% argument (varargin) is given a number of -1, which is not a valid figure
% number.
MAX_NARGIN = 4; % The largest Number of argument inputs to the function
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
        narginchk(1,MAX_NARGIN);
        
        % Check the velocities_omegas input
        fcn_DebugTools_checkInputsToFunctions(velocities_omegas, '2orMorecolumn_of_numbers');
    end
end

% Does user specify plotXvalues?
plotXvalues = []; % Default is to use indices
if 2 <= nargin
    temp = varargin{1};
    if ~isempty(temp)
        plotXvalues = temp; 

        if flag_check_inputs
            Ndata = length(velocities_omegas(:,1));
            % Check the plotXvalues input
            fcn_DebugTools_checkInputsToFunctions(plotXvalues, '1column_of_numbers',[Ndata Ndata]);
        end
    end
end

% Does user specify plotXlabelString?
plotXlabelString = [];
if 3 <= nargin
    temp = varargin{2};
    if ~isempty(temp)
        plotXlabelString = temp; 
        if flag_check_inputs
            % Check the test_options variable
            assert(isstring(plotXlabelString)||ischar(plotXlabelString));
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

rEff = velocities_omegas(:,1)./velocities_omegas(:,2:end);
rEff(isinf(rEff)) = nan;


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
    
    Ndata = length(velocities_omegas(:,1));
    if isempty(plotXvalues)
        plotXvalues = (1:Ndata)';
    end
    
    if isempty(plotXlabelString)
        plotXlabelString = 'Indices (unitless)';
    end

    % Prep the figure for plotting
    temp_h = figure(fig_num);
    flag_rescale_axis = 0;
    if isempty(get(temp_h,'Children'))
        flag_rescale_axis = 1;
    end

    % Is this 2D or 3D?
    dimension_of_points = 2;

    % Find size of plotting domain
    allPoints = [plotXvalues rEff(:,1)];
    max_plotValues = max(allPoints);
    min_plotValues = min(allPoints);
    sizePlot = max(max_plotValues) - min(min_plotValues);
    nudge = sizePlot*0.006; %#ok<NASGU>


    % Find size of plotting domain
    if flag_rescale_axis
        percent_larger = 0.1;
        axis_range = max_plotValues - min_plotValues;
        if (0==axis_range(1,1))
            axis_range(1,1) = 1/percent_larger;
        end
        if (0==axis_range(1,2))
            axis_range(1,2) = 1/percent_larger;
        end
        if dimension_of_points==3 && (0==axis_range(1,3))
            axis_range(1,3) = 1/percent_larger;
        end


        % Force the axis to be equal?
        if 1==0
            min_valuesInPlot = min(min_plotValues);
            max_valuesInPlot = max(max_plotValues);
        else
            min_valuesInPlot = min_plotValues;
            max_valuesInPlot = max_plotValues;

        end

        % Stretch the axes
        stretched_min_vertexValues = min_valuesInPlot - percent_larger.*axis_range;
        stretched_max_vertexValues = max_valuesInPlot + percent_larger.*axis_range;
        axesTogether = [stretched_min_vertexValues; stretched_max_vertexValues];
        newAxis = reshape(axesTogether, 1, []);
        axis(newAxis);

    end
    goodAxis = axis;

    % Check to see if hold is already on. If it is not, set a flag to turn it
    % off after this function is over so it doesn't affect future plotting
    flag_shut_hold_off = 0;
    if ~ishold
        flag_shut_hold_off = 1;
        hold on
    end

    Nencoders  = length(rEff(1,:));
    for ith_encoder = 1:Nencoders
        plot(plotXvalues(:,1),rEff(:,ith_encoder),'.-','LineWidth',3,'DisplayName',sprintf('Encoder %.0d',ith_encoder));
    end
    
    ylabel('rEff [m]');
    xlabel(plotXlabelString);

    axis(goodAxis);

    legend('Interpreter','none');

    % Shut the hold off?
    if flag_shut_hold_off
        hold off;
    end

end % Ends if statement to check if plotting should happen

if flag_do_debug
    fprintf(1,'ENDING function: %s, in file: %s\n\n',st(1).name,st(1).file);
end

end % Ends main function

%% Functions follow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ______                _   _
%  |  ____|              | | (_)
%  | |__ _   _ _ __   ___| |_ _  ___  _ __  ___
%  |  __| | | | '_ \ / __| __| |/ _ \| '_ \/ __|
%  | |  | |_| | | | | (__| |_| | (_) | | | \__ \
%  |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%§
