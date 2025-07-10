% script_test_fcn_TireRadiusEstimation_rEstVelFilteredTime
% tests fcn_TireRadiusEstimation_rEstVelFilteredTime

% Revision history
% 2025_07_10 by Sean Brennan
% -- first write of the code

%% Set up the workspace
close all

%% Code demos start here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   _____                              ____   __    _____          _
%  |  __ \                            / __ \ / _|  / ____|        | |
%  | |  | | ___ _ __ ___   ___  ___  | |  | | |_  | |     ___   __| | ___
%  | |  | |/ _ \ '_ ` _ \ / _ \/ __| | |  | |  _| | |    / _ \ / _` |/ _ \
%  | |__| |  __/ | | | | | (_) \__ \ | |__| | |   | |___| (_) | (_| |  __/
%  |_____/ \___|_| |_| |_|\___/|___/  \____/|_|    \_____\___/ \__,_|\___|
%
%
% See: https://patorjk.com/software/taag/#p=display&f=Big&t=Demos%20Of%20Code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 1

close all;
fprintf(1,'Figure: 1XXXXXX: DEMO cases\n');

%% DEMO case: Basic demo with perfect data
fig_num = 10001;
titleString = sprintf('DEMO case: Basic demo with perfect data');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load test data
specific_test_cases = 2;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = [];
plotXlabelString = [];

% Call the function
rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rEff));

% Check variable sizes
Npoints = length(velocities_omegas(:,1));
assert(isequal(size(rEff),[Npoints 1])); 

% Check variable values
assert(max(abs(rEff-1))<0.00001);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));


%% DEMO case: Basic demo using optional inputs, perfect data
fig_num = 10002;
titleString = sprintf('DEMO case: Basic demo using optional inputs, perfect data');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load test data
specific_test_cases = 2;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = (1:length(velocities_omegas(:,1)))'*0.1;
plotXlabelString = 'Time [sec]';

% Call the function
rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rEff));

% Check variable sizes
Npoints = length(velocities_omegas(:,1));
assert(isequal(size(rEff),[Npoints 1])); 

% Check variable values
assert(max(abs(rEff-1))<0.00001);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% FIELD DATA case: HSOV Reber 3 laps
fig_num = 7002;
titleString = sprintf('FIELD DATA case: HSOV Reber 3 laps');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Load test data
specific_test_cases = 7002;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = (1:length(velocities_omegas(:,1)))'*0.05;
plotXlabelString = 'Time [sec]';

% Call the function
rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), (fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(isnumeric(rEff));

% Check variable sizes
Npoints = length(velocities_omegas(:,1));
Nencoders = length(velocities_omegas(1,:))-1;
assert(isequal(size(rEff),[Npoints Nencoders])); 

% Check variable values
assert(max(abs(rEff-0.09),[],'all')<0.1);

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% Fast Mode Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ______        _     __  __           _        _______        _
% |  ____|      | |   |  \/  |         | |      |__   __|      | |
% | |__ __ _ ___| |_  | \  / | ___   __| | ___     | | ___  ___| |_ ___
% |  __/ _` / __| __| | |\/| |/ _ \ / _` |/ _ \    | |/ _ \/ __| __/ __|
% | | | (_| \__ \ |_  | |  | | (_) | (_| |  __/    | |  __/\__ \ |_\__ \
% |_|  \__,_|___/\__| |_|  |_|\___/ \__,_|\___|    |_|\___||___/\__|___/
%
%
% See: http://patorjk.com/software/taag/#p=display&f=Big&t=Fast%20Mode%20Tests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures start with 8

close all;
fprintf(1,'Figure: 8XXXXXX: FAST mode cases\n');

%% Basic example - NO FIGURE
fig_num = 80001;
fprintf(1,'Figure: %.0f: FAST mode, empty fig_num\n',fig_num);
figure(fig_num); close(fig_num);

% Load test data
specific_test_cases = 2;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = [];
plotXlabelString = [];

% Call the function
rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), ([]));

% Check variable types
assert(isnumeric(rEff));

% Check variable sizes
Npoints = length(velocities_omegas(:,1));
assert(isequal(size(rEff),[Npoints 1])); 

% Check variable values
assert(max(abs(rEff-1))<0.00001);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Load test data
specific_test_cases = 2;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = [];
plotXlabelString = [];

% Call the function
rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), (-1));

% Check variable types
assert(isnumeric(rEff));

% Check variable sizes
Npoints = length(velocities_omegas(:,1));
assert(isequal(size(rEff),[Npoints 1])); 

% Check variable values
assert(max(abs(rEff-1))<0.00001);

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Load test data
specific_test_cases = 2;
test_options = [];
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
velocities_omegas = cellArrayOf_sampleData{1};

% Set values
plotXvalues = [];
plotXlabelString = [];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), ([]));

end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    rEff = fcn_TireRadiusEstimation_rEstVelFilteredTime(velocities_omegas, (plotXvalues), (plotXlabelString), (-1));

end
fast_method = toc;

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));

% Plot results as bar chart
figure(373737);
clf;
hold on;

X = categorical({'Normal mode','Fast mode'});
X = reordercats(X,{'Normal mode','Fast mode'}); % Forces bars to appear in this exact order, not alphabetized
Y = [slow_method fast_method ]*1000/Niterations;
bar(X,Y)
ylabel('Execution time (Milliseconds)')


% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% BUG cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  ____  _    _  _____
% |  _ \| |  | |/ ____|
% | |_) | |  | | |  __    ___ __ _ ___  ___  ___
% |  _ <| |  | | | |_ |  / __/ _` / __|/ _ \/ __|
% | |_) | |__| | |__| | | (_| (_| \__ \  __/\__ \
% |____/ \____/ \_____|  \___\__,_|___/\___||___/
%
% See: http://patorjk.com/software/taag/#p=display&v=0&f=Big&t=BUG%20cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All bug case figures start with the number 9

% close all;

%% BUG 

%% Fail conditions
if 1==0
    %
    
end


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




