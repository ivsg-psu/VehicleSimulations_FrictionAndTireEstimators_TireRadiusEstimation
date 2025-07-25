% script_test_fcn_TireRadiusEstimation_fillSampleData
% tests fcn_TireRadiusEstimation_fillSampleData.m

% Revision history
% 2025_07_09 by Sean Brennan
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

%% DEMO case: Basic demo
fig_num = 10001;
titleString = sprintf('DEMO case: Basic demo, loading all');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set inputs
specific_test_cases = [];
test_options = [];

% Call the function
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(fig_num));

sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOf_sampleData));

% Check variable sizes
Ntests = 2;
assert(isequal(Ntests,length(cellArrayOf_sampleData))); 

% Check variable values
% (too many to test)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% FIELD DATA case: HSOV Reber 3 laps position
fig_num = 7001;
titleString = sprintf('FIELD DATA case: HSOV Reber 3 laps position');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set inputs
specific_test_cases = 7001;
test_options = [];

% Call the function
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(fig_num));

legend('Interpreter','none');
sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOf_sampleData));
assert(isnumeric(cellArrayOf_sampleData{1}));

% Check variable sizes
Ntests = 1;
Ndata = 8933;
assert(isequal(Ntests,length(cellArrayOf_sampleData))); 
assert(isequal(size(cellArrayOf_sampleData{1}),[Ndata 5]));

% Check variable values
% (too much data to test)

% Make sure plot opened up
assert(isequal(get(gcf,'Number'),fig_num));

%% FIELD DATA case: HSOV Reber 3 laps velocity
fig_num = 7002;
titleString = sprintf('FIELD DATA case: HSOV Reber 3 laps velocity');
fprintf(1,'Figure %.0f: %s\n',fig_num, titleString);
figure(fig_num); clf;

% Set inputs
specific_test_cases = 7002;
test_options = [];

% Call the function
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(fig_num));

legend('Interpreter','none');
sgtitle(titleString, 'Interpreter','none');

% Check variable types
assert(iscell(cellArrayOf_sampleData));
assert(isnumeric(cellArrayOf_sampleData{1}));

% Check variable sizes
Ntests = 1;
Ndata = 8935;
assert(isequal(Ntests,length(cellArrayOf_sampleData))); 
assert(isequal(size(cellArrayOf_sampleData{1}),[Ndata 5]));

% Check variable values
% (too much data to test)

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

% Set inputs
specific_test_cases = [];
test_options = [];

% Call the function
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),([]));

% Check variable types
assert(iscell(cellArrayOf_sampleData));

% Check variable sizes
Ntests = 2;
assert(isequal(Ntests,length(cellArrayOf_sampleData))); 

% Check variable values
% (too many to test)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Basic fast mode - NO FIGURE, FAST MODE
fig_num = 80002;
fprintf(1,'Figure: %.0f: FAST mode, fig_num=-1\n',fig_num);
figure(fig_num); close(fig_num);

% Set inputs
specific_test_cases = [];
test_options = [];

% Call the function
cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));

% Check variable types
assert(iscell(cellArrayOf_sampleData));

% Check variable sizes
Ntests = 2;
assert(isequal(Ntests,length(cellArrayOf_sampleData))); 

% Check variable values
% (too many to test)

% Make sure plot did NOT open up
figHandles = get(groot, 'Children');
assert(~any(figHandles==fig_num));


%% Compare speeds of pre-calculation versus post-calculation versus a fast variant
fig_num = 80003;
fprintf(1,'Figure: %.0f: FAST mode comparisons\n',fig_num);
figure(fig_num);
close(fig_num);

% Set inputs
specific_test_cases = [];
test_options = [];

Niterations = 20;

% Do calculation without pre-calculation
tic;
for ith_test = 1:Niterations
    % Call the function
    cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),([]));
end
slow_method = toc;

% Do calculation with pre-calculation, FAST_MODE on
tic;
for ith_test = 1:Niterations
    % Call the function
    cellArrayOf_sampleData = fcn_TireRadiusEstimation_fillSampleData((specific_test_cases),(test_options),(-1));
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




