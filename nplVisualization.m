%% Load data and initialize NPL algorithm parameters
load("exp_values.mat");

% NPL algorithm configuration
M = 40;                        % Average power scaling factor
S = IM2d/max(IM2d);            % Normalized distortion sensitivity
gamma = max(IM2d);             % Maximum distortion coefficient
Delta = 0.5;                   % Power increment step
Plim = N * M;                  % Total power limit
N1 = 20;                       % Step at which to create detailed visualization

% Initialize algorithm variables
PCurr = zeros(1, N);           % Current power distribution
Pnpl = zeros(1, N);            % Final NPL power distribution
C = zeros(1, N);               % Capacity increments for each subcarrier
G2 = H.^1;                     % Channel gains (for calculation convenience)

% Data storage for visualization
PCtime = [];                   % Power distribution over time
Ctime = [];                    % Capacity over time
ND2time = [];                  % Distortion over time
SNRprojtime = [];              % SINR over time
nntime = [];                   % Selected subcarrier over time
it = 1;                        % Iteration counter

% Store detailed information for step N1
N1_detailed = struct();
N1_detailed.pre_PCurr = [];
N1_detailed.pre_SNR = [];
N1_detailed.subcarrier_tests = struct();
N1_detailed.selected_subcarrier = 0;
N1_detailed.post_PCurr = [];

%% Run NPL algorithm
while sum(PCurr) < Plim
    % Store current power distribution
    Pnpl = PCurr;

    % Store more detailed info if this is step N1
    if it == N1
        N1_detailed.pre_PCurr = PCurr;
        N1_detailed.pre_SNR = PCurr .* G2 ./ (N0 + gamma * S * (sum(PCurr .* G2)).^2);
    end
    
    % Evaluate capacity increment for each subcarrier
    for n = 1:N
        % Project power distribution with increment on nth subcarrier
        Pproj = PCurr;
        Pproj(n) = Pproj(n) + Delta;
        
        % Compute projected SINR and capacity
        SNRproj = Pproj .* G2 ./ (N0 + gamma * S * (sum(Pproj .* G2)).^2);
        C(n) = sum(log2(1 + SNRproj));
        
        % Store each subcarrier's projection for N1
        if it == N1
            N1_detailed.subcarrier_tests(n).Pproj = Pproj;
            N1_detailed.subcarrier_tests(n).SNRproj = SNRproj;
            N1_detailed.subcarrier_tests(n).capacity = C(n);
        end
    end

    % Find subcarrier with maximum capacity increment
    [~, nn] = max(C);
    nntime(it) = nn;
    
    % Store selected subcarrier for N1
    if it == N1
        N1_detailed.selected_subcarrier = nn;
    end

    % Allocate additional power to selected subcarrier
    PCurr(nn) = PCurr(nn) + Delta;
    
    % Store post-step state for N1
    if it == N1
        N1_detailed.post_PCurr = PCurr;
    end
    
    % Store data for visualization
    PCtime(it,:) = PCurr;
    Ctime(it,:) = C;
    ND2time(it,:) = gamma * S * (sum(PCurr .* G2)).^2;
    SNRprojtime(it,:) = SNRproj;
    it = it + 1;
end

%% Create regular animation (Part 1: steps 1 to N1)
fprintf('Generating regular animation (Part 1)...\n');

% Setup figure and animation parameters
f1 = figure('Color', 'w', 'Position', [100, 100, 1280, 720]);

% Visual styling
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], ...
          [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
markerSize = 6;
lineWidth = 1.5;
maxMarkerSize = 10;
fontSize = 12;
titleFontSize = 14;
axisFontSize = 12;

% Initialize video file for part 1
video_filename_part1 = 'npl_animation_part1.mp4';
if exist('v1', 'var'); close(v1); end
v1 = VideoWriter(video_filename_part1, 'MPEG-4');
v1.FrameRate = 24;
v1.Quality = 100;
open(v1);

% Create subplot layout and initialize plot objects
ax1 = subplot(3, 2, [1, 3, 5]); % Capacity plot
capPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
maxCapPlot = plot(1, 0, 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('C(n) [line 9 in Alg. 1]', 'FontSize', axisFontSize);
title({'Total capacity after incrementing power at subcarrier by \Delta'}, 'FontSize', titleFontSize);
grid on;
text(0.05, 0.9, "'X' = Currently Selected Subcarrier", 'Units', 'normalized', ...
     'FontSize', 16, 'Color', colors{2});
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

ax2 = subplot(3, 2, 2); % Power Distribution
powerPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
maxPowerPlot = plot(1, 0, 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
title('Power distribution', 'FontSize', titleFontSize, 'FontWeight', 'bold');
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
text(0.4, 0.2, {"'X' = Subcarrier Getting", "Power Increment"}, 'Units', 'normalized', ...
     'FontSize', 16, 'Color', colors{2});
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

ax3 = subplot(3, 2, 4); % SINR
sinrPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
xlim([1, N]);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('SINR (dB)', 'FontSize', axisFontSize);
title({'Current SINR', 'based on assigned powers'}, 'FontSize', titleFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

ax4 = subplot(3, 2, 6); % Distortion
distPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
noisePlot = plot(1:N, ones(1,N)*mean(N0(10:end-10)), '-', 'Color', colors{3}, 'LineWidth', 2.5);
hold off;
xlim([1, N]);
ylim([0 max([ND2time(end,:) N0(10:end-10)])]);
legend({'ND2(n)','N_0'}, 'FontSize', fontSize-1, 'Location', 'best', 'Box', 'off');
title('Distortion', 'FontSize', titleFontSize, 'FontWeight', 'bold');
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Level', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Add title for the whole figure
sgtitle('NPL algorithm iterations (exp. data)', 'FontSize', 16, 'FontWeight', 'bold');

% Define slow-down range for end of Part 1
slowdown_start = max(1, N1-12);  % Start slowing down 12 steps before N1
steps_part1 = 1:4:N1;
steps_count = length(steps_part1);

% Animation loop - only update data instead of recreating plots
for idx = 1:steps_count
    i = steps_part1(idx);
    
    % Update C(n) subplot
    set(capPlot, 'YData', Ctime(i,:));
    [~, ix] = max(Ctime(i,:));
    set(maxCapPlot, 'XData', ix, 'YData', Ctime(i,ix));
    
    % Calculate appropriate y-limits
    y_mean = mean(Ctime(i,:));
    y_range = max(0.1, 2*max(abs(y_mean-minmax(Ctime(i,:)))));
    ylim(ax1, [y_mean-y_range/2, y_mean+y_range/2]);
    
    % Update Power Distribution subplot
    power_data = 10*log10(PCtime(i,:));
    set(powerPlot, 'YData', power_data);
    set(maxPowerPlot, 'XData', ix, 'YData', power_data(ix));
    
    % Calculate appropriate y-limits for power
    power_min = min(power_data);
    power_max = min(20, max(power_data));
    power_padding = 0.1 * max(0.1, power_max - power_min);
    ylim(ax2, [max(0, power_min-power_padding), min(20, power_max+power_padding)]);
    
    % Update SINR subplot
    snr_data = 10*log10(SNRprojtime(i,:));
    set(sinrPlot, 'YData', snr_data);
    
    % Calculate appropriate y-limits for SINR
    snr_min = min(snr_data);
    snr_max = max(snr_data);
    snr_padding = 0.1 * (snr_max - snr_min);
    ylim(ax3, [max(10, snr_min-snr_padding), snr_max+snr_padding]);
    
    % Update Distortion subplot
    set(distPlot, 'YData', ND2time(i,:));
    
    % Render changes and capture frame
    drawnow;
    
    % Calculate how many times to repeat the frame based on position
    % This creates a slowdown effect approaching the end of Part 1
    repeat_count = 1;  % Default is 1 (no repeat)
    
    % If we're in the slow-down range, gradually increase repeat_count
    if i >= slowdown_start
        % Calculate how far we are into the slowdown range (0 to 1)
        progress = (i - slowdown_start) / (N1 - slowdown_start);
        % Exponential slowdown: start with 1x speed, end with 8x slower
        repeat_count = ceil(1 + 7 * (progress^2));
    end
    
    % Repeat the frame to create slowdown effect
    frame = getframe(gcf);
    for r = 1:repeat_count
        writeVideo(v1, frame);
    end
    
    % For the very last frame, hold it longer
    if i == steps_part1(end)
        for r = 1:10
            writeVideo(v1, frame);
        end
    end
end

% Close video file
close(v1);
fprintf('Regular animation part 1 saved to: %s\n', video_filename_part1);

%% Create detailed animation for step N1
fprintf('Generating detailed animation for step N1...\n');

% Close previous figure and create a new one for detailed animation
close(f1);
f2 = figure('Color', 'w', 'Position', [100, 100, 1280, 720]);

% Animation parameters
detailed_video_filename = 'npl_detailed_animation.mp4';
if exist('v3', 'var'); close(v3); end
v3 = VideoWriter(detailed_video_filename, 'MPEG-4');
v3.FrameRate = 5;  % Slower framerate for detailed explanation
v3.Quality = 100;
open(v3);

% Create explanatory titles for each step
step_titles = {
    'Step 1: Start with current power distribution',
    'Step 2: For each subcarrier, calculate capacity increase after adding power',
    'Step 3: Find subcarrier with maximum capacity increase',
    'Step 4: Add power to the selected subcarrier',
    'Step 5: Update power distribution and proceed to next iteration'
};

% Get data for step N1
step_data = N1_detailed;

% --- STEP 1: Show initial state ---
% Create 2x2 subplot layout for detailed visualization
ax1 = subplot(2, 2, 1);  % Power Distribution
pDist = plot(1:N, 10*log10(step_data.pre_PCurr), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
xlim([1, N]);
title('Initial power distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

ax2 = subplot(2, 2, 2);  % SINR
sinrDist = plot(1:N, 10*log10(step_data.pre_SNR), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
xlim([1, N]);
title('Initial SINR distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('SINR (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Algorithm display - STEP 1
ax3 = subplot(2, 2, 3);
cla;
axis off;

% Draw algorithm text in a contained textbox
text_string = sprintf(['NPL algorithm:\n' ...
                       '1. P_{NPL}(n) = 0, n ∈ 1:N\n' ...
                       '2. While ∑P_{NPL} < P_{lim}\n' ...
                       'a. For each n ∈ 1:N:\n' ...
                       '   - P_{proj} = P_{NPL}\n' ...
                       '   - P_{proj}(n) += Δ\n' ...
                       '   - Calculate SNR and C(n)\n' ...
                       'b. Select m = argmax C(n)\n' ...
                       'c. P_{NPL}(m) += Δ\n' ...
                       '3. Return P_{NPL}']);

% First draw the main text
algText = text(0.05, 0.95, text_string, 'VerticalAlignment', 'top', 'FontSize', 12);

ax4 = subplot(2, 2, 4);  % Blank placeholder for capacity
axis off;
capacityInfo = text(0.5, 0.5, {'Capacity calculation will be shown','in the next steps'}, ...
    'HorizontalAlignment', 'center', 'FontSize', 14);

% Add explanation text
expText = annotation('textbox', [0.1, 0.01, 0.8, 0.05], 'String', ...
    'Step 1: Algorithm starts with the current power distribution and calculates initial SINR', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

sgtitle(sprintf('Step %d: Detailed NPL iterations\n%s', N1, step_titles{1}), 'FontSize', 16, 'FontWeight', 'bold');

% Capture frame
drawnow;
for i = 1:15  % Hold this frame longer
    frame = getframe(gcf);
    writeVideo(v3, frame);
end

% --- STEP 2: Calculate capacity for each subcarrier ---
temp_C = zeros(1, N);

% Create plots once
subplot(2, 2, 1);
origPower = plot(1:N, 10*log10(step_data.pre_PCurr), 'o', 'MarkerSize', markerSize, 'Color', [0.7, 0.7, 0.7], ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'none');
hold on;
projPower = plot(1:N, 10*log10(step_data.pre_PCurr), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
currSubc = plot(1, 10*log10(step_data.pre_PCurr(1)), 'o', 'MarkerSize', 15, ...
    'Color', 'r', 'LineWidth', 2);
hold off;
xlim([1, N]);
title('Projected power (Testing subcarrier 1)', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

subplot(2, 2, 2);
yyaxis left;
origSNR = plot(1:N, 10*log10(step_data.pre_SNR), 'o', 'MarkerSize', markerSize, 'Color', [0.7, 0.7, 0.7], ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'none');
hold on;
projSNR = plot(1:N, 10*log10(step_data.pre_SNR), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold off;
ylabel('SINR (dB)', 'FontSize', axisFontSize);
set(gca, 'YColor', 'k');

yyaxis right;
capacityPerCarrier = plot(1:N, zeros(1,N), '.', 'MarkerSize', 10, 'Color', colors{3});
ylabel('Capacity per subcarrier', 'FontSize', axisFontSize, 'Color', colors{3});
set(gca, 'YColor', colors{3});

xlim([1, N]);
ylim([floor(min(log2(1+[step_data.subcarrier_tests(:).SNRproj]))) ceil(max(log2(1+[step_data.subcarrier_tests(:).SNRproj])))]);
title('Projected SINR & capacity (Testing subcarrier 1)', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

subplot(2, 2, 3);
cla;
axis off;
algText2 = text(0.05, 0.95, '', 'VerticalAlignment', 'top', 'FontSize', 12);

subplot(2, 2, 4);
capProgress = plot(1:N, ones(1,N)*NaN, 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
currCap = plot(1, 0, 'o', 'MarkerSize', 15, 'Color', 'r', 'LineWidth', 2);
hold off;
xlim([1, N]);
title('Capacity after incrementing power', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('\Sigma(Capacity per subc.)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update plots for each subcarrier
for n = 1:N
    % Update temporary capacity array
    temp_C(n) = step_data.subcarrier_tests(n).capacity;
    
    % Update titles
    sgtitle(sprintf('Step %d: Detailed NPL iterations\n%s', N1, step_titles{2}), 'FontSize', 16, 'FontWeight', 'bold');
    subplot(2, 2, 1); title(sprintf('Projected Power (Testing Subcarrier %d)', n), 'FontSize', titleFontSize);
    subplot(2, 2, 2); title(sprintf('Projected SINR & Capacity (Testing Subcarrier %d)', n), 'FontSize', titleFontSize);
    
    % Update power plot
    set(projPower, 'YData', 10*log10(step_data.subcarrier_tests(n).Pproj));
    set(currSubc, 'XData', n, 'YData', 10*log10(step_data.subcarrier_tests(n).Pproj(n)));
    
    % Update SINR and capacity plots
    yyaxis left;
    set(projSNR, 'YData', 10*log10(step_data.subcarrier_tests(n).SNRproj));
    
    yyaxis right;
    capacity_per_carrier = log2(1 + step_data.subcarrier_tests(n).SNRproj);
    set(capacityPerCarrier, 'YData', capacity_per_carrier);
    
    % Update algorithm text
    text_string = sprintf(['NPL algorithm:\n' ...
                       '1. P_{NPL}(n) = 0, n ∈ 1:N\n' ...
                       '2. While ∑P_{NPL} < P_{lim}\n' ...
                       '\\bf a. For each n ∈ 1:N:\\rm\n' ...
                       '\\bf   - P_{proj} = P_{NPL}\\rm\n' ...
                       '\\bf   - P_{proj}(n) += Δ\\rm\n' ...
                       '\\bf   - Calculate SNR and C(n)\\rm\n' ...
                       'b. Select m = argmax C(n)\n' ...
                       'c. P_{NPL}(m) += Δ\n' ...
                       '3. Return P_{NPL}']);
    set(algText2, 'String', text_string);
    
    % Update capacity progress
    y_data = ones(1,N)*NaN;
    y_data(1:n) = temp_C(1:n);
    set(capProgress, 'YData', y_data);
    set(currCap, 'XData', n, 'YData', temp_C(n));
    
    % Update y-limits for capacity
    if n > 1
        c_min = min(temp_C(1:n));
        c_max = max(temp_C(1:n));
        c_padding = 0.1 * (c_max - c_min);
        if c_max - c_min < 0.1
            c_padding = 0.05;
        end
        ylim(subplot(2, 2, 4), [c_min-c_padding, c_max+c_padding]);
    end
    
    % Update explanation text
    set(expText, 'String', sprintf('Step 2: Testing subcarrier %d by adding power Δ.', n));
    
    % Render and capture frame
    drawnow;
    frame = getframe(gcf);
    writeVideo(v3, frame);
end

% --- STEP 3: Find maximum capacity ---
[max_cap, max_idx] = max(temp_C);

% Update title
sgtitle(sprintf('Step %d: Detailed NPL iterations\n%s', N1, step_titles{3}), 'FontSize', 16, 'FontWeight', 'bold');

% Update power plot
subplot(2, 2, 1);
cla;
plot(1:N, 10*log10(step_data.subcarrier_tests(max_idx).Pproj), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, 10*log10(step_data.subcarrier_tests(max_idx).Pproj(max_idx)), 'x', 'MarkerSize', 15, ...
    'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
title('Chosen power distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update SINR plot
subplot(2, 2, 2);
cla;
yyaxis left;
plot(1:N, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj(max_idx)), 'x', 'MarkerSize', 15, ...
    'Color', colors{2}, 'LineWidth', 3);
hold off;
ylabel('SINR (dB)', 'FontSize', axisFontSize);
set(gca, 'YColor', 'k');
ylim([floor(min(10*log10([step_data.subcarrier_tests(:).SNRproj]))) ceil(max(10*log10([step_data.subcarrier_tests(:).SNRproj])))]);

yyaxis right;
capacity_per_carrier = log2(1 + step_data.subcarrier_tests(max_idx).SNRproj);
plot(1:N, capacity_per_carrier, '.', 'MarkerSize', 10, 'Color', colors{3});
ylabel('Capacity per subcarrier', 'FontSize', axisFontSize, 'Color', colors{3});
set(gca, 'YColor', colors{3});
ylim([floor(min(log2(1+[step_data.subcarrier_tests(:).SNRproj]))) ceil(max(log2(1+[step_data.subcarrier_tests(:).SNRproj])))]);

xlim([1, N]);
title('Chosen SINR & capacity distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update algorithm text
subplot(2, 2, 3);
cla;
axis off;
text_string = sprintf(['NPL algorithm:\n' ...
                   '1. P_{NPL}(n) = 0, n ∈ 1:N\n' ...
                   '2. While ∑P_{NPL} < P_{lim}\n' ...
                   'a. For each n ∈ 1:N:\n' ...
                   '   - P_{proj} = P_{NPL}\n' ...
                   '   - P_{proj}(n) += Δ\n' ...
                   '   - Calculate SNR and C(n)\n' ...
                   '\\bfb. Select m = argmax C(n)\\rm\n' ...
                   'c. P_{NPL}(m) += Δ\n' ...
                   '3. Return P_{NPL}']);
text(0.05, 0.95, text_string, 'VerticalAlignment', 'top', 'FontSize', 12);

% Update capacity plot
subplot(2, 2, 4);
cla;
plot(1:N, temp_C, 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, max_cap, 'x', 'MarkerSize', 15, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
c_min = min(temp_C);
c_max = max(temp_C);
c_padding = 0.1 * (c_max - c_min);
ylim([c_min-c_padding, c_max+c_padding]);
title('Capacity after incrementing power', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('\Sigma(Capacity per subc.)', 'FontSize', axisFontSize);
grid on;
text(0.05, 0.9, sprintf("Max capacity = %.4f at subcarrier %d", max_cap, max_idx), ...
     'Units', 'normalized', 'FontSize', 11, 'Color', 'k');
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update explanation text
set(expText, 'String', sprintf('Step 3: Found maximum capacity %.4f at subcarrier %d', max_cap, max_idx));

% Render and capture frame
drawnow;
for i = 1:15  % Hold this frame longer
    frame = getframe(gcf);
    writeVideo(v3, frame);
end

% --- STEP 4: Add power to selected subcarrier ---
% Update title
sgtitle(sprintf('Step %d: Detailed NPL iterations\n%s', N1, step_titles{4}), 'FontSize', 16, 'FontWeight', 'bold');

% Update power plot with arrow
subplot(2, 2, 1);
cla;
plot(1:N, 10*log10(step_data.pre_PCurr), 'o', 'MarkerSize', markerSize, 'Color', [0.7, 0.7, 0.7], ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'none');
hold on;
plot(1:N, 10*log10(step_data.post_PCurr), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
y_data_old = 10*log10(step_data.pre_PCurr(max_idx));
y_data_new = 10*log10(step_data.post_PCurr(max_idx));
quiver(max_idx, y_data_old, 0, y_data_new-y_data_old, 0, ...
    'LineWidth', 2, 'Color', 'r', 'MaxHeadSize', 1);
text(max_idx+2, (y_data_old+y_data_new)/2, ['+', num2str(Delta), 'dB'], ...
    'Color', 'r', 'FontWeight', 'bold');
hold off;
xlim([1, N]);
title('Power increment applied', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update SINR plot
subplot(2, 2, 2);
cla;
yyaxis left;
plot(1:N, 10*log10(step_data.pre_SNR), 'o', 'MarkerSize', markerSize, 'Color', [0.7, 0.7, 0.7], ...
    'MarkerFaceColor', [0.7, 0.7, 0.7], 'MarkerEdgeColor', 'none');
hold on;
plot(1:N, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
plot(max_idx, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj(max_idx)), 'x', 'MarkerSize', 15, ...
    'Color', colors{2}, 'LineWidth', 3);
hold off;
ylabel('SINR (dB)', 'FontSize', axisFontSize);
set(gca, 'YColor', 'k');
ylim([floor(min(10*log10([step_data.subcarrier_tests(:).SNRproj]))) ceil(max(10*log10([step_data.subcarrier_tests(:).SNRproj])))]);

yyaxis right;
capacity_per_carrier = log2(1 + step_data.subcarrier_tests(max_idx).SNRproj);
plot(1:N, capacity_per_carrier, '.', 'MarkerSize', 10, 'Color', colors{3});
ylabel('Capacity per subcarrier', 'FontSize', axisFontSize, 'Color', colors{3});
set(gca, 'YColor', colors{3});
ylim([floor(min(log2(1+[step_data.subcarrier_tests(:).SNRproj]))) ceil(max(log2(1+[step_data.subcarrier_tests(:).SNRproj])))]);

xlim([1, N]);
title('Updated SINR & capacity', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update algorithm text
subplot(2, 2, 3);
cla;
axis off;
text_string = sprintf(['NPL algorithm:\n' ...
                   '1. P_{NPL}(n) = 0, n ∈ 1:N\n' ...
                   '2. While ∑P_{NPL} < P_{lim}\n' ...
                   'a. For each n ∈ 1:N:\n' ...
                   '   - P_{proj} = P_{NPL}\n' ...
                   '   - P_{proj}(n) += Δ\n' ...
                   '   - Calculate SNR and C(n)\n' ...
                   'b. Select m = argmax C(n)\n' ...
                   '\\bfc. P_{NPL}(m) += Δ\\rm\n' ...
                   '3. Return P_{NPL}']);
text(0.05, 0.95, text_string, 'VerticalAlignment', 'top', 'FontSize', 12);

% Update capacity plot
subplot(2, 2, 4);
cla;
plot(1:N, temp_C, 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, max_cap, 'x', 'MarkerSize', 15, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
title('Final capacity distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('\Sigma(Capacity per subc.)', 'FontSize', axisFontSize);
grid on;
text(0.05, 0.9, sprintf("Selected capacity = %.4f at subcarrier %d", max_cap, max_idx), ...
     'Units', 'normalized', 'FontSize', 11, 'Color', 'k');
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update explanation text
set(expText, 'String', sprintf('Step 4: Increasing power on subcarrier %d by Δ = %.2f', max_idx, Delta));

% Render and capture frame
drawnow;
for i = 1:15  % Hold this frame longer
    frame = getframe(gcf);
    writeVideo(v3, frame);
end

% --- STEP 5: Final State and Next Iteration ---
% Update title
sgtitle(sprintf('Step %d: Detailed NPL iterations\n%s', N1, step_titles{5}), 'FontSize', 16, 'FontWeight', 'bold');

% Update power plot
subplot(2, 2, 1);
cla;
plot(1:N, 10*log10(step_data.post_PCurr), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, 10*log10(step_data.post_PCurr(max_idx)), 'x', 'MarkerSize', 15, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
title('Updated power distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update SINR plot
subplot(2, 2, 2);
cla;
yyaxis left;
plot(1:N, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, 10*log10(step_data.subcarrier_tests(max_idx).SNRproj(max_idx)), 'x', 'MarkerSize', 15, 'Color', colors{2}, 'LineWidth', 3);
hold off;
ylabel('SINR (dB)', 'FontSize', axisFontSize);
set(gca, 'YColor', 'k');
ylim([floor(min(10*log10([step_data.subcarrier_tests(:).SNRproj]))) ceil(max(10*log10([step_data.subcarrier_tests(:).SNRproj])))]);

yyaxis right;
capacity_per_carrier = log2(1 + step_data.subcarrier_tests(max_idx).SNRproj);
plot(1:N, capacity_per_carrier, '.', 'MarkerSize', 10, 'Color', colors{3});
ylabel('Capacity per subcarrier', 'FontSize', axisFontSize, 'Color', colors{3});
set(gca, 'YColor', colors{3});
ylim([floor(min(log2(1+[step_data.subcarrier_tests(:).SNRproj]))) ceil(max(log2(1+[step_data.subcarrier_tests(:).SNRproj])))]);

xlim([1, N]);
title('Final SINR & capacity distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update algorithm text
subplot(2, 2, 3);
cla;
axis off;
text_string = sprintf(['NPL algorithm:\n' ...
                   '1. P_{NPL}(n) = 0, n ∈ 1:N\n' ...
                   '\\bf2. While ∑P_{NPL} < P_{lim}\\rm\n' ...
                   'a. For each n ∈ 1:N:\n' ...
                   '   - P_{proj} = P_{NPL}\n' ...
                   '   - P_{proj}(n) += Δ\n' ...
                   '   - Calculate SNR and C(n)\n' ...
                   'b. Select m = argmax C(n)\n' ...
                   'c. P_{NPL}(m) += Δ\n' ...
                   '3. Return P_{NPL}']);
text(0.05, 0.95, text_string, 'VerticalAlignment', 'top', 'FontSize', 12);

% Update capacity plot
subplot(2, 2, 4);
cla;
plot(1:N, temp_C, 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none');
hold on;
plot(max_idx, max_cap, 'x', 'MarkerSize', 15, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
title('Final capacity distribution', 'FontSize', titleFontSize);
xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
ylabel('\Sigma(Capacity per subc.)', 'FontSize', axisFontSize);
grid on;
set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Update explanation text
set(expText, 'String', sprintf('Step 5: Iteration %d completed. NPL algorithm continues to next iteration...', N1));

% Render and capture frame
drawnow;
for i = 1:20  % Hold this frame longer
    frame = getframe(gcf);
    writeVideo(v3, frame);
end

% Close detailed video
close(v3);
fprintf('Detailed animation saved to: %s\n', detailed_video_filename);

% Clean up figures
close(f2);

%% Create regular animation (Part 2: steps N1+1 to end)
if N1 < size(PCtime,1)
    fprintf('Generating regular animation (Part 2)...\n');
    
    % Re-create the figure for Part 2 (same layout as Part 1)
    f1 = figure('Color', 'w', 'Position', [100, 100, 1280, 720]);
    
    % Initialize video file for part 2
    video_filename_part2 = 'npl_animation_part2.mp4';
    if exist('v2', 'var'); close(v2); end
    v2 = VideoWriter(video_filename_part2, 'MPEG-4');
    v2.FrameRate = 24;
    v2.Quality = 100;
    open(v2);
    
    % Set up the same plots as in Part 1
    ax1 = subplot(3, 2, [1, 3, 5]); % Capacity plot
    capPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
        'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
    hold on;
    maxCapPlot = plot(1, 0, 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
    hold off;
    xlim([1, N]);
    xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
    ylabel('C(n) [line 9 in Alg. 1]', 'FontSize', axisFontSize);
    title({'Total capacity after incrementing power at subcarrier by \Delta'}, 'FontSize', titleFontSize);
    grid on;
    text(0.05, 0.9, "'X' = Currently Selected Subcarrier", 'Units', 'normalized', ...
         'FontSize', 16, 'Color', colors{2});
    set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

    ax2 = subplot(3, 2, 2); % Power Distribution
    powerPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
        'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
    hold on;
    maxPowerPlot = plot(1, 0, 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
    hold off;
    xlim([1, N]);
    title('Power distribution', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
    ylabel('Power (dB)', 'FontSize', axisFontSize);
    grid on;
    text(0.4, 0.2, {"'X' = Subcarrier Getting", "Power Increment"}, 'Units', 'normalized', ...
         'FontSize', 16, 'Color', colors{2});
    set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

    ax3 = subplot(3, 2, 4); % SINR
    sinrPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
        'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
    xlim([1, N]);
    xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
    ylabel('SINR (dB)', 'FontSize', axisFontSize);
    title({'Current SINR', 'based on assigned powers'}, 'FontSize', titleFontSize);
    grid on;
    set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

    ax4 = subplot(3, 2, 6); % Distortion
    distPlot = plot(1:N, zeros(1,N), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
        'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
    hold on;
    noisePlot = plot(1:N, ones(1,N)*mean(N0(10:end-10)), '-', 'Color', colors{3}, 'LineWidth', 2.5);
    hold off;
    xlim([1, N]);
    ylim([0 max([ND2time(end,:) N0(10:end-10)])]);
    legend({'ND2(n)','N_0'}, 'FontSize', fontSize-1, 'Location', 'best', 'Box', 'off');
    title('Distortion', 'FontSize', titleFontSize, 'FontWeight', 'bold');
    xlabel('Subcarrier number [#]', 'FontSize', axisFontSize);
    ylabel('Level', 'FontSize', axisFontSize);
    grid on;
    set(gca, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

    % Add title for the whole figure
    sgtitle('NPL algorithm iterations (exp. data)', 'FontSize', 16, 'FontWeight', 'bold');
    
    % Define steps for part 2 animation
    steps_part2 = N1+1:4:size(PCtime,1);
    steps_count = length(steps_part2);
    
    % Define acceleration range for beginning of Part 2
    speedup_end = min(steps_count, 4);  % How many iterations to gradually speed up
    
    % Animation loop with slowdown at start
    for idx = 1:steps_count
        i = steps_part2(idx);
        
        % Update C(n) subplot
        set(capPlot, 'YData', Ctime(i,:));
        [~, ix] = max(Ctime(i,:));
        set(maxCapPlot, 'XData', ix, 'YData', Ctime(i,ix));
        
        % Calculate appropriate y-limits
        y_mean = mean(Ctime(i,:));
        y_range = max(0.1, 2*max(abs(y_mean-minmax(Ctime(i,:)))));
        ylim(ax1, [y_mean-y_range/2, y_mean+y_range/2]);
        
        % Update Power Distribution subplot
        power_data = 10*log10(PCtime(i,:));
        set(powerPlot, 'YData', power_data);
        set(maxPowerPlot, 'XData', ix, 'YData', power_data(ix));
        
        % Calculate appropriate y-limits for power
        power_min = min(power_data);
        power_max = min(20, max(power_data));
        power_padding = 0.1 * max(0.1, power_max - power_min);
        ylim(ax2, [max(0, power_min-power_padding), min(20, power_max+power_padding)]);
        
        % Update SINR subplot
        snr_data = 10*log10(SNRprojtime(i,:));
        set(sinrPlot, 'YData', snr_data);
        
        % Calculate appropriate y-limits for SINR
        snr_min = min(snr_data);
        snr_max = max(snr_data);
        snr_padding = 0.1 * (snr_max - snr_min);
        ylim(ax3, [max(10, snr_min-snr_padding), snr_max+snr_padding]);
        
        % Update Distortion subplot
        set(distPlot, 'YData', ND2time(i,:));
        
        % Render changes and capture frame
        drawnow;
        
        % Calculate slowdown factor for the beginning
        repeat_count = 1; % Default is 1 (normal speed)
        
        if idx <= speedup_end
            % Start slow and gradually speed up
            progress = idx / speedup_end;
            repeat_count = ceil(8 * (1 - progress^2));  % Start at 8x slower, exponentially approach normal speed
        end
        
        % For the final frames, hold longer
        if idx >= steps_count - 2
            repeat_count = max(repeat_count, 5);  % At least 5x longer for the last frames
        end
        
        % Repeat the frame to create slowdown/lingering effect
        frame = getframe(gcf);
        for r = 1:repeat_count
            writeVideo(v2, frame);
        end
        
        % For the very last frame, hold it even longer
        if idx == steps_count
            for r = 1:20
                writeVideo(v2, frame);
            end
        end
    end
    
    % Close video file
    close(v2);
    fprintf('Regular animation part 2 saved to: %s\n', video_filename_part2);
    close(f1);
end

fprintf('All animations have been generated successfully!\n');