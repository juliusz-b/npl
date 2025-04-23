%% Load data and initialize NPL algorithm parameters
load("exp_values.mat");

% NPL algorithm configuration
M = 40;                        % Average power scaling factor
S = IM2d/max(IM2d);            % Normalized distortion sensitivity
gamma = max(IM2d);             % Maximum distortion coefficient
Delta = 0.5;                   % Power increment step
Plim = N * M;                  % Total power limit

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

%% Run NPL algorithm
while sum(PCurr) < Plim
    % Store current power distribution
    Pnpl = PCurr;

    % Evaluate capacity increment for each subcarrier
    for n = 1:N
        % Project power distribution with increment on nth subcarrier
        Pproj = PCurr;
        Pproj(n) = Pproj(n) + Delta;
        
        % Compute projected SINR and capacity
        SNRproj = Pproj .* G2 ./ (N0 + gamma * S * (sum(Pproj .* G2)).^2);
        C(n) = sum(log2(1 + SNRproj));
    end

    % Find subcarrier with maximum capacity increment
    [~, nn] = max(C);
    nntime(it) = nn;

    % Allocate additional power to selected subcarrier
    PCurr(nn) = PCurr(nn) + Delta;
    
    % Store data for visualization
    PCtime(it,:) = PCurr;
    Ctime(it,:) = C;
    ND2time(it,:) = gamma * S * (sum(PCurr .* G2)).^2;
    SNRprojtime(it,:) = SNRproj;
    it = it + 1;
end

%% Create visualization
% Setup figure and animation parameters
figure('Color', 'w');
animation_length_seconds = 25;
enable_sound = true;

% Visual styling
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], ...
          [0.9290, 0.6940, 0.1250], [0.4940, 0.1840, 0.5560]};
markerSize = 6;
lineWidth = 1.5;
maxMarkerSize = 10;
fontSize = 12;
titleFontSize = 14;
axisFontSize = 12;

% Adaptive y-axis configuration
smoothing_factor_snr = 0.85;
smoothing_factor_power = 0.75;
min_power_range = 3.0;

% Sound parameters
sample_rate = 22050;
base_freq = 220;
freq_range = 1320;
tone_duration = 0.10;
tone_time = 0:1/sample_rate:tone_duration;
volume_scale = 0.5;
prev_ix = 0;

% Create tiled layout
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(t, 'NPL Algorithm Iterations', 'FontSize', 16, 'FontWeight', 'bold');

% Tile 1: Capacity
ax1 = nexttile;
p1 = plot(1:N, Ctime(1,:), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
[~, ix] = max(Ctime(1,:));
p1_max = plot(ix, Ctime(1,ix), 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);
xlabel('Subcarrier Number [#]', 'FontSize', axisFontSize);
ylabel('C(n) [line 9 in Alg. 1]', 'FontSize', axisFontSize);
grid on;
text(0.05, 0.9, "'X' = Currently Selected Subcarrier", 'Units', 'normalized', ...
     'FontSize', 9, 'Color', colors{2});
set(ax1, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);
cap_ax = ax1;

% Tile 2: SNR
ax2 = nexttile;
p2 = plot(1:N, 10*log10(SNRprojtime(1,:)), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
xlim([1, N]);

% Initialize SNR y-limits
snr_data_1 = 10*log10(SNRprojtime(1,:));
snr_min = min(snr_data_1);
snr_max = max(snr_data_1);
snr_padding = 0.1 * (snr_max - snr_min);
current_snr_min = max(10, snr_min - snr_padding);
current_snr_max = snr_max + snr_padding;
ylim(ax2, [current_snr_min, current_snr_max]);

xlabel('Subcarrier Number [#]', 'FontSize', axisFontSize);
ylabel('SINR (dB)', 'FontSize', axisFontSize);
grid on;
set(ax2, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Tile 3: Power
ax3 = nexttile;
p3 = plot(1:N, 10*log10(PCtime(1,:)), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
p3_max = plot(ix, 10*log10(PCtime(1,ix)), 'x', 'MarkerSize', maxMarkerSize, 'Color', colors{2}, 'LineWidth', 3);
hold off;
xlim([1, N]);

% Initialize Power y-limits
power_data_1 = 10*log10(PCtime(1,:));
power_min = min(power_data_1);
power_max = min(20, max(power_data_1));

% Ensure minimum range
if (power_max - power_min) < min_power_range
    power_center = (power_max + power_min) / 2;
    power_min = power_center - min_power_range/2;
    power_max = power_center + min_power_range/2;
end

power_padding = 0.1 * (power_max - power_min);
current_power_min = max(0, power_min - power_padding);
current_power_max = min(20, power_max + power_padding);

if (current_power_max - current_power_min) < min_power_range
    power_center = (current_power_max + current_power_min) / 2;
    current_power_min = max(0, power_center - min_power_range/2);
    current_power_max = min(20, power_center + min_power_range/2);
end

ylim(ax3, [current_power_min, current_power_max]);
title('Power Distribution', 'FontSize', titleFontSize, 'FontWeight', 'bold');
xlabel('Subcarrier Number [#]', 'FontSize', axisFontSize);
ylabel('Power (dB)', 'FontSize', axisFontSize);
grid on;
text(0.4, 0.1, {"'X' = Subcarrier Getting", "Power Increment"}, 'Units', 'normalized', ...
     'FontSize', 9, 'Color', colors{2});
set(ax3, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

% Tile 4: Distortion
ax4 = nexttile;
p4 = plot(1:N, ND2time(1,:), 'o', 'MarkerSize', markerSize, 'Color', colors{1}, ...
    'MarkerFaceColor', colors{1}, 'MarkerEdgeColor', 'none', 'LineWidth', lineWidth);
hold on;
p4_mean = plot(1:N, ones(1,N)*mean(N0(10:end-10)), '-', 'Color', colors{3}, 'LineWidth', 2.5);
hold off;
xlim([1, N]);
ylim([0 max([ND2time(end,:) N0(10:end-10)])]);
legend({'ND2(n)','N_0'}, 'FontSize', fontSize-1, 'Location', 'best', 'Box', 'off');
title('Distortion', 'FontSize', titleFontSize, 'FontWeight', 'bold');
xlabel('Subcarrier Number [#]', 'FontSize', axisFontSize);
ylabel('Level [dB]', 'FontSize', axisFontSize);
grid on;
set(ax4, 'FontSize', fontSize, 'Box', 'on', 'LineWidth', 1.2);

%% Generate animation with audio
% Calculate frame rate
total_frames = floor(size(PCtime,1)/4);
frame_rate = total_frames / animation_length_seconds;

% Initialize audio array
audio_buffer = [];
audio_samples_per_frame = ceil(sample_rate / frame_rate);

% Create temporary audio file
temp_audio_file = 'temp_audio.wav';
if exist(temp_audio_file, 'file')
    delete(temp_audio_file);
end
audiowrite(temp_audio_file, 0, sample_rate);

% Initialize video file
video_filename = 'power_allocation_animation.mp4';
v = VideoWriter(video_filename, 'MPEG-4');
v.FrameRate = frame_rate;
v.Quality = 100;
open(v);

% First pass: generate audio track
fprintf('Generating audio track...\n');
for i = 1:4:size(PCtime,1)
    [~, ix] = max(Ctime(i,:));
    
    % Create silent audio frame by default
    frame_audio = zeros(1, audio_samples_per_frame);
    
    % Add sound effect when power is added to a different subcarrier
    if enable_sound && (i > 1) && (ix ~= prev_ix)
        % Calculate frequency based on subcarrier position and power level
        position_factor = (ix - 1) / (N - 1);
        power_at_x = 10*log10(PCtime(i,ix));
        power_normalized = min(max(power_at_x / 20, 0), 1);
        
        % Calculate tone frequency
        frequency = base_freq + freq_range * (0.7*position_factor + 0.3*power_normalized);
        
        % Create tone with envelope
        tone_samples = round(tone_duration * sample_rate);
        tone_time = (0:tone_samples-1) / sample_rate;
        envelope = sin(pi * tone_time / tone_duration).^0.5;
        tone = sin(2 * pi * frequency * tone_time) .* envelope * volume_scale;
        
        % Add tone to frame audio
        if length(tone) <= length(frame_audio)
            start_idx = floor((length(frame_audio) - length(tone)) / 2) + 1;
            frame_audio(start_idx:start_idx+length(tone)-1) = tone;
        else
            frame_audio = [tone(1:length(frame_audio)) zeros(1, max(0, length(frame_audio) - length(tone)))];
        end
    end
    
    % Append to audio buffer
    audio_buffer = [audio_buffer frame_audio];
    prev_ix = ix;
end

% Save audio to file
fprintf('Saving audio track...\n');
audiowrite(temp_audio_file, audio_buffer, sample_rate);

% Reset for animation
prev_ix = 0;

% Second pass: generate video frames
fprintf('Generating video...\n');
for i = 1:4:size(PCtime,1)
    % Update capacity plot
    set(p1, 'YData', Ctime(i,:));
    [~, ix] = max(Ctime(i,:));
    set(p1_max, 'XData', ix, 'YData', Ctime(i,ix));
    
    % Update capacity y-limits and ticks
    y_mean = mean(Ctime(i,:));
    y_range = 2*max(abs(y_mean-minmax(Ctime(i,:))));
    set(cap_ax, 'YLim', [-1 +1]*(y_range/2) + y_mean, ...
                'YTick', y_mean, ...
                'YTickLabel', sprintf('%.2f', y_mean));
    
    % Update SNR plot
    snr_data = 10*log10(SNRprojtime(i,:));
    set(p2, 'YData', snr_data);
    
    % Update SNR y-limits smoothly
    snr_min_target = min(snr_data);
    snr_max_target = max(snr_data);
    snr_padding = 0.1 * (snr_max_target - snr_min_target);
    snr_min_target = max(10, snr_min_target - snr_padding);
    snr_max_target = snr_max_target + snr_padding;
    
    current_snr_min = smoothing_factor_snr * current_snr_min + (1 - smoothing_factor_snr) * snr_min_target;
    current_snr_max = smoothing_factor_snr * current_snr_max + (1 - smoothing_factor_snr) * snr_max_target;
    set(ax2, 'YLim', [current_snr_min, current_snr_max]);
    
    % Update power plot
    power_data = 10*log10(PCtime(i,:));
    set(p3, 'YData', power_data);
    set(p3_max, 'XData', ix, 'YData', 10*log10(PCtime(i,ix)));
    
    % Check if power point needs immediate adjustment
    curr_ylim = get(ax3, 'YLim');
    power_at_x = 10*log10(PCtime(i,ix));
    needs_immediate_adjustment = (power_at_x > (curr_ylim(2) - 1));
    
    % Update power y-limits smoothly
    power_min_target = min(power_data);
    power_max_target = min(20, max(power_data));
    
    if (power_max_target - power_min_target) < min_power_range
        power_center = (power_max_target + power_min_target) / 2;
        power_min_target = power_center - min_power_range/2;
        power_max_target = power_center + min_power_range/2;
    end
    
    power_padding = 0.1 * (power_max_target - power_min_target);
    power_min_target = max(0, power_min_target - power_padding);
    power_max_target = min(20, power_max_target + power_padding);
    
    if needs_immediate_adjustment
        local_smoothing = 0.3;
    else
        local_smoothing = smoothing_factor_power;
    end
    current_power_min = local_smoothing * current_power_min + (1 - local_smoothing) * power_min_target;
    current_power_max = local_smoothing * current_power_max + (1 - local_smoothing) * power_max_target;
    
    if (current_power_max - current_power_min) < min_power_range
        power_center = (current_power_max + current_power_min) / 2;
        current_power_min = max(0, power_center - min_power_range/2);
        current_power_max = min(20, power_center + min_power_range/2);
    end
    
    set(ax3, 'YLim', [current_power_min, current_power_max]);
    
    % Update distortion plot
    set(p4, 'YData', ND2time(i,:));
    
    % Ensure consistent x-limits
    set(ax1, 'XLim', [1, N]);
    set(ax2, 'XLim', [1, N]);
    set(ax3, 'XLim', [1, N]);
    set(ax4, 'XLim', [1, N]);
    
    % Update display and capture frame
    drawnow;
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close video file
close(v);

%% Combine video with audio using FFmpeg
fprintf('Combining video with audio...\n');

% FFmpeg command for video-audio combination
if ispc
    ffmpeg_cmd = sprintf('ffmpeg -y -i %s -i %s -c:v copy -c:a aac -map 0:v:0 -map 1:a:0 %s_with_audio.mp4', ...
        video_filename, temp_audio_file, video_filename(1:end-4));
else
    ffmpeg_cmd = sprintf('ffmpeg -y -i %s -i %s -c:v copy -c:a aac -map 0:v:0 -map 1:a:0 %s_with_audio.mp4', ...
        video_filename, temp_audio_file, video_filename(1:end-4));
end

% Execute FFmpeg
status = system(ffmpeg_cmd);

if status == 0
    % Clean up and rename files
    delete(temp_audio_file);
    delete(video_filename);
    movefile(sprintf('%s_with_audio.mp4', video_filename(1:end-4)), video_filename);
    fprintf('Successfully created animation with audio: %s\n', video_filename);
else
    fprintf('WARNING: Could not add audio to video. FFmpeg may not be installed.\n');
    fprintf('Video without audio is available at: %s\n', video_filename);
    fprintf('Audio is available at: %s\n', temp_audio_file);
end