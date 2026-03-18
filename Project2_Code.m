clear all; close all; clc;
%Part 1
% Electrical Properties
Re = 0.350;             
Ri = 0.110;            
                       
Cm = 2.5;             
% Geometric Properties
axon_radius = 7e-4;     
axon_diameter = 2 * axon_radius;  
node_width = 1e-4;       
internodal_gap = 0.1;   
% Ionic Conductances
GNa = 1445;              
Gl = 128;                
% Nernst Potentials 
ENa = 115;             
El = -0.01;            
%  parameters
params.Re = Re;
params.Ri = Ri;
params.Cm = Cm;
params.d = axon_diameter;
params.L = node_width;
params.dx = internodal_gap;
params.GNa = GNa;
params.Gl = Gl;
params.ENa = ENa;
params.El = El;
%% Spatial Setup
% Spatial discretization (9 cm axon centered at x = 0)
axon_length = 9;       
x = -axon_length/2 : internodal_gap : axon_length/2;  
num_nodes = length(x);
dx = internodal_gap;
% Temporal discretization
dt = 0.001;              
t_end = 3;               
t = 0:dt:t_end;          
num_timesteps = length(t);
% Display simulation info
fprintf('=========================================================\n');
fprintf('   BME 471: CRRSS Action Potential Simulation\n');
fprintf('=========================================================\n');
fprintf('Axon length: %.1f cm\n', axon_length);
fprintf('Number of nodes: %d\n', num_nodes);
fprintf('Internodal gap: %.1f mm\n', internodal_gap * 10);
fprintf('Time step: %.4f ms\n', dt);
fprintf('Simulation duration: %.1f ms\n', t_end);
fprintf('=========================================================\n\n');
%% Part 2 Bidirectionel
fprintf('PART 2: Bidirectional Action Potential Simulation\n');
fprintf('---------------------------------------------------------\n');
% Stimulus timing parameters
t_stim_start = 0.5;      
t_anodic_dur = 0.8;      
t_cathodic_dur = 0.2;   
% Electrode distances and  threshold currents 
distances = [0.1, 0.2, 0.4, 0.8];  % cm (1, 2, 4, 8 mm)
distance_labels = {'1 mm', '2 mm', '4 mm', '8 mm'};
% Threshold currents
threshold_currents = [-209, -782, -3653, -21089];  % uA
anodic_amp = 100;
fprintf('\nMinimum Cathodic Currents for Bidirectional AP:\n');
for d_idx = 1:length(distances)
    fprintf('  %s: Cathodic = %d uA, Anodic = %d uA\n', ...
            distance_labels{d_idx}, threshold_currents(d_idx), anodic_amp);
end
%% Run simulation for 1mm case and create surface plot
fprintf('\nRunning simulation for 1 mm case...\n');
electrode_dist = 0.1;  % cm (1 mm)
cath_amp = threshold_currents(1);  % -209 uA
anod_amp = anodic_amp;  % 100 uA
% Generate stimulus waveform
I_stim = generate_biphasic_stimulus(t, t_stim_start, t_anodic_dur, ...
                                    t_cathodic_dur, anod_amp, cath_amp);
% Calculate extracellular potential from monopole
Ve = calculate_Ve_monopole(x, t, I_stim, electrode_dist, 0, Re);
% Run CRRSS simulation
[Vm, m_gate, h_gate] = run_CRRSS(x, t, Ve, dt, params);
% Create surface plot for bidirectional AP (1mm)
figure('Name', 'Bidirectional AP - 1mm', 'Position', [100, 100, 900, 700]);
surf(t, x, Vm);
shading flat;
colormap(jet);
colorbar;
xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Position on Axon (cm)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Membrane Voltage (mV)', 'FontSize', 12, 'FontWeight', 'bold');
title({'Bidirectional Action Potential Propagation with a monopole 1 mm from axon'}, ...
      'FontSize', 14);
xlim([0 t_end]);
ylim([min(x) max(x)]);
zlim([-20 100]);
view([-37.5, 30]);
saveas(gcf, 'Bidirectional_AP_1mm_surface.fig');
saveas(gcf, 'Bidirectional_AP_1mm_surface.png');
fprintf('Surface plot saved: Bidirectional_AP_1mm_surface.png\n');
%% Create AVI movie for bidirectional AP
fprintf('Creating AVI movie for bidirectional AP...\n');
create_AP_movie(t, x, Vm, 'Bidirectional_AP_1mm.avi', ...
                'Action Potential Bidirectional Propagation (1mm away from the axon)', ...
                [-20 100]);
fprintf('Movie saved: Bidirectional_AP_1mm.avi\n');
%% Create 2x2 subplot showing all four distances
fprintf('\nGenerating comparison plot for all distances...\n');
figure('Name', 'All Distances', 'Position', [50, 50, 1200, 900]);
for d_idx = 1:length(distances)
    electrode_dist = distances(d_idx);
    cath_amp = threshold_currents(d_idx);
    
    % Generate stimulus and run simulation
    I_stim_d = generate_biphasic_stimulus(t, t_stim_start, t_anodic_dur, ...
                                          t_cathodic_dur, anod_amp, cath_amp);
    Ve_d = calculate_Ve_monopole(x, t, I_stim_d, electrode_dist, 0, Re);
    [Vm_d, ~, ~] = run_CRRSS(x, t, Ve_d, dt, params);
    
    subplot(2, 2, d_idx);
    surf(t, x, Vm_d);
    shading flat;
    colormap(jet);
    xlabel('Time (ms)');
    ylabel('Position (cm)');
    zlabel('Vm (mV)');
    title(sprintf('%s: I_{cath} = %d uA', distance_labels{d_idx}, cath_amp));
    xlim([0 t_end]);
    zlim([-20 100]);
    view([-37.5, 30]);
end
sgtitle('Bidirectional AP at Different Electrode Distances', 'FontSize', 14, 'FontWeight', 'bold');
saveas(gcf, 'Bidirectional_AP_all_distances.fig');
saveas(gcf, 'Bidirectional_AP_all_distances.png');
fprintf('All distances plot saved: Bidirectional_AP_all_distances.png\n\n');
%% Part 3 Unidirectionional
fprintf('=========================================================\n');
fprintf('PART 3: Unidirectional Action Potential Generation\n');
fprintf('=========================================================\n');
% Extended simulation 
t_end_uni = 3.5;         % ms
t_uni = 0:dt:t_end_uni;
num_timesteps_uni = length(t_uni);
% Electrode configuration for unidirectional block
electrode_z = 0.1;        
cathode_x = 0;            
anode_x = 0.2;            
% Stimulus parameters
cath_amp_uni = -300;      
anod_amp_uni = 200;       
% Timing
t_stim_start_uni = 0.3;   
cath_duration = 0.15;     
anod_duration = 0.3;     
fprintf('\nElectrode configuration:\n');
fprintf('  Cathode position: x = %.1f cm (center)\n', cathode_x);
fprintf('  Anode position: x = %.1f cm\n', anode_x);
fprintf('  Distance from axon: %.1f mm\n', electrode_z * 10);
fprintf('  Cathode amplitude: %d uA, duration: %.2f ms\n', cath_amp_uni, cath_duration);
fprintf('  Anode amplitude: %d uA, duration: %.2f ms\n', anod_amp_uni, anod_duration);
% Generate stimulus currents for each electrode
I_cathode = zeros(1, num_timesteps_uni);
I_anode = zeros(1, num_timesteps_uni);
for i = 1:num_timesteps_uni
    % Cathode pulse
    if t_uni(i) >= t_stim_start_uni && t_uni(i) < t_stim_start_uni + cath_duration
        I_cathode(i) = cath_amp_uni;
    end
    % Anode pulse (starts at same time, lasts longer)
    if t_uni(i) >= t_stim_start_uni && t_uni(i) < t_stim_start_uni + anod_duration
        I_anode(i) = anod_amp_uni;
    end
end
% Calculate combined extracellular potential from both electrodes
Ve_uni = zeros(num_nodes, num_timesteps_uni);
for i = 1:num_nodes
    for j = 1:num_timesteps_uni
        % Distance from node to cathode
        r_cath = sqrt((x(i) - cathode_x)^2 + electrode_z^2);
        % Distance from node to anode
        r_anod = sqrt((x(i) - anode_x)^2 + electrode_z^2);
        
        % Superposition of potentials from both electrodes
        Ve_uni(i,j) = (Re * I_cathode(j)) / (4 * pi * r_cath) + ...
                      (Re * I_anode(j)) / (4 * pi * r_anod);
    end
end
% Run CRRSS simulation
[Vm_uni, m_uni, h_uni] = run_CRRSS(x, t_uni, Ve_uni, dt, params);
% Verify unidirectional propagation
left_max = max(max(Vm_uni(1:round(num_nodes/4), :)));
right_max = max(max(Vm_uni(round(3*num_nodes/4):end, :)));
fprintf('\nVerification:\n');
fprintf('  Max Vm on left side (x < -2.25 cm): %.1f mV\n', left_max);
fprintf('  Max Vm on right side (x > 2.25 cm): %.1f mV\n', right_max);
if left_max > 50 && right_max < 50
    fprintf('  Result: SUCCESS - Unidirectional propagation to the LEFT!\n');
elseif right_max > 50 && left_max < 50
    fprintf('  Result: SUCCESS - Unidirectional propagation to the RIGHT!\n');
else
    fprintf('  Note: Adjust parameters for better directional control\n');
end
% Create surface plot for unidirectional AP
figure('Name', 'Unidirectional AP', 'Position', [100, 100, 900, 700]);
surf(t_uni, x, Vm_uni);
shading flat;
colormap(jet);
colorbar;
xlabel('Time (ms)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('Position on Axon (cm)', 'FontSize', 12, 'FontWeight', 'bold');
zlabel('Membrane Voltage (mV)', 'FontSize', 12, 'FontWeight', 'bold');
title('Unidirectional Action Potential Propagation 1mm from axon', 'FontSize', 14);
xlim([0 t_end_uni]);
ylim([min(x) max(x)]);
zlim([-50 100]);
view([45, 30]);
saveas(gcf, 'Unidirectional_AP_surface.fig');
saveas(gcf, 'Unidirectional_AP_surface.png');
fprintf('\nSurface plot saved: Unidirectional_AP_surface.png\n');
% Create AVI movie for unidirectional AP
fprintf('Creating unidirectional AVI movie...\n');
create_AP_movie(t_uni, x, Vm_uni, 'Unidirectional_AP.avi', ...
                'Unidirectional Action Potential Propagation (1mm from axon)', [-50 100]);
fprintf('Movie saved: Unidirectional_AP.avi\n');
%% Create stimulus waveform figure
figure('Name', 'Stimulus Analysis', 'Position', [100, 100, 1000, 400]);
subplot(1, 2, 1);
plot(t_uni, I_cathode, 'b-', 'LineWidth', 2);
hold on;
plot(t_uni, I_anode, 'r-', 'LineWidth', 2);
hold off;
xlabel('Time (ms)');
ylabel('Current (uA)');
title('Unidirectional Stimulus Waveforms');
legend('Cathode', 'Anode', 'Location', 'best');
grid on;
xlim([0 1]);
subplot(1, 2, 2);
% Plot voltage traces at different positions
positions = [1, round(num_nodes/4), round(num_nodes/2), round(3*num_nodes/4), num_nodes];
colors = lines(length(positions));
for p = 1:length(positions)
    plot(t_uni, Vm_uni(positions(p), :), 'Color', colors(p,:), 'LineWidth', 1.5);
    hold on;
end
hold off;
xlabel('Time (ms)');
ylabel('Membrane Voltage (mV)');
title('Vm at Different Positions');
legend_labels = cell(1, length(positions));
for p = 1:length(positions)
    legend_labels{p} = sprintf('x = %.1f cm', x(positions(p)));
end
legend(legend_labels, 'Location', 'best');
grid on;
xlim([0 t_end_uni]);
saveas(gcf, 'Stimulus_Analysis.fig');
saveas(gcf, 'Stimulus_Analysis.png');
%% Final Summary
fprintf('\n=========================================================\n');
fprintf('   SIMULATION COMPLETE!\n');
fprintf('=========================================================\n');
fprintf('Output files generated:\n');
fprintf('  - Bidirectional_AP_1mm_surface.png/.fig\n');
fprintf('  - Bidirectional_AP_1mm.avi\n');
fprintf('  - Bidirectional_AP_all_distances.png/.fig\n');
fprintf('  - Unidirectional_AP_surface.png/.fig\n');
fprintf('  - Unidirectional_AP.avi\n');
fprintf('  - Stimulus_Analysis.png/.fig\n');
fprintf('=========================================================\n');
%% Functions
function I_stim = generate_biphasic_stimulus(t, t_start, t_anod_dur, t_cath_dur, anod_amp, cath_amp)
    I_stim = zeros(1, length(t));
    
    for i = 1:length(t)
        if t(i) >= t_start && t(i) < t_start + t_anod_dur
            I_stim(i) = anod_amp;
        end
        if t(i) >= t_start + t_anod_dur && t(i) < t_start + t_anod_dur + t_cath_dur
            I_stim(i) = cath_amp;
        end
    end
end
function Ve = calculate_Ve_monopole(x, t, I_stim, z, x_electrode, Re)
    
    num_nodes = length(x);
    num_timesteps = length(t);
    Ve = zeros(num_nodes, num_timesteps);
    
    for i = 1:num_nodes
        r = sqrt((x(i) - x_electrode)^2 + z^2);
        
        Ve(i, :) = (Re * I_stim) / (4 * pi * r);
    end
end
function [Vm, m, h] = run_CRRSS(x, t, Ve, dt, params)
    % Extract parameters
    d = params.d;       % axon diameter
    L = params.L;       % node width
    dx = params.dx;     % internodal distance
    Ri = params.Ri;     % intracellular resistance
    Cm = params.Cm;     % membrane capacitance
    GNa = params.GNa;   % sodium conductance
    Gl = params.Gl;     % leak conductance
    ENa = params.ENa;   % sodium Nernst potential
    El = params.El;     % leak Nernst potential
    
    num_nodes = length(x);
    num_timesteps = length(t);
    
    % Initialize state variables
    Vm = zeros(num_nodes, num_timesteps);
    m = zeros(num_nodes, num_timesteps);
    h = zeros(num_nodes, num_timesteps);
    
    % Set initial conditions for gating variables at rest (Vm = 0)
    [alpha_m0, beta_m0] = mgate(0);
    [alpha_h0, beta_h0] = hgate(0);
    m_inf = alpha_m0 / (alpha_m0 + beta_m0);
    h_inf = alpha_h0 / (alpha_h0 + beta_h0);
    m(:, 1) = m_inf;
    h(:, 1) = h_inf;
    
    % Main simulation loop
    for j = 1:(num_timesteps - 1)
        for i = 2:(num_nodes - 1)
            [alpha_m, beta_m] = mgate(Vm(i, j));
            [alpha_h, beta_h] = hgate(Vm(i, j));
            
            dmdt = -(alpha_m + beta_m) * m(i,j) + alpha_m;
            m(i, j+1) = m(i, j) + dmdt * dt;
            
            dhdt = -(alpha_h + beta_h) * h(i,j) + alpha_h;
            h(i, j+1) = h(i, j) + dhdt * dt;
            
            m(i, j+1) = max(0, min(1, m(i, j+1)));
            h(i, j+1) = max(0, min(1, h(i, j+1)));
            
            INa = GNa * m(i,j)^3 * h(i,j) * (Vm(i,j) - ENa);  % Sodium current
            Il = Gl * (Vm(i,j) - El);                          % Leak current
            Im = INa + Il;                                      % Total membrane current
            
            d2Ve = (Ve(i-1,j) - 2*Ve(i,j) + Ve(i+1,j)) / (dx^2);
            d2Vm = (Vm(i-1,j) - 2*Vm(i,j) + Vm(i+1,j)) / (dx^2);
            
        
            If = (d * dx) / (4 * Ri * L) * (d2Ve + d2Vm);
            
            dVmdt = (If - Im) / Cm;
            Vm(i, j+1) = Vm(i, j) + dVmdt * dt;
        end
        
        Vm(1, j+1) = Vm(2, j+1);
        Vm(num_nodes, j+1) = Vm(num_nodes-1, j+1);
        m(1, j+1) = m(2, j+1);
        m(num_nodes, j+1) = m(num_nodes-1, j+1);
        h(1, j+1) = h(2, j+1);
        h(num_nodes, j+1) = h(num_nodes-1, j+1);
    end
end
function [alpha_m, beta_m] = mgate(Vm)
    
    alpha_m = (97 + 0.363 * Vm) / (1 + exp((31 - Vm) / 5.3));
    beta_m = alpha_m / exp((Vm - 23.8) / 4.17);
end
function [alpha_h, beta_h] = hgate(Vm)
   
    
    beta_h = 15.6 / (1 + exp((24 - Vm) / 10));
    alpha_h = beta_h / exp((Vm - 5.5) / 5);
end
function create_AP_movie(t, x, Vm, filename, plot_title, zlimits)
 
    video = VideoWriter(filename, 'Motion JPEG AVI');
    video.FrameRate = 30;
    video.Quality = 90;
    open(video);
    
    fig = figure('Position', [100, 100, 800, 600], 'Visible', 'on');
    
    num_frames = min(300, length(t));
    frame_indices = round(linspace(2, length(t), num_frames));  % Start from 2 to ensure matrix
    
    for idx = 1:length(frame_indices)
        j = frame_indices(idx);
        
        if j >= 2
            surf(t(1:j), x, Vm(:, 1:j));
            shading flat;
            colormap(jet);
            xlabel('Time (ms)', 'FontSize', 11);
            ylabel('Position on Axon (cm)', 'FontSize', 11);
            zlabel('Membrane Voltage (mV)', 'FontSize', 11);
            title(plot_title, 'FontSize', 12);
            xlim([0 max(t)]);
            ylim([min(x) max(x)]);
            zlim(zlimits);
            view([-37.5, 30]);
            drawnow;
            
            frame = getframe(fig);
            writeVideo(video, frame);
        end
    end
    
    for extra = 1:15
        frame = getframe(fig);
        writeVideo(video, frame);
    end
    
    close(video);
    close(fig);
end