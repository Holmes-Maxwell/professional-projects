% hh_main.m
% Main script for BME 471 Homework #1 — Hodgkin-Huxley using ODE45
clear; close all; clc;

%% -----------------------
% Section A: Parameters
% ------------------------
params.gNa_max = 120;    % mS/cm^2
params.gK_max  = 36;     % mS/cm^2
params.gL      = 0.3;    % mS/cm^2
params.Cm      = 1.0;    % uF/cm^2

% Concentrations (given) -> compute Nernst using 61.54*log10([out]/[in]) at 37C
RFa = 61.54; % mV per decade at 37 C
Na_out = 490; Na_in = 50;
K_out  = 20;  K_in  = 400;
params.ENa = RFa * log10(Na_out/Na_in);   % mV
params.EK  = RFa * log10(K_out/K_in);     % mV
params.EL  = -50;                          % Leak Nernst potential (given)

% Display reversal potentials
fprintf('ENa = %.2f mV, EK = %.2f mV, EL = %.2f mV\n', params.ENa, params.EK, params.EL);

%% -----------------------
% Section B: initial conditions & ODE options
% ------------------------
V0 = -65;    % mV (typical resting)
% compute steady-state gating variables at V0 using functions in hh_diff_eq (we call helper)
[~, infs] = init_gates(V0);
y0 = [V0; infs.n; infs.m; infs.h];

tspan = [0 20]; % ms, adjust later for longer runs
options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',0.01);

%% -----------------------
% Section C: Stimulus definitions & helper functions
% ------------------------
% Stimulus shape function handle: Is = stim(t, amplitude, tstart, dur)
stim = @(t, amp, tstart, dur) amp .* (t>=tstart & t < tstart+dur);  % uA/cm^2

%% -----------------------
% Section D: Find threshold amplitude for depolarizing square pulse (0.35 ms)
% ------------------------
pulse_dur = 0.35;   % ms (given)
tstart = 1.0;       % ms, start pulse at 1 ms to allow rest
% we will sweep amplitudes from small to larger until an AP occurs
amps = 0:1:200; % uA/cm^2 search range (coarse -> step 1)
threshold_amp = NaN;

for A = amps
    % create stimulus function for this amplitude
    params.stimfun = @(t) stim(t,A,tstart,pulse_dur);
    % integrate
    [t,y] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params), tspan, y0, options);
    % detect spike: peak > 0 mV (conservative)
    if max(y(:,1)) > 0
        threshold_amp = A;
        break;
    end
end

if isnan(threshold_amp)
    warning('No action potential found in amplitude range. Try larger amplitudes or increase tspan.');
else
    fprintf('Found threshold depolarizing amplitude = %g uA/cm^2 (pulse dur = %.3g ms)\n', threshold_amp, pulse_dur);
end

% Save the found amplitude for later use
A_thresh = threshold_amp;

%% -----------------------
% Section E: Plot a typical action potential with gating & currents (use threshold amp)
% ------------------------
if ~isnan(A_thresh)
    params.stimfun = @(t) stim(t,A_thresh,tstart,pulse_dur);
    tspan_long = [0 20];
    [t_ap,y_ap] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params), tspan_long, y0, options);
    % compute currents and stimulus at solution times
    Isq = arrayfun(params.stimfun, t_ap);
    [INaq, IKq, ILq] = currents_from_solution(y_ap, params); % helper below
    
    % Main AP figure
    figure('Name','Typical Action Potential','NumberTitle','off','Position',[100 100 1400 900]);
    
    % Stimulus and Membrane Potential
    subplot(3,3,1);
    plot(t_ap, Isq, 'k', 'LineWidth', 1.5);
    xlabel('Time (ms)'); ylabel('I_{stim} (\muA/cm^2)');
    title('Stimulus Current'); grid on;
    
    subplot(3,3,2);
    plot(t_ap, y_ap(:,1), 'b', 'LineWidth', 2);
    xlabel('Time (ms)'); ylabel('V_m (mV)');
    title('Membrane Potential'); grid on;
    
    % Gating variables - all on one plot
    subplot(3,3,3);
    plot(t_ap, y_ap(:,2), 'r', 'LineWidth', 2); hold on;
    plot(t_ap, y_ap(:,3), 'g', 'LineWidth', 2);
    plot(t_ap, y_ap(:,4), 'm', 'LineWidth', 2);
    xlabel('Time (ms)'); ylabel('Gating Variable');
    title('Gating Variables');
    legend('n (K^+)', 'm (Na^+)', 'h (Na^+)', 'Location', 'best'); grid on;
    
    % All currents on one plot
    subplot(3,3,7);
    plot(t_ap, INaq, 'r', 'LineWidth', 2); hold on;
    plot(t_ap, IKq, 'b', 'LineWidth', 2);
    plot(t_ap, ILq, 'g', 'LineWidth', 2);
    xlabel('Time (ms)'); ylabel('Current (\muA/cm^2)');
    title('Ionic Currents');
    legend('I_{Na}', 'I_{K}', 'I_{L}', 'Location', 'best'); grid on;
end

%% -----------------------
% Section F: REFRACTORY PERIOD DEBUGGING
% Let's figure out what's actually happening step by step
% -----------------------

fprintf('\n=== COMPREHENSIVE REFRACTORY PERIOD DEBUGGING ===\n');

% First, let's verify our threshold amplitude is correct
fprintf('1. VERIFYING THRESHOLD AMPLITUDE:\n');
fprintf('   Current threshold: %.1f uA/cm²\n', A_thresh);

% Test if this really produces an AP
params_test = params;
params_test.stimfun = @(t) stim(t, A_thresh, 1.0, pulse_dur);
[t_single, y_single] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_test), [0 10], y0, options);
max_V_single = max(y_single(:,1));
fprintf('   Single pulse test: max V = %.1f mV\n', max_V_single);
if max_V_single > 0
    fprintf('   ✓ Single AP confirmed\n');
else
    fprintf('   ⚠ NO AP with threshold amplitude! This is the problem!\n');
    % Let's find the real threshold
    fprintf('   Finding actual threshold...\n');
    for A_test = A_thresh:1:50
        params_test.stimfun = @(t) stim(t, A_test, 1.0, pulse_dur);
        [~, y_test] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_test), [0 10], y0, options);
        if max(y_test(:,1)) > 0
            fprintf('   Real threshold: %.1f uA/cm²\n', A_test);
            A_thresh = A_test; % Update threshold
            break;
        end
    end
end

% 2. Test what happens with two pulses at different delays
fprintf('\n2. TESTING TWO-PULSE RESPONSE:\n');
test_delays = [1.0, 2.0, 3.0, 5.0, 8.0, 12.0, 15.0];

figure('Name','Refractory Period Debugging','Color','w','Position',[100 100 1400 900]);

for i = 1:length(test_delays)
    delay = test_delays(i);
    
    % Test with threshold amplitude
    params_local = params;
    params_local.stimfun = @(t) stim(t, A_thresh, 1.0, pulse_dur) + ...
                                 stim(t, A_thresh, 1.0 + delay, pulse_dur);
    [t_debug, y_debug] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_local), [0 25], y0, options);
    
    subplot(3,3,i);
    plot(t_debug, y_debug(:,1), 'b-', 'LineWidth', 2);
    hold on;
    
    % Mark stimulus periods
    yl = ylim;
    plot([1.0, 1.0], yl, 'r--', 'LineWidth', 1);
    plot([1.0+pulse_dur, 1.0+pulse_dur], yl, 'r--', 'LineWidth', 1);
    plot([1.0+delay, 1.0+delay], yl, 'g--', 'LineWidth', 1);
    plot([1.0+delay+pulse_dur, 1.0+delay+pulse_dur], yl, 'g--', 'LineWidth', 1);
    
    % Analyze the response
    second_stim_start = 1.0 + delay;
    search_after_second = (t_debug >= second_stim_start + 1.0) & (t_debug <= second_stim_start + 10);
    max_V_after_second = max(y_debug(search_after_second, 1));
    
    % Check gating variables at time of second stimulus
    [~, idx] = min(abs(t_debug - second_stim_start));
    h_val = y_debug(idx, 4); % h-gate (sodium inactivation)
    n_val = y_debug(idx, 2); % n-gate (potassium activation)
    
    title(sprintf('Delay %.1f ms\nMax after: %.1f mV\nh=%.3f, n=%.3f', ...
          delay, max_V_after_second, h_val, n_val));
    xlabel('Time (ms)');
    ylabel('V_m (mV)');
    grid on;
    
    fprintf('   Delay %.1f ms: max after 2nd = %.1f mV, h=%.3f, n=%.3f\n', ...
            delay, max_V_after_second, h_val, n_val);
end

% 3. Let's check the gating variable dynamics more carefully
fprintf('\n3. GATING VARIABLE ANALYSIS:\n');
fprintf('   At rest (V = %.1f mV): n=%.3f, m=%.3f, h=%.3f\n', V0, infs.n, infs.m, infs.h);

% Test what happens during and after an AP
params_test = params;
params_test.stimfun = @(t) stim(t, A_thresh, 1.0, pulse_dur);
[t_gating, y_gating] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_test), [0 20], y0, options);

% Find when h-gate recovers
h_recovery_threshold = 0.5; % Typical recovery level
recovery_time = NaN;
for i = 1:length(t_gating)
    if t_gating(i) > 5 && y_gating(i,4) > h_recovery_threshold
        recovery_time = t_gating(i);
        break;
    end
end

fprintf('   h-gate recovery (>0.5) at: %.1f ms\n', recovery_time);

% 4. Simple refractory period test with clear criteria
fprintf('\n4. SIMPLE REFRACTORY PERIOD TEST:\n');

delays_to_test = 1:1:20;
second_ap_found = false(size(delays_to_test));

for i = 1:length(delays_to_test)
    delay = delays_to_test(i);
    
    params_local = params;
    params_local.stimfun = @(t) stim(t, A_thresh, 1.0, pulse_dur) + ...
                                 stim(t, A_thresh, 1.0 + delay, pulse_dur);
    [t_simple, y_simple] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_local), [0 30], y0, options);
    
    % Look for two clear AP peaks
    V_simple = y_simple(:,1);
    ap_peaks = 0;
    last_peak_time = -10;
    
    for j = 3:length(V_simple)-2
        if V_simple(j) > V_simple(j-1) && V_simple(j) > V_simple(j+1) && ...
           V_simple(j) > V_simple(j-2) && V_simple(j) > V_simple(j+2)
            if V_simple(j) > 0 && (t_simple(j) - last_peak_time) > 2.0
                ap_peaks = ap_peaks + 1;
                last_peak_time = t_simple(j);
            end
        end
    end
    
    second_ap_found(i) = (ap_peaks >= 2);
    
    if second_ap_found(i)
        fprintf('   ✓ Second AP found at delay %.1f ms (%d peaks total)\n', delay, ap_peaks);
        min_gap_simple = delay;
        break;
    end
end

if any(second_ap_found)
    fprintf('   Minimum gap for second AP: %.1f ms\n', min_gap_simple);
else
    fprintf('   ⚠ No second AP found at any delay up to %.1f ms\n', delays_to_test(end));
    
    % Let's try with stronger stimulus
    fprintf('\n5. TESTING WITH STRONGER STIMULUS:\n');
    strong_amp = 2 * A_thresh;
    for delay = [3.0, 5.0, 8.0]
        params_local = params;
        params_local.stimfun = @(t) stim(t, strong_amp, 1.0, pulse_dur) + ...
                                     stim(t, strong_amp, 1.0 + delay, pulse_dur);
        [t_strong, y_strong] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_local), [0 30], y0, options);
        
        max_V_strong = max(y_strong(:,1));
        fprintf('   Strong amp %.1f uA/cm², delay %.1f ms: max V = %.1f mV\n', ...
                strong_amp, delay, max_V_strong);
    end
end

% 6. Final analysis with the working approach
fprintf('\n6. FINAL REFRACTORY PERIOD ANALYSIS:\n');

% Use the simple peak counting method
function min_gap = find_refractory_period(amp, params, y0, pulse_dur)
    delays = 1:0.5:20;
    for i = 1:length(delays)
        delay = delays(i);
        params_local = params;
        params_local.stimfun = @(t) stim(t, amp, 1.0, pulse_dur) + ...
                                     stim(t, amp, 1.0 + delay, pulse_dur);
        [t_test, y_test] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_local), [0 30], y0, options);
        
        % Count AP peaks
        V_test = y_test(:,1);
        peaks = 0;
        for j = 3:length(V_test)-2
            if V_test(j) > 20 && V_test(j) > V_test(j-1) && V_test(j) > V_test(j+1) && ...
               V_test(j) > V_test(j-2) && V_test(j) > V_test(j+2)
                peaks = peaks + 1;
            end
        end
        
        if peaks >= 2
            min_gap = delay;
            return;
        end
    end
    min_gap = NaN;
end

% Test both amplitudes
min_gap_same = find_refractory_period(A_thresh, params, y0, pulse_dur);
min_gap_double = find_refractory_period(2*A_thresh, params, y0, pulse_dur);

fprintf('   Same amplitude (%.1f uA/cm²): ', A_thresh);
if ~isnan(min_gap_same)
    fprintf('refractory period = %.1f ms\n', min_gap_same);
else
    fprintf('no second AP found\n');
end

fprintf('   Double amplitude (%.1f uA/cm²): ', 2*A_thresh);
if ~isnan(min_gap_double)
    fprintf('refractory period = %.1f ms\n', min_gap_double);
else
    fprintf('no second AP found\n');
end

%% Create final summary plot
figure('Name','Final Refractory Analysis','Color','w','Position',[200 200 1000 600]);

if ~isnan(min_gap_same)
    % Plot the successful case
    params_final = params;
    params_final.stimfun = @(t) stim(t, A_thresh, 1.0, pulse_dur) + ...
                                 stim(t, A_thresh, 1.0 + min_gap_same, pulse_dur);
    [t_final, y_final] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_final), [0 25], y0, options);
    I_final = arrayfun(params_final.stimfun, t_final);
    
    subplot(1,2,1);
    plot(t_final, I_final, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_final, y_final(:,1), 'r-', 'LineWidth', 2);
    xlabel('Time (ms)');
    title(sprintf('Same Amp: Gap = %.1f ms', min_gap_same));
    legend('Stimulus', 'V_m', 'Location', 'best');
    grid on;
end

if ~isnan(min_gap_double)
    % Plot the successful case for double amplitude
    params_final = params;
    params_final.stimfun = @(t) stim(t, 2*A_thresh, 1.0, pulse_dur) + ...
                                 stim(t, 2*A_thresh, 1.0 + min_gap_double, pulse_dur);
    [t_final2, y_final2] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_final), [0 25], y0, options);
    I_final2 = arrayfun(params_final.stimfun, t_final2);
    
    subplot(1,2,2);
    plot(t_final2, I_final2, 'b-', 'LineWidth', 2);
    hold on;
    plot(t_final2, y_final2(:,1), 'r-', 'LineWidth', 2);
    xlabel('Time (ms)');
    title(sprintf('Double Amp: Gap = %.1f ms', min_gap_double));
    legend('Stimulus', 'V_m', 'Location', 'best');
    grid on;
end

fprintf('\n=== DEBUGGING COMPLETE ===\n');
%% -----------------------
% Section G: Anode-break analysis (part 4)
% Using your function-based approach
% -----------------------

% Create anode_break_search function
function [minDur, output] = anode_break_search(ode_func, params, y0, tstart, amp, tspan, maxStep)
    % Search for minimum hyperpolarization duration that produces anode-break spike
    dur_vec = 0.1:0.05:20; % ms - search from 0.1 to 20ms
    
    options = odeset('RelTol',1e-4,'AbsTol',[1e-8 1e-8 1e-8 1e-8],'MaxStep',maxStep);
    minDur = NaN;
    output = [];
    
    for i = 1:length(dur_vec)
        dur_test = dur_vec(i);
        
        % Create hyperpolarizing pulse
        stim_fun = @(t) (t >= tstart & t < tstart + dur_test) * (-amp); % Negative current
        params_local = params;
        params_local.stimfun = stim_fun;
        
        % Solve ODE
        [t, y] = ode45(@(tt,yy) ode_func(tt,yy,params_local), tspan, y0, options);
        
        % Detect if anode-break spike occurred after pulse ends
        pulse_end = tstart + dur_test;
        search_mask = (t >= pulse_end + 0.5) & (t <= pulse_end + 15);
        if sum(search_mask) > 0
            V_after_pulse = y(search_mask, 1);
            if max(V_after_pulse) > 20 % Anode-break spike detected
                minDur = dur_test;
                
                % Store output for plotting
                output.t = t;
                output.V = y(:,1);
                output.n = y(:,2);
                output.m = y(:,3);
                output.h = y(:,4);
                output.I = arrayfun(stim_fun, t);
                break;
            end
        end
    end
end

% Now run anode-break analysis
tspanAB = [0 50]; % ms - longer time span for anode-break
[minDurAB, out4] = anode_break_search(@hh_diff_eq, params, y0, tstart_ref, amp_test, tspanAB, maxStep);

if isnan(minDurAB)
    warning('Problem 4: No anode-break spike found up to 20 ms hyperpolarization (|A|=%.3f).', amp_test);
else
    fprintf('Problem 4: Minimal hyperpolarizing duration for anode-break (|A|=%.3f) = %.3f ms\n', ...
            amp_test, minDurAB);
end

figure('Name','P4 Anode-break','Color','w');
subplot(3,1,1); 
plot(out4.t, out4.I, 'LineWidth',1.5); 
ylabel('I_{stim} (\muA/cm^2)');
title(sprintf('Problem 4: Anode-break (dur=%.3f ms, amp=-%.3f)', minDurAB, amp_test));
grid on;

subplot(3,1,2); 
plot(out4.t, out4.V, 'LineWidth',1.5); 
ylabel('V_m (mV)');
grid on;

subplot(3,1,3); 
plot(out4.t, [out4.n, out4.m, out4.h], 'LineWidth',1.2);
legend('n','m','h','Location','best'); 
xlabel('Time (ms)'); 
ylabel('Gates');
grid on;
%% -----------------------
% Section G: Anode-break (hyperpolarizing) test (part 4)
% Use same magnitude but negative (hyperpolarizing) and increase duration until
% anode-break spike occurs upon release.
% ------------------------
A_hyp = -A_thresh; % negative current with same magnitude
dur_vec = 0.05:0.05:20; % ms durations to test
anode_found = false;
anode_dur = NaN;
for i = 1:length(dur_vec)
    dur_test = dur_vec(i);
    params_local = params;
    % stimulus: hyperpolarizing pulse starting at 1 ms, duration dur_test
    params_local.stimfun = @(t) stim(t, A_hyp, 1.0, dur_test);
    tspan_test = [0 50];
    [t_test,y_test] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params_local), tspan_test, y0, options);
    % check for spike after release (after 1.0+dur_test)
    t_after = t_test >= (1.0 + dur_test);
    V_after = y_test(t_after,1);
    if max(V_after) > 0
        anode_found = true;
        anode_dur = dur_test;
        break;
    end
end

if anode_found
    fprintf('Anode-break spike found for hyperpolarizing amplitude %g uA/cm^2 at duration = %.3g ms\n', A_hyp, anode_dur);
else
    fprintf('No anode-break spike found up to %.2f ms duration for hyperpolarizing amplitude %g\n', dur_vec(end), A_hyp);
end

% For visualization, run one case with the anode_dur (if found) and plot
if anode_found
    params.stimfun = @(t) stim(t, A_hyp, 1.0, anode_dur);
    tspan_vis = [0 50];
    [t_ab,y_ab] = ode45(@(tt,yy) hh_diff_eq(tt,yy,params), tspan_vis, y0, options);
    Is_ab = arrayfun(params.stimfun, t_ab);
    [INa_ab, IK_ab, IL_ab] = currents_from_solution(y_ab, params);
    
    figure('Name','Anode-Break Excitation','NumberTitle','off','Position',[150 150 1200 800]);
    
    subplot(2,2,1);
    plot(t_ab, Is_ab, 'k', 'LineWidth', 2);
    ylabel('I_{stim} (\muA/cm^2)');
    title('Hyperpolarizing Stimulus Pulse');
    grid on;
    
    subplot(2,2,2);
    plot(t_ab, y_ab(:,1), 'b', 'LineWidth', 2);
    ylabel('V_m (mV)');
    title('Membrane Potential (Anode-Break)');
    grid on;
    
    subplot(2,2,3);
    plot(t_ab, INa_ab, 'r', 'LineWidth', 2); hold on;
    plot(t_ab, IK_ab, 'b', 'LineWidth', 2);
    plot(t_ab, IL_ab, 'g', 'LineWidth', 2);
    ylabel('Current (\muA/cm^2)');
    title('Ionic Currents');
    legend('I_{Na}', 'I_{K}', 'I_{L}', 'Location', 'best');
    grid on;
    
    subplot(2,2,4);
    plot(t_ab, y_ab(:,2), 'r', 'LineWidth', 2); hold on;
    plot(t_ab, y_ab(:,3), 'g', 'LineWidth', 2);
    plot(t_ab, y_ab(:,4), 'm', 'LineWidth', 2);
    xlabel('Time (ms)'); ylabel('Gating Variable');
    title('Gating Variables');
    legend('n (K^+)', 'm (Na^+)', 'h (Na^+)', 'Location', 'best');
    grid on;
end

%% Helper: compute currents from solution (returns INa, IK, IL)
function [INa, IK, IL] = currents_from_solution(y, params_local)
    Vvec = y(:,1);
    nvec = y(:,2);
    mvec = y(:,3);
    hvec = y(:,4);
    gNa = params_local.gNa_max .* (mvec.^3 .* hvec);
    gK  = params_local.gK_max  .* (nvec.^4);
    INa = gNa .* (Vvec - params_local.ENa);  % mS/cm2 * mV -> uA/cm^2
    IK  = gK  .* (Vvec - params_local.EK);
    IL  = params_local.gL .* (Vvec - params_local.EL);
end

%% Helper: initial gate steady-state at V (returns alpha/beta and steady-state)
function [rates, infs] = init_gates(V)
    % Use classic HH expressions (V in mV)
    % alpha_n
    if abs(V+55) < 1e-6
        alphan = 0.01/0.1; % limit
    else
        alphan = 0.01*(V+55)/(1 - exp(-(V+55)/10));
    end
    betan = 0.125*exp(-(V+65)/80);
    
    % alpha_m
    if abs(V+40) < 1e-6
        alpha_m = 0.1/0.1; % limit
    else
        alpha_m = 0.1*(V+40)/(1 - exp(-(V+40)/10));
    end
    beta_m  = 4*exp(-(V+65)/18);
    
    alpha_h = 0.07*exp(-(V+65)/20);
    beta_h  = 1/(1+exp(-(V+35)/10));
    
    rates = struct('an',alphan,'bn',betan,'am',alpha_m,'bm',beta_m,'ah',alpha_h,'bh',beta_h);
    infs.n = alphan/(alphan+betan);
    infs.m = alpha_m/(alpha_m+beta_m);
    infs.h = alpha_h/(alpha_h+beta_h);
end

function dy = hh_diff_eq(t, y, params)
% hh_diff_eq implements Hodgkin-Huxley ODEs for V, n, m, h
% y = [V; n; m; h]
% params: struct containing gNa_max, gK_max, gL, ENa, EK, EL, Cm, stimfun

V = y(1);
n = y(2);
m = y(3);
h = y(4);

% ---- HH alpha/beta (classic forms) ----
% alpha_n
if abs(V+55) < 1e-6
    an = 0.01/0.1;
else
    an = 0.01*(V+55)/(1 - exp(-(V+55)/10));
end
bn = 0.125*exp(-(V+65)/80);

% alpha_m
if abs(V+40) < 1e-6
    am = 0.1/0.1;
else
    am = 0.1*(V+40)/(1 - exp(-(V+40)/10));
end
bm = 4*exp(-(V+65)/18);

% alpha_h
ah = 0.07*exp(-(V+65)/20);
bh = 1/(1 + exp(-(V+35)/10));

% gating variable derivatives
dn = an*(1-n) - bn*n;
dm = am*(1-m) - bm*m;
dh = ah*(1-h) - bh*h;

% conductances
gNa = params.gNa_max * (m^3) * h;
gK  = params.gK_max  * (n^4);
gL  = params.gL;

% ionic currents (ionic currents positive outward)
INa = gNa * (V - params.ENa);  % uA/cm^2
IK  = gK  * (V - params.EK);
IL  = gL  * (V - params.EL);

% stimulus current
if isfield(params,'stimfun') && ~isempty(params.stimfun)
    Iext = params.stimfun(t); % uA/cm^2
else
    Iext = 0;
end

% Membrane equation: Cm dV/dt = - (INa + IK + IL) + Iext
% Note: Sign convention - positive stimulus current depolarizes membrane
dV = ( - (INa + IK + IL) + Iext ) / params.Cm;

dy = [dV; dn; dm; dh];
end