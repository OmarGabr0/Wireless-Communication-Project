%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part A %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;
clear;
clc;
pkg load communications;
%%%%%%%%%%%%%%%%%%%% Inverse Erlang function %%%%%%%%%%%%%%%%%%%%%%%%%%

function A = inverlangb(c, gos)
    fun = @(A) gos - (A^c/factorial(c)) / sum(A.^((0:c))./factorial(0:c));
    A = fzero(fun, [0, 1000]);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Hata function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function L = Hata(f,hm,hb,d)
    CH = 0.8 + (1.1 * log10(f) - 0.7) * hm - 1.5 * log10(f);
    L = 69.55 + 26.16 * log10(f) - 13.82 * log10(hb) ...
     - CH + (44.9 - 6.55 * log10(hb)) * log10(d);
end

%%%%%%%%%%%%%%%%%%%%%%% Valid N function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function N = valid_N(n)
    limit = 100; % limit of i & k
    i = 0:limit;
    valid = [];
    for k = 0:limit
        valid = [valid (i.^2 + k^2 + i.*k)];
    end
    valid = unique(sort(valid));
    while ~(ismember(n,valid))
        n = n + 1;
        if n > 30000
            print ('N is too large');
            break;
        end
    end
    N = n;
end

%%%%%%%%%%%%%%%%%%%%% Cluster Size function %%%%%%%%%%%%%%%%%%%%%%%%%%

function N = cluster_size(S, interference, SIR_dB, n)
    SIR = 10^(SIR_dB/10);
    N = ceil((1/3)*((interference*SIR)^(1/n) + 1)^2);
    N = valid_N(N);
end

%%%%%%%%%%%%%%%%%%%%%%%%% Acell function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Acell, Asector] = traffic_intensity(S, N, sectors, GOS)
    K = floor(S / (N * sectors)); % number of channels per cell
    Asector = inverlangb(K, GOS);
    Acell = Asector * sectors;
end

%%%%%%%%%%%%%%%%%%%%%%% No. of Cells function %%%%%%%%%%%%%%%%%%%%%%%%

function Cells = no_of_cells(Area, User_Den, Acell, Au)
    Users_Per_cell = Acell / Au; % number of users per cell
    Total_Users = User_Den * Area; % total number of users
    Cells = ceil(Total_Users / Users_Per_cell); % total number of cells
end

%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R = radius(Area, Cells)
    Area_Per_cell = Area / Cells;
    R = sqrt(Area_Per_cell / (3 * sqrt(3) / 2));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%


S = 340;  % Total number of channels : S= N*K
freq = 900; % Frequency in MHz
sensitivity = -95; % in db
Au = 0.025; % in erlangs
n = 4; % path loss exponent
h_BS = 20; % Base Station height
h_MS = 1.5; % Mobile Station height
Area = 100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s_sweep = [1, 3, 6];    % sectorization sweep
i_sweep = [6, 2, 1];    % interference sweep

%%% 1) SIR Sweep

sir_sweep = 1:0.1:30;
N_sweep = [];
figure(2);
for j=1:3
    for i=1:length(sir_sweep)
        N_sweep(i) = cluster_size(S, i_sweep(j), sir_sweep(i), n);
    end
    hold on;
    plot(sir_sweep, N_sweep);
    xlabel('SIRmin (dB)', 'FontSize', 14);
    ylabel('Cluster Size (cells)', 'FontSize', 14);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    title('SIRmin vs Cluster Size', 'FontSize', 16);
    grid on;
end

%%%  2) 𝑆𝐼𝑅𝑚𝑖𝑛 = 19𝑑𝐵 & user density = 1400 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2


% Plot the number of cells & traffic intensity per cell versus GOS (1% to 30%).

gos_sweep = 1:0.1:30;
gos_sweep = gos_sweep / 100;
cells_sweep = [];
Acell_sweep = [];
Asector_sweep = [];
User_Den = 1400;

for j=1:3
    N = cluster_size(S, i_sweep(j), 19, n);
    for i=1:length(gos_sweep)
        [Acell_sweep(i), Asector_sweep(i)] = traffic_intensity(S, N, s_sweep(j), gos_sweep(i));
        cells_sweep(i) = no_of_cells(Area, User_Den, Acell_sweep(i), Au);
    end

    figure(3);
    hold on;
    plot(gos_sweep, cells_sweep);
    xlabel('GOS (%)', 'FontSize', 14);
    ylabel('Number of Cells (cells)', 'FontSize', 14);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    title('GOS vs Number of Cells (SIRmin = 19 dB)', 'FontSize', 16);
    grid on;

    figure(4);
    hold on;
    plot(gos_sweep, Acell_sweep);
    xlabel('GOS (%)', 'FontSize', 14);
    ylabel('Traffic Intensity per Cell (Erlang)', 'FontSize', 14);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    title('GOS vs Traffic Intensity per Cell (SIRmin = 19 dB)', 'FontSize', 16);
    grid on;

end


% 3) At 𝑆𝐼𝑅𝑚𝑖𝑛 = 14𝑑𝐵 & user density= 1400 𝑢𝑠𝑒𝑟𝑠/𝑘𝑚2

gos_sweep = 1:0.1:30;
gos_sweep = gos_sweep / 100;
cells_sweep = [];
Acell_sweep = [];
Asector_sweep = [];
User_Den = 1400;

for j=1:3
    N = cluster_size(S, i_sweep(j), 14, n);
    for i=1:length(gos_sweep)
        [Acell_sweep(i), Asector_sweep(i)] = traffic_intensity(S, N, s_sweep(j), gos_sweep(i));
        cells_sweep(i) = no_of_cells(Area, User_Den, Acell_sweep(i), Au);
    end

    figure(5);
    hold on;
    plot(gos_sweep, cells_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('GOS (%)', 'FontSize', 14);
    ylabel('Number of Cells (cells)', 'FontSize', 14);
    title('GOS vs Number of Cells (SIRmin = 14 dB)', 'FontSize', 16);
    grid on;

    figure(6);
    hold on;
    plot(gos_sweep, Acell_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('GOS (%)', 'FontSize', 14);
    ylabel('Traffic Intensity per Cell (Erlang)', 'FontSize', 14);
    title('GOS vs Traffic Intensity per Cell (SIRmin = 14 dB)', 'FontSize', 16);
    grid on;
end

% 4) At 𝑆𝐼𝑅𝑚𝑖𝑛 = 14𝑑𝐵 & GOS= 2%

User_Den_sweep = 100:100:2000;
GOS = 2/100;
cells_sweep = [];
Acell_sweep = [];
Asector_sweep = [];
radius_sweep = [];


for j=1:3
    N = cluster_size(S, i_sweep(j), 14, n);
    for i=1:length(User_Den_sweep)
        [Acell_sweep(i), Asector_sweep(i)] = traffic_intensity(S, N, s_sweep(j), GOS);
        cells_sweep(i) = no_of_cells(Area, User_Den_sweep(i), Acell_sweep(i), Au);
        radius_sweep(i) = radius(Area, cells_sweep(i));
    end

    figure(7);
    hold on;
    plot(User_Den_sweep, cells_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('User Density (users/km^2)', 'FontSize', 14);
    ylabel('Number of Cells (cells)', 'FontSize', 14);
    title('User Density vs Number of Cells (SIRmin = 14 dB)', 'FontSize', 16);
    grid on;

    figure(8);
    hold on;
    plot(User_Den_sweep, radius_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('User Density (users/km^2)', 'FontSize', 14);
    ylabel('Cell Radius (km)', 'FontSize', 14);
    title('User Density vs Cell Radius (SIRmin = 14 dB)', 'FontSize', 16);
    grid on;
end

% 5) At 𝑆𝐼𝑅𝑚𝑖𝑛 = 19𝑑𝐵 & GOS= 2%

User_Den_sweep = 100:100:2000;
GOS = 2/100;
cells_sweep = [];
Acell_sweep = [];
Asector_sweep = [];
radius_sweep = [];


for j=1:3
    N = cluster_size(S, i_sweep(j), 19, n);
    for i=1:length(User_Den_sweep)
        [Acell_sweep(i), Asector_sweep(i)] = traffic_intensity(S, N, s_sweep(j), GOS);
        cells_sweep(i) = no_of_cells(Area, User_Den_sweep(i), Acell_sweep(i), Au);
        radius_sweep(i) = radius(Area, cells_sweep(i));
    end

    figure(9);
    hold on;
    plot(User_Den_sweep, cells_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('User Density (users/km^2)', 'FontSize', 14);
    ylabel('Number of Cells (cells)', 'FontSize', 14);
    title('User Density vs Number of Cells (SIRmin = 19 dB)', 'FontSize', 16);
    grid on;

    figure(10);
    hold on;
    plot(User_Den_sweep, radius_sweep);
    if j == 3
        legend('Omni-directional', '120° sectorization', '60° sectorization', 'FontSize', 12);
    end
    xlabel('User Density (users/km^2)', 'FontSize', 14);
    ylabel('Cell Radius (km)', 'FontSize', 14);
    title('User Density vs Cell Radius (SIRmin = 19 dB)', 'FontSize', 16);
    grid on;
end
