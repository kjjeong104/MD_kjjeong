function density_gyroid_levelset(filename)

[v1,v2,grid,cell] = grid_data(filename);        % Read old rho_rgrid file

f_A = 0.35;                           		% Volume fraction fo A block that 
                                            % is forming the Network structure
phi_A = v1; phi_B = zeros(length(phi_A),1);
grid_points = grid(1)*grid(2)*grid(3);

file_id = fopen('rho_rgrid','w');
fprintf(file_id,'  format  1  0\n'); 			% Summary of the system in
fprintf(file_id,'dim\n'); 						% in.rho_kgrid file
fprintf(file_id,' \t       %d\n',3);
fprintf(file_id,'crystal_system\n');
fprintf(file_id,' \t      ''%s'' \n','cubic');
fprintf(file_id,'N_cell_param\n');
fprintf(file_id,' \t       %d \n',1);
fprintf(file_id,'cell_param\n');
fprintf(file_id,' \t       %6.4f \n',cell(1));
fprintf(file_id,'group_name\n');
fprintf(file_id,' \t      ''%s'' \n','I a -3 d');
fprintf(file_id,'N_monomer\n');
fprintf(file_id,' \t       %d \n',2);
fprintf(file_id,'ngrid\n');
fprintf(file_id,' \t      %i \t\t %i \t\t %i \n',grid(1),grid(2),grid(3));

N_bins = 1000;                    % Number for bins for the historgram
Edges = min(v1):(max(v1)-min(v1))/N_bins:max(v1); 
subplot(1,2,1); h = histogram(v1,Edges);    % Making and ploting histogram
title('Historgram of fraction of points having density within the bin');
xlabel('density');
ylabel('fraction of points');

f = h.Values; f_norm = f/sum(f);       % Computing fraction of points
f_cumm = cumsum(f_norm);               % Calculating cummulative fraction
subplot(1,2,2);bar(Edges(2:N_bins+1),f_cumm);  % Ploting cummulative histogram
title('Cumulative historgram');
xlabel('density');
ylabel('fraction of points');

[K,I] = min(abs(f_cumm-f_A));        % Finding the level-set value which gives
                                     % fraction of points most close to the 
                                     % desired volume fraction (f_A)
s = Edges(I);

points = find(v1<s);

vol_frac = length(points)/grid_points; % Fraction of points below level-set (s)

for i=1:grid_points
    
    if(phi_A(i,1) < s)
        phi_A(i,1) = 1;                   % Volume fraction < level-set = 1
    else
        phi_A(i,1) = 0;                   % Volume fraction > level-set = 0
    end
    
    phi_B(i,1) = 1 - phi_A(i,1);
    
    fprintf(file_id,'      %9.6f       %9.6f\n',phi_A(i,1),phi_B(i,1));
    
end

fprintf('\nLevel-set value: %6.4f \n', s) 
fprintf('Fraction of points below level-set: %6.4f \n',vol_frac);

end

function [v1,v2,grid,cell] = grid_data(filename)

C = textread(filename, '%s','delimiter', '\n');
dim = str2num(C{3});                    % Read dimension
cell = str2num(C{9});                   % Read cell parameter
grid = str2num(C{15});                  % Read grid

if(length(cell)==1)
   new_cell = ones(1,3)*cell;           % Cubic crystals 
elseif(length(cell)==2)
    new_cell(1:2) = cell(1);            % Tetragonal crystals
    new_cell(3)   = cell(2);
else
    new_cell = cell;                    % Orthorhombic crystals
end
 
if(length(grid)==1)
   grid(2) = grid(1);                   % 3D grid for 1D crystals
   grid(3) = grid(1);
elseif(length(grid)==2)
    grid(3) = grid(1);                  % 3D grid for 2D crystals
end

clear cell;cell = new_cell;
for i =16:length(C)    
    A(i-15,:) = str2num(C{i});
end
v1 = A(:,1); v2 = A(:,2);

end