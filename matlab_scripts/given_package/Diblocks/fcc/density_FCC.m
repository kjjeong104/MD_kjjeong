function [rho_A, rho_B, wave] = density_FCC(filename)


C = textread(filename, '%s','delimiter', '\n'); % Read model paramters
N_monomer = str2num(C{5});                      % No. of Monomers
f = str2num(C{11});                             % Block lengths
dim = str2num(C{15});                           % No. of Monomers
crystal_system = C{17};                         % Crystal system
N_cell_param = str2num(C{19});                  % No of unit-cell parameters
cell_param = str2num(C{21});                    % Unit-cell dimensions 
grid = str2num(C{25});                          % Read grid
group_name = C{29};                             % Space group
output_filename = C{32}(2:end-1);               % Output filename(rho_kgrid)
N_particles = str2num(C{36});                   % Number of particles
pos = zeros(N_particles,3);
for i=38:38+N_particles-1
   pos(i-38+1,:) = str2num(C{i});               % Positions of particles
end


file_id = fopen(output_filename,'w');
fprintf(file_id,'  format  1  0\n'); 			% Summary of the system in
fprintf(file_id,'dim\n'); 						% output (rho_kgrid) file
fprintf(file_id,' \t       %d\n',dim);
fprintf(file_id,'crystal_system\n');
fprintf(file_id,' \t       %s \n',crystal_system);
fprintf(file_id,'N_cell_param\n');
fprintf(file_id,' \t       %i \n',N_cell_param);
fprintf(file_id,'cell_param\n');
fprintf(file_id,' \t       %6.4f \n',cell_param);
fprintf(file_id,'group_name\n');
fprintf(file_id,' \t       %s \n',group_name);
fprintf(file_id,'N_monomer\n');
fprintf(file_id,' \t       %d \n',N_monomer);
fprintf(file_id,'ngrid\n');
fprintf(file_id,' \t       %i \t\t %i \t\t %i \n',grid(1),grid(2),grid(3));


const = f(1)/N_particles; sigma_smear = 0.32;   % const = f/(No. of particles)              
Rsph = (3*f(1)/(16*pi))^(1/3);                  % Particle size
form_factor =@(x) 3*(sin(x) - x*cos(x))/x^3;    % x = qR
nodes = (grid(1)/2+1)*grid(2)*grid(3);          % Number of plane waves
rho_A=zeros(nodes,1); rho_B=zeros(nodes,1);
wave = zeros(nodes,3);

t = 0;
for i=0:grid(1)/2
    for j=0:grid(2)-1  
        for k=0:grid(3)-1
           
           G = [i j k];
           G = G_to_bz(G,grid,dim);             % Tranforming waves to first
                                                % brizilion zone (Aliasing)
            %First plane wave corresponds to the volume fraction, f.
          if(G(1)==0 && G(2)==0 && G(3)==0)
             t = t+1;
             wave(t,:)  = G;
             rho_A(t,1) = f(1);
             rho_B(t,1) = f(2);
            
            % Reflection condition for the FCC phase: All i,j,k should be even.
          elseif(mod(G(1),2)==0 && mod(G(2),2)==0 && mod(G(3),2)==0)
                
             t = t+1;
             wave(t,:) = G;
             qR = (sum(G.*G))^(0.5)*2*pi*Rsph;
             rho_A(t,1) = const*4*form_factor(qR)*exp(-(sigma_smear)^2*qR^2/2);
             rho_B(t,1) = -rho_A(t);
            
            % Reflection condition for the FCC phase: All i,j,k should be odd.
          elseif(mod(G(1),2)==1 && mod(G(2),2)==1 && mod(G(3),2)==1)    
                
             t = t+1;
             wave(t,:) = G;
             qR = (sum(G.*G))^(0.5)*2*pi*Rsph;
             rho_A(t,1) = const*4*form_factor(qR)*exp(-(sigma_smear)^2*qR^2/2);
             rho_B(t,1) = -rho_A(t);
            
            % Waves not satisfying the reflection conditions have rho = 0.
          else
              t = t+1;
              wave(t,:)  = G;
              rho_A(t,1) = 0;
              rho_B(t,1) = 0;
          end
            
          fprintf(file_id,'(%6.4E,%6.4E)  (%6.4E,%6.4E)\n', ...
                                           rho_A(t,1),0.0,rho_B(t,1),0.0);
            
       end
   end
end

end

% The function transform the wave to the first brizillion cell (Aliasing)
function Gout = G_to_bz(G,grid,dim)

Gout = G;
if(dim==3)
    if(G(2) > grid(2)/2)
        Gout(2) = G(2) - grid(2);
    end
    
    if(G(3) > grid(3)/2)
        Gout(3) = G(3) - grid(3);
    end
end

if(dim==2)
    if(G(2) > grid(2)/2)
        Gout(2) = G(2) - grid(2);
    end
end

end
