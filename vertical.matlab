%
% 2D Elastic wave propagtaion problem with cracks
% Staggered grid
%      
%
%     (j+1/2) --- Txz-----uz-----Txz----uz-----Txz
%                  |      |      |      |      |
%                  |      |      |      |      |
%                  |      |      |      |      |
%         (j) --- ux---Txx,Tzz---ux---Txx,Tzz---ux
%                  |      |      |      |      |
%                  |      |      |      |      |
%                  |      |      |      |      |
%     (j-1/2) --- Txz-----uz-----Txz-----uz-----Txz
%                         |      |      |
%                        (i)  (i+1/2) (i+1)
% In the code:
% Txx,Tzz(i,j) = Txx,Tzz(i,j);
% Txz(i,j) = Txz(i+1/2,j+1/2);
% ux(i,j) = ux(i+1/2,j);
% uz(i,j) = uz(i,j+1/2);

% paper  Simulations of P-SV wave scattering due to cracks by the 2-D finite difference method

clear all;

%Box dimensions
L       =   6000;      %   Width of box    [m]     
H       =   6000;      %   Height of box   [m] 




% Numerical parameters
nx      =   300;         %   # gridpoints in x-direction
dh      =   L/(nx-1);    %   Spacing of grid
%nz      =   floor(H/dh/2)*2+1;         %   # gridpoints in z-direction
nz      =   300; 

%%% choose nt to avoid the reflections from the boundaries
nt      =   2000;          %   Number of timesteps to compute

vp = 5000;
vs = 2000;
r = 2000;

% The phasical properties grid
Vp      =   vp*ones(nx,nz);      %   Compressional wave velocity [m/s]
Vs      =   vs*ones(nx,nz);      %   Shear wave velocity [m/s]


Vs(2:4,:)    = zeros(3,nz); 
%Vp(:,nx / 2 - 5 : nx/2 +5) = Vp(:,nx / 2 - 5 : nx/2 +5) * 0.8

% Add heterogeneity here
% ????
rho     =   r*ones(nx,nz);      %   Density [kg/m^3]


miu     =   Vs.^2*r;  %   MU
lambda  =   (Vp.^2*r-2*miu);   %   Lambda

% %think layer of water
% Vs(2:30,:) = zeros(29,nz);
% miu(2:30,:) = zeros(29,nz);


% create grid
%[z2d,x2d] = meshgrid(0:dh:H, -L/2:dh:L/2);  
[z2d,x2d] = meshgrid(0:dh:H, 0:dh:L);  


% Compute stable timestep -- need prove 
dt   = 0.8*dh/(vp*sqrt(2))


% Source time function
half_dur = 0.2;                              % Source half duration [s]


% Wavelength
wl = min(min(Vs))*2*half_dur
% Setup initial velocity and stress profile
ux = zeros(nx,nz);
uz = zeros(nx,nz);
uz_above = zeros(nx,nz);
uz_below = zeros(nx,nz);


Txx = zeros(nx,nz);
Tzz = zeros(nx,nz);
Txz = zeros(nx,nz);

% Source location
src_nx  =   floor(nx/2)+1;      % source node at x-direction
%src_nz  =   floor(nz/2)+1;      % source node at z-direction

src_nz  = 3;


%crack dimension 
a = 150;
crack_tip = 30; 
crack_tip_end = 80;

crack_plane = linspace(crack_tip, crack_tip_end, crack_tip_end-crack_tip + 1);

% 2 cracks
crack = linspace(crack_tip, crack_tip_end, crack_tip_end-crack_tip + 1);


% Station location
% Add more stations
sta_x = [500 750 1000 1250;];
sta_z = [3000 3000 3000 3000];
sta_num = length(sta_x);
sta_nx = floor((sta_x+L/2)/dh)+1;
sta_nz = floor(sta_z/dh)+1;

 %Txz = 0 on a sequence of horizontally arrayed grid points no matter
offset = src_nx;


% Bouyancy and other parameter
b = 1./rho.*ones(nx,nz);
a = dt/dh;

BU = b*a;
LAM = lambda*a;
MU = miu*a;
GAMMA = LAM + 2*MU;

time = 0;
for n=1:nt
    %Add the source term
    if(time<=2*half_dur)
        % P type source
        % typically use this one if you want S type source, you can compare
        % with the P type source
        %UZ = uz(src_nx,src_nz);
        Txz(src_nx,src_nz) =  Txz(src_nx,src_nz) + exp(-100*(time-half_dur).^2);      
    end
    % Update stress from velocity
    % In the code:
    % ux(i,j) = ux(i+1/2,j);
    % uz(i,j) = uz(i,j+1/2);
    for i=2:nx-1
        for j =2:nz-1
            % below crack
             if i == offset && ismember(j,crack_plane) && mode(i,2) == 0
               
                Txx(i,j) = Txx(i,j) + GAMMA(i,j)*(uz_below(i-1,j)-uz_below(i,j)) + LAM(i,j)*(uz(i,j)-uz(i,j-1));
                Tzz(i,j) = Tzz(i,j) + GAMMA(i,j)*(uz(i,j)- uz(i,j-1)) + LAM(i,j)*(uz_below(i-1,j)-uz_below(i,j));
                
            % above crack
            elseif i == offset && ismember(j,crack_plane) && mode(i,2) == 1
                Txx(i,j) = Txx(i,j) + GAMMA(i,j)*(uz_above(i-1,j)-uz_above(i,j)) + LAM(i,j)*(uz(i,j)-uz(i,j-1));
                Tzz(i,j) = Tzz(i,j) + GAMMA(i,j)*(uz(i,j)- uz(i,j-1)) + LAM(i,j)*(uz_abovex(i-1,j)-uz_above(i,j));
                
            else
                Txx(i,j) = Txx(i,j) + GAMMA(i,j)*(ux(i,j)-ux(i-1,j)) + LAM(i,j)*(uz(i,j)-uz(i,j-1));
                Tzz(i,j) = Tzz(i,j) + GAMMA(i,j)*(uz(i,j)- uz(i,j-1)) + LAM(i,j)*(ux(i,j)-ux(i-1,j));
             end
            Txz(i,j) = Txz(i,j) + MU(i,j) * (uz(i+1,j)- uz(i,j)) + MU(i,j)*(ux(i,j+1)-ux(i,j)) ;
        end
    end
    
    % all Txz is zero in the crack
    %Txz(crack_tip_left:crack_tip_right, offset) = zeros(length(crack_plane), 1) ;
    Txz(offset, crack_tip: crack_tip_end) = zeros(1,length(crack_plane));



    % Set boundary conditions  Free Boundary:
    Txx([1 nx],:)   =  0;
    Txx(:,[1 nz])   =  0;
    Tzz([1 nx],:)   =  0;
    Tzz(:,[1 nz])   =  0;
    Txz([1 nx],:)   =  0;
    Txz(:,[1 nz])   =  0;    
    

    % Update velocity from stress
    % Txx,Tzz(i,j) = Txx,Tzz(i,j);
    % Txz(i,j) = Txz(i+1/2,j+1/2);
    
    for i=2:nx-1
        for j=2:nz-1
            uz(i,j) = uz(i,j) + BU(i,j) * (Txz(i,j) - Txz(i - 1,j)) + BU(i,j) * (Tzz(i,j + 1) -  Tzz(i,j));
                
             %below crack
            if i == offset && ( ismember(j,crack))
               uz_below(i,j)= uz_below(i,j) + (1/3) * BU(i,j) *(9 * Txx(i-1, j)- Txx(i, j));
               uz_above(i,j) = uz_above(i,j) - (1/3) * BU(i,j) *(9 * Txx(i+1, j)-Txx(i+2, j)); 
%             elseif j == offset  && ismember(i,crack_plane)
%                uz(i,j) = 0;
            else
                ux(i,j) = ux(i,j) + BU(i,j) * (Txx(i+1,j)-Txx(i,j)) + BU(i,j) * (Txz(i,j) - Txz(i,j-1)) ;
            end        
        end
    end
      
    
    
    % Set boundary conditions  fixed boundary
    ux([1 nx],:)   =  0;
    ux(:,[1 nz])   =  0;
    uz([1 nx],:)   =  0;
    uz(:,[1 nz])   =  0;    
    % Time info
    time        =   time+dt;

    
    % For pcolor plot
     % For pcolor plot
    PP = ux;
    
    % Save the synthetics at station
    time_save(n) = time;
    
    % save the ux at stations
    %U_save(n) = ?;
    % save the uz at stations
    %V_save = ?;
    
    % Plot solution every 50 timesteps
    if (mod(n,5)==0)
        figure(1), clf
        pcolor(x2d,z2d,PP); shading interp,  colorbar
        colormap winter;

        hold on
        %crack position
        %plot(x2d(crack_tip_left:crack_tip_right - 4,offset),z2d(crack_tip_left +4 :crack_tip_right,offset),'w-','LineWidth',2);
        plot(x2d(offset,crack_tip:crack_tip_end),z2d(offset,crack_tip:crack_tip_end),'w-','LineWidth',2);
                %source position
        plot(x2d(src_nx,src_nz),z2d(src_nx,src_nz),'rp');
        xlabel('x [m]')
        ylabel('z [m]')
        zlabel('Pressure [Pa]')
        title(['Wave propagation after T = ',num2str(time),' sec'])
        axis equal, axis tight
        drawnow
    end
    %pause
end