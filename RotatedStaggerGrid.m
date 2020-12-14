
%      (j+1/2) -   uzux1-----------------------uzux2
%                  |                           |
%                  |                           |       
%                  |                           |      
%          (j) --  | ------- Txx,zz,xz         |      
%                  |             |             |    
%                  |             |             |   
%                  |             |             |
%     (j-1/2) -     uzux3 ---------------------- uzux4
%                  |             |             |
%                 (i-1/2)        i          (i+1/2)
% In the code:

% Txx,Tzz,Txz(i,j) = Txx,Tzz,Txz(i,j);
% ux,uz(i,j) = ux(i+1/2,j + 1/2);
% --------nEW ------        
%  Dx(u) = [u( x + dx/2 , z + dz/2 , t) -  u( x - dx/2 , z - dz/2 , t)]
        % =  2 -3
%  Dz(u) = [u( x - dx/2 , z + dz/2 , t) -  u( x + dx/2 , z - dz/2 , t)] 
        % =  1 - 4

clear all;

%Box dimensions
L       =   6000;      %   Width of box    [m]     
H       =   6000;      %   Height of box   [m] 

% Numerical parameters
nx      =   1000;         %   # gridpoints in x-direction
dh      =   L/(nx-1);    %   Spacing of grid
%nz      =   floor(H/dh/2)*2+1;         %   # gridpoints in z-direction
nz      =   1000; 

%%% choose nt to avoid the reflections from the boundaries
nt      =   6000;          %   Number of timesteps to compute

vp =2000;
vs = 2000;
r = 2120;

% The phasical properties grid
Vp      =   vp*ones(nx,nz);      %   Compressional wave velocity [m/s]
Vs      =   vs*ones(nx,nz);      %   Shear wave velocity [m/s]

rho     =   r*ones(nx,nz);      %   Density [kg/m^3] 


miu     =   Vs.^2*r;  %   MU
lambda  =   (Vp.^2*r-2*miu);   %   Lambda

% %think layer of water
% Vs(2:30,:) = zeros(29,nz);
% miu(2:30,:) = zeros(29,nz);


% create grid
[z2d,x2d] = meshgrid(0:dh:H, -L/2:dh:L/2);  



% Compute stable timestep -- need prove 
dt   = 0.9 * dh/(vp*sqrt(2));

dr  = sqrt(dh^2 + dh^2);



% Source time function
half_dur = 0.2;                              % Source half duration [s]


% Wavelength
wl = min(min(Vs))*2*half_dur
% Setup initial velocity and stress profile
ux = zeros(nx,nz);
uz = zeros(nx,nz);
Txx = zeros(nx,nz);
Tzz = zeros(nx,nz);
Txz = zeros(nx,nz);

% Source location
src_nx  =   floor(nx/2)+1;      % source node at x-direction
src_nz  =   floor(nz/2)+1;      % source node at z-direction


% Bouyancy and other parameter
b = 1./rho.*ones(nx,nz);
a = dt/dh;

BU = b*a;
LAM = lambda*a;
MU = miu*a;
GAMMA = LAM + 2*MU;
R = dr/(2*dh);
K = 0.7* dh;

time = 0;
for n=1:nt
    %Add the source term
    if(time<=2*half_dur)
        % P type source
        % typically use this one if you want S type source, you can compare
        % with the P type source
        Txz(src_nx,src_nz) =  Txz(src_nx,src_nz)+source_time(time,half_dur);      
    end
    % Update stress from velocity
    % In the code:
    % ux(i,j) = ux(i+1/2,j);
    % uz(i,j) = uz(i,j+1/2);
    for i=2:nx-1
        for j=2:nz-1
            
            OpUx_x = (ux(i+1,j+1) - ux(i,j))/dr;
            OpUx_z = (ux(i,j + 1 ) - ux(i+1,j))/dr;
            
            %horizontal rotation
            dUx_dx = R * (OpUx_x  - OpUx_z) ;
            %vertical
            dUx_dz = R * (OpUx_x  + OpUx_z);
            
            OpUz_x = (uz(i + 1,j + 1) - uz(i,j))/dr;
            OpUz_z = (uz(i,j + 1) - uz(i + 1,j))/dr;


            %horizontal
            dUz_dx = R * (OpUz_x - OpUz_z) ;
            %vertical
            dUz_dz = R * (OpUz_x + OpUz_z);
         
        
            
            Txx(i,j) = Txx(i,j) + GAMMA(i,j)* dUx_dx + LAM(i,j)* dUz_dz;
            Tzz(i,j) = Tzz(i,j) + GAMMA(i,j)* dUz_dz + LAM(i,j)* dUx_dx;
            Txz(i,j) = Txz(i,j) + MU(i,j) * dUx_dz + MU(i,j)* dUz_dx ;
        end
    end

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
            %roatated grid
            %OpUx = (ux(i,j) - ux(i-1,j-1))/dr;
            %OpUz = (uz(i-1,j) - uz(i,j-1))/dr;
        
         
            
            OpTxx_X = (Txx(i,j) - Txx(i-1,j-1))/dr;
            OpTxx_Z = (Txx(i-1,j) - Txx(i,j-1))/dr;

            OpTxz_X = (Txz(i,j) - Txz(i-1,j-1))/dr;
            OpTxz_Z = (Txz(i-1,j) - Txz(i,j-1))/dr;
            
            OpTzz_X = (Tzz(i,j) - Tzz(i-1,j-1))/dr;
            OpTzz_Z = (Tzz(i-1,j) - Tzz(i,j-1))/dr;
            
           
            
            dTxx_dx = R * (OpTxx_X  - OpTxx_Z);
            dTzz_dz = R * (OpTzz_Z + OpTzz_X);
            
            dTxz_dx = R * (OpTxz_X  - OpTxz_Z);
            dTxz_dz = R * (OpTxz_Z + OpTxz_X);
         
            ux(i,j) = ux(i,j) + BU(i,j) *  dTxx_dx + BU(i,j) *  dTxz_dz + K *(ux(i-1,j) + ux(i+1,j) + ux(i,j-1) + ux(i,j+1) -4*ux(i,j)) / dh^2;
            uz(i,j) = uz(i,j) + BU(i,j) * dTxz_dx + BU(i,j) * dTzz_dz + K *(uz(i-1,j) + uz(i+1,j) + uz(i,j-1) + uz(i,j+1) -4*uz(i,j)) / dh^2;
            
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
    PP = ux;
    
    % Save the synthetics at station
    time_save(n) = time;
    
    
    % Plot solution every 50 timesteps
    if (mod(n,5)==0)
        figure(1), clf
        pcolor(x2d,z2d,PP); shading interp,  colorbar
        colormap winter
        hold on
        plot(x2d(src_nx,src_nz),z2d(src_nx,src_nz),'rp','MarkerSize',12,'MarkerFaceColor','red');
        xlabel('x [m]')
        ylabel('z [m]')
        zlabel('Pressure [Pa]')
        title(['Wave propagation after T = ',num2str(time),' sec'])
        axis equal, axis tight
        drawnow
    end
    %pause
end
figure;
% Plot the station records
%subplot(2,1,1), plot(time_save,U_save), xlabel('Time [sec]'), ylabel('Raidal Amp');
%subplot(2,1,1), plot(time_save,V_save), xlabel('Time [sec]'), ylabel('Vertical Amp');