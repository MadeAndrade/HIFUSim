%% Authored by Joshua Soneson 2018, %%%%%%%%%%%%%%%%%%%
function[Grid,Layer,Q] = WAKZK_mat()
%% Implementation of wide-angle parabolic method for axisymmetric HITU
%% beams.
tic;


%% Transducer %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tx.f = 2e6;		% frequency (Hz)
Tx.a1 = 1;		% inner radius (cm)
Tx.a2 = 3;		% outer radius (cm)
Tx.d = 5;		% geometric focus (cm) [= Inf if planar array]
Tx.P = 200;		% total acoustic power (W)



%% Computational Grid %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Tx.d==Inf
    z_start = 0;
else
    z_start = Tx.d-sqrt(Tx.d^2-Tx.a2^2);	% axial location of equivalent source (cm)
end
Grid.Z = 7;		% max axial distance (cm)
Grid.KK = 128;		% max number of harmonics in simulation (use power of 2)
ppw_r = 30;		% grid resolution in r-direction; points per wavelength
ppw_z = 20;		% and z-direction


%% Spatial averaging %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hd = 0;		% diameter of hydrophone element (mm).  [hd = 0 no averaging]


%% Graphical output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z_output = [4.65 4.95 5.25];	% locations on z-axis where plots are produced
LL = length(z_output);	% code determines number of plot locations
ll = 1;			% initialize index


%% Layered media %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II = 3;			% number of layers

% material 1 parameters:
Layer(1).z = 0;		% material transition distance (cm)
Layer(1).c = 1482;	% small-signal sound speed (m/s)
Layer(1).rho = 1000;	% mass density (kg/m^3)
Layer(1).alpha = 0.217;	% attenuation at 1MHz (dB/m)
Layer(1).fraction = 0;	% fraction of attenuation due to absorption
Layer(1).eta = 2;	% exponent in attenuation power law
Layer(1).beta = 3.5;	% nonlinear parameter
Layer(1).Cp = 4180;	% heat capacity (J/kg/K)
Layer(1).kappa = 0.6;	% thermal conductivity (W/m/K)
Layer(1).w = 0;		% perfusion rate (kg/m^3/s)

% material 2 parameters:
Layer(2).z = 4;		% material transition distance (cm)
Layer(2).c = 1629;
Layer(2).rho = 1000;
Layer(2).alpha = 58;
Layer(2).fraction = 0.9;
Layer(2).eta = 1;
Layer(2).beta = 4.5;
Layer(2).Cp = 4180;
Layer(2).kappa = 0.6;
Layer(2).w = 20;

% material 3 parameters:
Layer(3).z = 6;		% material transition distance (cm)
Layer(3).c = 1482;
Layer(3).rho = 1000;
Layer(3).alpha = 0.217;
Layer(3).fraction = 0;
Layer(3).eta = 2;
Layer(3).beta = 3.5;
Layer(3).Cp = 4180;
Layer(3).kappa = 0.6;
Layer(3).w = 0;

layrho(1:II) = Layer.rho;
layc(1:II) = Layer.c;
layfrac(1:II) = Layer.fraction;
laybeta(1:II) = Layer.beta;
layeta(1:II) = Layer.eta;

% User may add more Layers if necessary: Layer(4), Layer(5), etc..
% Be sure to change the value of II above.

% calculate wavenumber in each layer:
for ii=1:II
    Layer(ii).k = 2*pi*Tx.f/(100*layc(ii));	% wavenumber (cm^-1)
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% More computational domain stuff %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Tx.d==Inf
    a2p = Tx.a2;
else
    a2p = Tx.a2*(Tx.d-z_start)/sqrt(Tx.d^2-Tx.a2^2);	% radius of equivalent source (a_2')
end

w = 1.05*a2p;		% width (radius) of physical domain (Tx radius + 5% pad)
lambda = 2*pi/Layer(1).k;	% wavelength (cm)
th = 2*lambda;			% PML thickness
Grid.R = w + th;		% max radius of computational domain (cm)
Grid.JJ = ceil(ppw_r*Grid.KK^0.35*Grid.R/lambda);	% Gridpoints in r-dir
Grid.NN = ceil(ppw_z*(Grid.Z-z_start)/lambda);      % Gridpoints in z-direction
% node vectors
Grid.r = linspace(0,Grid.R,Grid.JJ)';
Grid.z = linspace(z_start,Grid.Z,Grid.NN);
dr = Grid.r(2);
dz = Grid.z(2)-Grid.z(1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Equivalent source %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This source is a converging spherical wave at z=z_start bounded by a2p.
% Note: it's dimensionless:
if Tx.d==Inf
    A = ones(Grid.JJ,1).*Grid.r<a2p;
    a1p = Tx.a1;
else
    A = Tx.d*exp(-1i*Layer(1).k*sqrt(Grid.r.^2+(Tx.d-z_start)^2)).*(Grid.r<a2p) ...
        ./sqrt(Grid.r.^2+(Tx.d-z_start)^2);
    if(Tx.a1 ~= 0)		% if there's a hole in the Tx
        a1p = Tx.a1*(Tx.d-z_start)/sqrt(Tx.d^2-Tx.a1^2);
        A = A.*(Grid.r>a1p);
    end
end



% The user could specify a custom source here [any complex function A=A(r)].
%A = ;



% Apply a low-pass filter to the source:
A = SourceFilterH(Grid.r,A,Layer(1).k);

% Next scale the source by the appropriate pressure coefficient so that it has
% the proper total acoustic power:
integral = 2*pi*dr*trapz(abs(A).^2.*Grid.r);
integral = 1e-4*integral;		% convert units from cm^2 to m^2
p0 = sqrt(2*Layer(1).rho*Layer(1).c*Tx.P/integral);
A = p0*A;	% dimensionalize the boundary condition (units of Pa)



%%%%%%%%%%%%%%%%%%%%%%%%
%% Spatial averaging %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hr = 0.1*hd/2;		% convert to radius in cm


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate attenuation, dispersion %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

v = Tx.f*(1:Grid.KK)'/1e6;		% vector of frequencies
for ii=1:II
    Layer(ii).alpha = Layer(ii).alpha/8.686/100;% convert to Np/cm
    if layeta(ii) == 1                                 % linear media
        Layer(ii).alpha = Layer(ii).alpha*v.*(1+2*1i*log(v)/pi);
    else                                                  % everything else
        Layer(ii).alpha = Layer(ii).alpha*(v.^layeta(ii) ...
            - 1i*tan(pi*layeta(ii)/2)*(v.^layeta(ii) - v));
    end
end


%% some reporting
fprintf('\n\tWavelength = %2.2f mm\n',10*lambda)
fprintf('\tNode count\n')
fprintf('\t\tAxial %d\n',Grid.NN)
fprintf('\t\tRadial %d\n',Grid.JJ)
fprintf('\tGrid stepsize\n')
fprintf('\t\tdz = %2.2f mm\n',10*dz)
fprintf('\t\tdr = %2.2f mm\n',10*dr)



%% dependent variable (pressure) matrices:
p = zeros(Grid.JJ,Grid.KK);	% new pressure (n+1 th step)
q = zeros(Grid.JJ,Grid.KK);	% old pressure (n th step)
q(:,1) = A;			% apply boundary condition (source)

dr_max = Grid.r(Grid.JJ)-Grid.r(Grid.JJ-1);         % mesh spacing near PML
JJ_ = Grid.JJ - ceil(hr/dr_max);      % Index of radial limit where spatial averaging occurs
p_r = zeros(1,Grid.NN);	% peak rarefactional pressure
p_c = zeros(1,Grid.NN);	% peak compressional pressure
[p_r(1),p_c(1),p5(:,1)] = SynthAxScan(Grid.r,q,hr,JJ_);

I = zeros(Grid.JJ,Grid.NN);	% intensity
Q = zeros(Grid.JJ,Grid.NN);	% power density
% I_mock = zeros(Grid.JJ,1);
% Q_mock = I_mock;



%% PML stuff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
phi = 2*(Grid.r-Grid.R+th)/th.*(Grid.r>Grid.R-th);
phi = phi+(1-phi).*(Grid.r>Grid.R-th/2);
u = exp(-1i*pi*phi/4);
% Du2 = zeros(Grid.JJ,1);
Du2 = (-1i*pi/2/th).*u.*(Grid.r>Grid.R);
Du2 = Du2.*(Grid.r<Grid.R+th/2);
ur = sparse(diag(u./Grid.r));
ur(1,1) = 0;
u = sparse(diag(u));
Du = sparse(diag(Du2));

%% Build transverse Laplacian operator w/PML:
clear A;
e = ones(Grid.JJ,1);
D1 = spdiags([-e e],[-1 1],Grid.JJ,Grid.JJ)/2/dr;
D2 = spdiags([e -2*e e],-1:1,Grid.JJ,Grid.JJ)/dr/dr;
A = u*((ur+Du)*D1 + u*D2);
A(1,2) = 2*A(1,2);			% zero flux BC at r=0;


%% peripherals for nonlinear integrator:
Ppos = zeros(Grid.JJ,Grid.NN);
Ppos(:,1) = abs(q(:,1));
Pneg = zeros(Grid.JJ,Grid.NN);
Pneg(:,1) = -abs(q(:,1));
p5 = zeros(min(Grid.KK,5),Grid.NN);
p5(1,1) = abs(q(1,1));
w = zeros(Grid.NN,2*Grid.KK);                % waveform data vectors
Y = zeros(1,2*Grid.KK);
I_td = zeros(Grid.JJ,2);			% change in intensity
dt = 1/Tx.f/(2*Grid.KK-1);			% in seconds
t = 1e6*(0:dt:1/Tx.f);				% in us

for kk = 1:Grid.KK
    M(kk).P1 = zeros(size(speye(Grid.JJ)));
    M(kk).P2 = zeros(size(speye(Grid.JJ)));
end

%% more reporting:
fprintf('\t\tdt = %2.2f ns\n',1e9*dt)


%% find indices of first Gridpoint in each Layer:
Layer(1).index = 1;
for ii=2:II
    Layer(ii).index = ceil((Layer(ii).z-z_start)/dz)+1;
end
Layer(II+1).index = Grid.NN;



%% integration loop:
disp('Starting integration loop: ')
lj = Grid.JJ;
p11_I = speye(Grid.JJ);
p11_A = A;

for ii=1:II
    disp(['Initialising simulations for layer ' num2str(ii)])
    %build operators for Layer ii:
    lk = Layer(ii).k;

    % OBSOLETE FORMULATION, BELOW IS FASTER
    %     for kk=1:Grid.KK		 %this loop should not be parallelised
    %         [M(kk).P1,M(kk).P2,M(kk).P3] = ...
    %          BuildPade12operators(A,kk,dz,Layer(ii).k,lj);
    %
    % %         [M(kk).P1,M(kk).P2] = ...
    % %             BuildPade11operators(A,kk,dz,lk,lj);
    %     end

    %this Pade11 wins against FDA formulation
        for kk = 1:Grid.KK
            M(kk).P1 = p11_I + 0.25*(1-(lk*kk)*dz*(1i))*(p11_A/((kk*lk)^2));
            M(kk).P2 = p11_I + 0.25*(1+(lk*kk)*dz*(1i))*(p11_A/((kk*lk)^2));
        end

    %     this Pade12 wins against FDA formulation
%     for kk = 1:Grid.KK
%         s = 1i*kk*lk*dz;
%         A_12 = p11_A/((kk*lk)^2);
%         muplus = (3-2*(s^2)+1i*sqrt((((2*s+6)*s-6)*s-18)*s-9))/12/(1+s);
%         muminus = (3-2*(s^2)-1i*sqrt((((2*s+6)*s-6)*s-18)*s-9))/12/(1+s);
%         epsilon = ((s+3)*s+3)/6/(1+s);
%         M(kk).P1 = p11_I + muplus*A_12;
%         M(kk).P2 = p11_I + muminus*A_12;
%         M(kk).P3 = p11_I + epsilon*A_12;
%     end

    mu = (laybeta(ii)/2/layrho(ii)/layc(ii)^3)*(0.01*dz/dt);

    cutoff = Layer(ii).alpha(1)*layrho(ii)*layc(ii)^2 ...
        / laybeta(ii)/Layer(ii).k; % cutoff for nonlinearity




    %% PARALELLISE
    for nn=Layer(ii).index:Layer(ii+1).index-1

        %integrate nonlinear term:
        [p,w(nn+1,:),Ppos(:,nn+1),Pneg(:,nn+1),I_td(:,1)] = ...
            TDNL_mat(q,w(nn+1,:),Grid.KK,Grid.JJ,mu,cutoff,Ppos(:,nn),Pneg(:,nn),I_td(:,1));
        %attenuation/dispersion term and diffraction term:

        parfor kk=1:Grid.KK % 0.5 speedup
            p(:,kk) = p(:,kk).*exp(-Layer(ii).alpha(kk)*dz);
            %             p(:,kk) = M(kk).P1 \ (M(kk).P2 \ (M(kk).P3*p(:,kk)));	% for Pade 12
            p(:,kk) = M(kk).P1 \ (M(kk).P2*p(:,kk));			% for Pade 11
        end


        %%
        Norm = norm(p(:,1));
        if(isnan(Norm) || isinf(Norm))	% stop if something goes wrong
            fprintf('\tNaN or Inf detected! Simulation stopped at z = %d cm.\n',Grid.z(nn))
        end
        q = p;
        % this loop should not be paralellised

        for jj=1:Grid.JJ	% calculate intensity I and power density H
            I_next(jj,:) = sum(abs(p(jj,:)).^2)/2/layrho(ii)/layc(ii);
            Q_next(jj,:) = (sum(layfrac(ii)*real(Layer(ii).alpha).*abs(p(jj,:)').^2) ...
                + sum(I_td(jj,:))/dz/(2*Grid.KK-1))/layrho(ii)/layc(ii);
        end
        I(:,nn+1) = I_next;
        Q(:,nn+1) = Q_next;

        %collect/process data:

        if hd==0
            p5(:,nn+1) = p(1,1:min(Grid.KK,5))';
        else
            [p_r(nn+1),p_c(nn+1),p5(:,nn+1)] = SynthAxScan(Grid.r,p,hr,JJ_);
        end
        if ll<=LL
            if Grid.z(nn+1)>z_output(ll) && Grid.z(nn)<=z_output(ll)  % find special output locations
                if hd==0
                    SpecOut(ll).pr = Pneg(:,nn+1);
                    SpecOut(ll).pc = Ppos(:,nn+1);
                    SpecOut(ll).p5 = abs(p(:,1:min(5,Grid.KK)));
                    SpecOut(ll).w = w(nn+1,:);
                else
                    SpecOut(ll) = SynthRadScan(Grid.r,p,hr,JJ_);
                end
                SpecOut(ll).I = I(:,nn+1);
                ll = ll + 1;
            end
        end
    end

    %rescale pressure due to transmission at interface between layers ii and ii+1:
    if ii < II
        q = 2*Layer(ii+1).rho*Layer(ii+1).c*q ...
            / (layrho(ii)*layc(ii) + Layer(ii+1).rho*Layer(ii+1).c);
    end
    disp('Layer complete -------- ')
end
fprintf('\tTook %2.1f minutes.\n\n',toc/60)

%LinearHeating(Layer(2),Grid,I);	% assumes a focused beam and that
% the focus is in layer 2

Layer = Layer(1:II);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot routine %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure	% 3d plot showing axial and radial pressure amplitudes
r_ones = ones(1,length(SpecOut(1).p5));
z_zeros = zeros(1,length(Grid.z));
plot3(Grid.z,z_zeros,abs(p5)/1e6,'LineWidth',2)
hold on
for ll=1:LL
    plot3(z_output(ll)*r_ones,Grid.r(1:length(SpecOut(ll).p5)),SpecOut(ll).p5/1e6,'LineWidth',2)
end
xlabel('z (cm)')
ylabel('r (cm)')
zlabel('|p| (MPa)')
grid on

figure	% axial plots of amplitude of first 5 harmonics and intensity
subplot(211)
plot(Grid.z,abs(p5)/1e6,'LineWidth',2)
ylabel('|p| (MPa)')
grid on
subplot(212)
plot(Grid.z,I(1,:)/1e4,'LineWidth',2)
xlabel('z (cm)')
ylabel('I (W/cm^2)')
grid on

for ll=1:LL		% build plot label
    V(ll).str = strcat('z= ',num2str(z_output(ll)),' cm');
end

if(Grid.KK>1)
    figure	% temporal waveforms at specified axial locations
    hold on
    for ll=1:LL
        plot(t,SpecOut(ll).w/1e6,'LineWidth',2)
    end
    xlim([t(1) t(length(t))])
    legend(V.str)
    grid on
    xlabel('t (us)')
    ylabel('p (MPa)')
    figure	% axial plots of compressional and rarefactional pressure
    if(hd==0)
        p_c = Ppos(1,:);
        p_r = Pneg(1,:);
    end
    plot(Grid.z,p_c/1e6,Grid.z,p_r/1e6,'LineWidth',2)
    xlabel('z (cm)')
    ylabel('p (MPa)')
    grid on
end

for ll=1:LL
    figure	% radial plots of amplitude of first 5 harmonics and intensity
    subplot(211)	% at specified axial locations
    plot(Grid.r(1:length(SpecOut(ll).p5)),SpecOut(ll).p5/1e6,'LineWidth',2)
    ylabel('|p| (MPa)')
    title(V(ll).str)
    grid on
    subplot(212)
    plot(Grid.r,SpecOut(ll).I/1e4,'LineWidth',2)
    xlabel('r (cm)')
    ylabel('I (W/cm^2)')
    grid on
end

figure	% spatial distribution of field emphasizing low-amplitude variations
r = [-Grid.r(Grid.JJ:-1:2);Grid.r];
I = [I(Grid.JJ:-1:2,:);I];
imagesc(Grid.z,r,I.^0.2)
xlabel('z (cm)')
ylabel('r (cm)')
grid
