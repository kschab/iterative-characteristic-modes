close all
clear
clc

% options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kdex = 1      ka = 0.1\pi     (electrically small)
% kdex = 2      ka = \pi        (electrically medium)
% kdex = 3      ka = 10\pi      (electrically large)
%
% undersampling reduces sampling of incident and scattered planewaves in
% the azimuthal direction
%
% undersampling = 1         1 degree increments
% undersampling = 5         5 degree increments
% undersampling = 10        10 degree increments
%
% increasing undersampling generally...
%       1 -- reduces the accuracy of characteristic modes calculated by the
%       scattering dyadic (as compared to the XR formulation).  this is
%       particularly noticable for electrically large objects
%       2 -- decreases the speed-up realized by the iterative algorithm by
%       reducing the runtime of the full scattering dyadic calculation.
%       this, however, comes at the cost of accuracy (see note 1).
%
% fastflag precalculates the system matrix inverse for the fullwave
% simulation.  This is done purely for speed in running this demonstration
% code, though in practice such a matrix inverse might not be used in large
% problems where matrix-free methods are required.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin user settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kdex = 2;
undersampling = 1;
fastflag = 1;
plotting = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% end user settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% unpack data
load precalculated_mom_data.mat
Z = squeeze(Z_(:,:,kdex));
V = squeeze(V_(:,1:undersampling:end,kdex));
f = squeeze(f_(1:undersampling:end,:,kdex));
phi = phi(1:undersampling:end);
dphi = phi(2)-phi(1);
k0 = k0_(kdex);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% A -- IMPEDANCE-BASED CHARACTERISTIC MODE FORMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NMODES = [10,30,70];
Nmodes = NMODES(kdex);
[I,lam] = eigs(imag(Z),real(Z),Nmodes,'sm');
lam = diag(lam);
tlam = -1./(1+1j*lam);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B -- FULL SCATTERING DYADIC CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic();
pdex = 0;
S = zeros(length(phi));
for phi_ = phi
    pdex = pdex+1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % begin emulation of full-wave solver using precalculated data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate single plane wave excitation
    Vp = squeeze(V(:,pdex));

    % solve for currents using the impedance matrix
    if fastflag
        if pdex==1
            Zi = inv(Z);
        end
        I = Zi*Vp;
    else
        I = Z\Vp;
    end

    % generating far-field map and radiation patterns
    S(:,pdex) = f*I*dphi;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end emulation of full-wave solver using precalculated data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % plot
    if plotting
        if mod(pdex,5)==0
            figure(99)
            subplot(1,2,1)
            h = polarplot(phi.',abs(S(:,pdex)));
            set(h,'linewidth',2)
            hold on
            title({'full S construction';'scattering response';['run ',num2str(pdex),' / ',num2str(length(phi))]})
        end
    end
end

[~,t] = eigs(S,Nmodes,'lm');
t = diag(t)/sqrt(1j)*sqrt(k0)/(4*sqrt(pi/2));
s = 2*t+1;
TFULL = t;

if plotting
    figure(99)
    subplot(1,2,2)
    plot((cos(phi)-1)/2,(sin(phi))/2,'k:')
    hold on
    h1 = scatter(real(tlam),imag(tlam),'bo');
    axis equal
    h2 = scatter(real(t),imag(t),'r*');
    legend([h1,h2],{'XR formulation','full scattering dyadic'})
    xlabel('Re t_n')
    ylabel('Im t_n')
    title('eigenvalues t_n')
    disp('full S construction complete, press any key to begin iterative solution')
    waitforbuttonpress
else
    if fastflag == 0
        fullTime = toc();
        disp(['full calculation time: ',num2str(fullTime)])
        disp(['full calculation calls: ',num2str(pdex)])
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C -- ITERATIVE SCATTERING DYADIC APPROXIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic();
figure(101)
maxIter = size(S,1);   % max number of iterations
err = zeros(Nmodes,maxIter)*NaN;

% -- initialize iterative algorithm
m = 1;
a_ = [];
a_(:,m) = rand(size(S,1),1);
nA = 1;

% loop until max iterations reached or the algorithm is sufficiently
% converged.
sigModeError = 1;
figure(101)
while (m<maxIter) && sigModeError>1e-4

    % -- normalize most recent excitation
    a_(:,m) = a_(:,m)/sqrt(a_(:,m)'*a_(:,m)*dphi);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % begin emulation of full-wave solver using precalculated data
    % F(:,m) = L(a_(:,m))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % -- construct new excitation vector for simulation
    Vp = 0;
    pdex = 0;
    for phi_ = phi
        pdex = pdex+1;
        Vp = Vp+a_(pdex,m)*V(:,pdex)*dphi;
    end

    % -- invert system matrix to solve for currents
    if fastflag
        if pdex==1
            Zi = inv(Z);
        end
        I = Zi*Vp;
    else
        I = Z\Vp;
    end

    % -- calculate scattered field pattern due to most recent excitation
    F(:,m) = f*I;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % end emulation of full-wave solver using precalculated data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % -- create approximations of scattering dyadic and projection matrix
    Sm = 0;
    Pm = 0;
    for p = 1:m
        Sm = Sm+F(:,p)*a_(:,p)'*dphi;
        Pm = Pm+a_(:,p)*a_(:,p)';
    end

    % -- estimate eigenvalues
    [~,t_] = eigs(Sm,Nmodes,'lm');
    t_ = diag(t_)/sqrt(1j)*sqrt(k0)/(4*sqrt(pi/2));
    tm(:,m) = t_;

    % --  calculate indices of significant modes
    sigModes1EM2 = abs(t_)>0.01;
    numSigModes1EM2 = sum(sigModes1EM2);

    % -- calculate relative error
    if m>1
        err(:,m) = abs(tm(:,m)-tm(:,m-1))./abs(tm(:,m));
        sigModeError = max(abs(err(sigModes1EM2,m)));
    end

    % -- create excitation for next iteration
    a_(:,m+1) = F(:,m) - Pm*F(:,m)*dphi;
    nA = norm(a_(:,m+1));

    % -- update plots and output status
    if plotting
        if mod(m,1)==0
            cols = get(gca,'colororder');
            set(gcf,'position',[65 333 1274 464])
            subplot(1,3,2)
            cla
            h1 = scatter(real(tlam),imag(tlam),'bo');
            hold on
            h3 = plot(real(tm(:,m)),imag(tm(:,m)),'r*');
            plot((cos(phi)-1)/2,(sin(phi))/2,'k--')
            axis equal
            xlabel('Re t_n')
            ylabel('Im t_n')
            legend([h1,h3],{'XR formulation','latest iteration'})
            title('eigenvalues t_n')
            subplot(1,3,1)
            h = polarplot(phi.',abs(F(:,m)));
            set(h,'linewidth',2)
            hold on
            title({'iterative S estimation';'scattering response';['run ',num2str(m)]})
            drawnow()
            disp(['m = ',num2str(m),'    |a_m| = ',num2str(nA),'    sigModeError = ',num2str(sigModeError)])
            subplot(1,3,3)
            cla()
            imagesc(1:maxIter,1:Nmodes,log10(err))
            hold on
            plot(1:Nmodes,1:Nmodes,'r','linewidth',2)
            text(32,30,'m=n','color','r')
            yline(numSigModes1EM2,'r','t_n<0.01','linewidth',2)
            shading flat
            caxis([-10,0])
            colorbar
            ylabel('mode index')
            xlabel('iteration number')
            title('t_n estimate convergence')
            xlim([0,70])
        end
    end
    m = m+1;
    TITER = t_;
end

if (plotting==0) && (fastflag == 0)
    iterTime = toc();
    disp(['iterative calculation time: ',num2str(iterTime)])
    disp(['iterative calculation calls: ',num2str(m-1)])
    disp(['iterative calculation #modes: ',num2str(numSigModes1EM2)])

    figure()
    semilogy(abs(TFULL),'bo')
    hold on
    semilogy(abs(TITER),'r*')
    ylim([1e-3,1])
    yline(1e-2,'k--')
    xlabel('mode index')
    ylabel('modal significance')
end





