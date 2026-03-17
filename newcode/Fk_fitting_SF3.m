function [SpecFlux,Vt,ebs,kf,lf,kfo]=Fk_fitting_SF3(SF3, ...
    dist_axis,mindist,maxdist,kftype,inv_style,lambda,kf1);
% mindist=2,maxdist=500e3,kftype='log'
% SF3=(SF3lll_time+SF3ltt_time)';
nsamps = size(SF3,2);
disp(size(SF3));
disp('SF3 should be r × time ');
disp(size(dist_axis));
disp('dist_axis should be 1 × r');


r=dist_axis;
Nr=length(r);

% select part of the r axis that we think has reasonable
% data.
% modified by zsm
% ns=find(dist_axis>=2e3 ,1);
% ne=find(dist_axis<=500e3 ,1,'last');

ns=find(dist_axis>= mindist ,1);
ne=find(dist_axis<= maxdist ,1,'last');

R=r(ns:ne); 

NR=length(R);
lf=R;

if kftype=='log'
    % kf = logspace( log10(1/max(R)), log10(1/min(R)), 8).*2.*pi ;
    kf = logspace( log10(1/3.32), log10(1/0.04), 8).*2.*pi ;
else
    % error('what kf?');
    % kf=fliplr(1./lf).*2.*pi;
    kf=kf1;
end
kfo=kf;
dk = diff(kf);
kf = 0.5*(kf(1:end-1) + kf(2:end));

% kf=kf1;
% Code is written in a way that it is assumed that
% kf is decreasing. This matters when doing cumsum
kf = fliplr(kf); 
dk = fliplr(dk);

Nk=length(kf);
%%

% Define matrices
ebs = zeros(Nk+1,nsamps);
S = zeros(Nr, nsamps);
Vt = zeros(NR,nsamps);
SpecFlux = zeros(Nk, nsamps); 

norm_flag=1;

for n=1:nsamps
    
    % SF3 from the data
    %S(:,n) =s3lll(:,n)' +s3ltt(:,n)';
    S(:,n) = SF3(:,n);
    
    % SF3 over the selected range od points
    V=S(ns:ne, n)';
    
    % We will solve Ax=b problem
    % where x is the parameters (energy fluxes), b is SF3
    
    % build the matrix A
    A=zeros(NR,Nk+1);
    for j=1:Nk
        A(:,j)=-4/kf(j)*besselj(1,kf(j)*R)'*dk(j);
    end
    A(:,end)=2*R';
    
    % %% divide Ax=b by r on both sides to remove log dependence
    % normalize the magnitude
    if norm_flag == 1
        for j=1:NR
            A(j,:)=A(j,:)./abs(R(j));
        end
        V=V./abs(R);
    end
    
    % estimate the epsilons using least squares
    if strcmp(inv_style, 'NNLS')
        ebs(:,n) = lsqnonneg(A,V');
    elseif strcmp(inv_style, 'LS')
        ebs(:,n) = A\(V');

    elseif strcmp(inv_style, 'RLS')
        % 确保V是列向量
        V = V(:);  % 将行向量转为列向量
        
        % 获取A的列数（未知数个数）
        n_cols = size(A, 2);  % 应该是Nk+1 (27)
        
        % 构建正则化系统
        A_aug = [A; sqrt(lambda) * eye(n_cols)];
        V_aug = [V; zeros(n_cols, 1)];
        
        % 求解并确保结果是列向量
        ebs(:, n) = A_aug \ V_aug;
    end
    
    % What does the reconstructed V look like
    Vt(:,n) =2*ebs(end,n)*R;
    for j=1:Nk
        Vt(:,n)=Vt(:,n)-4*ebs(j,n)./kf(j).*besselj(1,kf(j).*R)'*dk(j);
    end
    
    % Convert V to spectral fluxes   
    SpecFlux(1,n) = -ebs(end,n); 
    
    for j = 2:Nk
        SpecFlux(j,n) = SpecFlux(j-1,n) + ebs(end-j+1, n)*dk(end-j+2);
    end
    
end
SpecFlux = flipud(SpecFlux); 
% subplot(222)
% semilogx(dist_axis./1e3,SF3'); hold on
% semilogx(lf./1e3,Vt)
% title('SF3')
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% grid on
% 
% 
% subplot(221)
% SpecFlux = flipud(SpecFlux); 
% semilogx(1./kf./1e3,SpecFlux,'Marker','x')
% title('F(k)')
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% grid on
% 
% 
% 
% subplot(223)
% semilogx(1./kf./1e3,ebs(1:end-1),'Marker','x')
% title('ebs')
% xlabel('km')
% grid on
% 
% 
% eval(['load ',Case,'_pars_P',num2str(nparticles),'T',num2str(days),'daysCG_Lag_tukey.mat'])
% fname='s2sflux_spec_smooth.0002.nc';
% ncdisp(fname)
% Thm_Eulerian=ncread(fname,'Thm');
% filtscale=ncread(fname,'filtscale');
% 
% 
% subplot(224)
% semilogx(1./kf./1e3,SpecFlux,'Marker','x'); hold on
% semilogx(filtscale./1e3,Th')
% semilogx(filtscale./1e3,Thm_Eulerian,'LineWidth',1.5,'Color',[.7,.7,.7]);
% 
% title('F(k) vs CG flux')
% xlabel('km')
% ylabel('m^{2}/s^{3}')
% grid on