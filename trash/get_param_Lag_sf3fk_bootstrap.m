
function [SpecFlux_Lagm,Vt_Lagm,ebs_Lagm,kf_Lag,...
    lf_Lag,CI_ebs,CI_Vt,...
    CI_SpecFlux,Th_Lag]=get_param_Lag_sf3fk_bootstrap(Case,nparticles,lambda,timerange);

% fname=[Case,'_pars_P',num2str(nparticles),'T480bootstrap.mat'];
% fname=[Case,'_pars_P',num2str(nparticles),'T720bootstrap.mat'];
% fname=[Case,'_pars_P',num2str(nparticles),'T960bootstrap.mat'];
% fname=[Case,'_pars_P',num2str(nparticles),'T1440bootstrap.mat'];
fname=[Case,'_pars_P',num2str(nparticles),'T',num2str(timerange(end)),'bootstrap.mat'];
% load filtscale.mat
load(fname);
Th_Lag=nanmean(Th_all,2);
SF3_Lag=nanmean(SF3,2);
SF3_Lag_all=(SF3);

% kf1=1./filtscale.*2.*pi;
kf1=1;
for jj=1:size(SF3,2)
[SpecFlux_Lag(:,jj),Vt_Lag(:,jj),ebs_Lag(:,jj),kf_Lag,...
    lf_Lag]=Fk_fitting_SF3_Lcurve(SF3_Lag_all(:,jj), ...
    dist_axis,2.5e3,500e3,'log','RLS',lambda,kf1,0);
end

% clear CI_ebs CI_Vt CI_SpecFlux CI_SF3
% ebs_varp = ebs;
% ebs_varp(1:end-1,:) = ebs(1:end-1,:).*kf';

for i = 1:size(ebs_Lag,1)
    CI_ebs(:,i) = prctile(ebs_Lag(i,:), [95, 5]); 
end

for i = 1:size(Vt_Lag,1)
    CI_Vt(:,i) = prctile(Vt_Lag(i,:), [95,5]);
end

% for i = 1:size(SF3,1)
%     CI_SF3(:,i) = prctile(SF3(i,:), [95,5]);
% end

for i = 1:size(SpecFlux_Lag,1)    
    CI_SpecFlux(:,i) = prctile(SpecFlux_Lag(i,:), [95,5]);
end

SpecFlux_Lagm=nanmean(SpecFlux_Lag,2);
Vt_Lagm=nanmean(Vt_Lag,2);
ebs_Lagm=nanmean(ebs_Lag,2);
% [SpecFlux_Lag1time,~,~,~,...
%     ~]=Fk_fitting_SF3_Lcurve(nanmean(SF3_Lag_all,2), ...
%     dist_axis,2.5e3,500e3,'log','RLS',lambda,kf1,0);
% disp(['process  ',fname,' done']);
return