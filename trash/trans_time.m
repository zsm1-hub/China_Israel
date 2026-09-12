function Var=trans_time(var,dt,Time)

Var=zeros(size(var));
for npart=1:size(var,2)
    TT=Time(:,npart);
    midvar=var(:,npart);
    
    % dt=3600;
    indice=round(TT(1,:)./dt);
    midvar2=zeros(size(midvar)).*nan;
    midvar2(1+indice:end)=midvar(1:size(midvar,1)-indice);
    Var(:,npart)=midvar2;
end
return