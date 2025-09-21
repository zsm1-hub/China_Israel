function [rr,St_azimuth,Sl_azimuth,S2]=test2_calc_SF2(u,v,N);

% u=u1w;v=v1w;



[X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);
R = sqrt(X.^2 + Y.^2);
cost=X./R;sint=Y./R;

uh = (fft2(u)); % 二维傅里叶变换
vh = (fft2(v));

Cuu=fftshift(ifft2(uh.*conj(uh)./(N.^2)));
Cvv=fftshift(ifft2(vh.*conj(vh)./(N.^2)));

Cuv=fftshift(ifft2(uh.*conj(vh)./(N.^2)));
Cvu=fftshift(ifft2(vh.*conj(uh)./(N.^2)));

% Ruu=Cuu.*cost.^2+Cvv.*sint.^2+sint.*cost.*(Cuv+Cvu);
% Rvv=Cvv.*sint.^2+Cvv.*cost.^2+sint.*cost.*(Cuv+Cvu);
% Ruv=(Cuu-Cvv).*cost.*sint+Cvu.*sint.^2-Cuv.*cost.^2;
% Rvu=(Cuu-Cvv).*cost.*sint-Cvu.*cost.^2+Cuv.*sint.^2;

u2=mean(u(:).^2);v2=mean(v(:).^2);uv=mean(u(:)).*mean(v(:));

S_transverse=2.*(u2.*sint.^2+v2.*cost.^2-2.*uv.*sint.*cost)-...
    2.*(Cuu.*sint.^2-(Cuv+Cvu).*sint.*cost+Cvv.*cost.^2);

S_longitudinal=2.*(u2.*cost.^2+v2.*sint.^2+2.*uv.*sint.*cost)-...
    2.*(Cuu.*cost.^2+(Cuv+Cvu).*sint.*cost+Cvv.*sint.^2);


rx=-N/2:N/2-1;ry=-N/2:N/2-1;
rr=unique(sqrt(rx.^2+ry.^2));
dr=sqrt(abs(rx(2)-rx(1)).^2+abs(ry(2)-ry(1)).^2);

for ii=1:length(rr)
    pos=(R>=rr(ii)-dr/2) & (R<rr(ii)+dr/2);

    St_azimuth(ii)=mean(S_transverse(pos));
    Sl_azimuth(ii)=mean(S_longitudinal(pos));
end

S2=1/2.*(St_azimuth+Sl_azimuth);