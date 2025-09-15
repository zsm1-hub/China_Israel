function [x1,y1,S3L1,S3T1]=test4_calc_SF3(u,v,N);
[X, Y] = meshgrid(-N/2:N/2-1, -N/2:N/2-1);

R = sqrt(X.^2 + Y.^2);
cost=X./R;sint=Y./R;

u2=u.^2;v2=v.^2;uv=u.*v;

% why not uh = fft2(u,2*N-1,2*N-1)? that will not account for aliasing
uh = fft2(u); % 
vh = fft2(v);
uuh= fft2(u2);
vvh= fft2(v2);
uvh= fft2(uv);

%%%%%% Cu_uu like u_hat.* conj((u.^2)_hat)
Cu_uu=fftshift(ifft2(uh.*conj(uuh)./N.^2));
Cuu_u=fftshift(ifft2(uuh.*conj(uh)./N.^2));
Cv_vv=fftshift(ifft2(vh.*conj(vvh)./N.^2));
Cvv_v=fftshift(ifft2(vvh.*conj(vh)./N.^2));

Cuv_u=fftshift(ifft2(uvh.*conj(uh)./N.^2));
Cuv_v=fftshift(ifft2(uvh.*conj(vh)./N.^2));
Cu_uv=fftshift(ifft2(uh.*conj(uvh)./N.^2));
Cv_uv=fftshift(ifft2(vh.*conj(uvh)./N.^2));

Cuu_v=fftshift(ifft2(uuh.*conj(vh)./N.^2));
Cvv_u=fftshift(ifft2(vvh.*conj(uh)./N.^2));
Cv_uu=fftshift(ifft2(vh.*conj(uuh)./N.^2));
Cu_vv=fftshift(ifft2(uh.*conj(vvh)./N.^2));
% 
Suuu=3.*Cu_uu-3.*Cuu_u;
Suvv=2.*Cv_uv-2.*Cuv_v-Cvv_u+Cu_vv;
Svvv=3.*Cv_vv-3.*Cvv_v;
Svuu=2.*Cu_uv-2.*Cuv_u-Cuu_v+Cv_uu;
%%%%%%%%%%%%%%%%%%% Jin-Han's code
% hSuuu_sum = hSuuu_sum - 3*huu.*conj(hu) + 3*conj(huu).*hu;
% 
%     hSuuw_sum = hSuuw_sum - 2*huw.*conj(hu) + 2*conj(huw).*hu - huu.*conj(hw) + conj(huu).*hw; 
% 
%     hSuww_sum = hSuww_sum - 2*huw.*conj(hw) + 2*conj(huw).*hw - hww.*conj(hu) + conj(hww).*hu; 
% 
%     hSwww_sum = hSwww_sum - 3*hww.*conj(hw) + 3*conj(hww).*hw;
%%%%%%%%%%%%%%%%%%%%% direction transfer
% \delta uL=\delta u \cos\theta + \delta v \sin\theta 
% \delta uT=-\delta u \sin\theta + \delta v \cos\theta 
S3L=cost.^3.*(Suuu)+.3.*cost.^2.*sint.*(Svuu)+3.*cost.*sint.^2.*(Suvv)+sint.^3.*(Svvv);



S3T=sint.^2.*cost.*(Suuu)+(sint.^3-2.*sint.*cost.^2).*Svuu+...
(cost.^3-2.*cost.*sint.^2).*Suvv+sint.*cost.^2.*(Svvv);

% only positive frequency part
S3L1=S3L(end/2+1:end,end/2+1:end);
S3T1=S3T(end/2+1:end,end/2+1:end);

x1=X(end/2+1:end,end/2+1:end);
y1=Y(end/2+1:end,end/2+1:end);


% [~, S3L1_iso] = calc_ispec2(x1(1,:), y1(:,1), S3L1, 2);
% [rbin, S3T1_iso] = calc_ispec2(x1(1,:), y1(:,1), S3T1, 2);

