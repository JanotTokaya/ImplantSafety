function [sigma_est1, sigma_est1l ,sigma_est1u ] = sigma_from_Ptrans_homogeneous_region(PT , X_PT, Y_PT, Z_PT)
%% J.Tokaya 23-11-2017
% This function is used to etimate the conductivity using phase only EPT 
% [Wen H. doi: 10.1117/12.480000]. A second degree poly fit is performed,
% i.e. Ptrans(X,Y,Z) = C_0+C_1*(X-X0).^2+C_2*(Y-Y0).^2+C_3*(Z-Z0).^2. From 
% C_1, C_2, C_3 the conductivity is estimated.
% INPUT - PT, a KxLxM matrix containing the measured transceive phase.
%       - X_PT,a KxLxM matrix containing x-coordinates of PT.
%       - Y_PT,a KxLxM matrix containing y-coordinates of PT.
%       - Z_PT,a KxLxM matrix containing z-coordinates of PT.
% OUTPUT- sigma_est1, the estimated conductivity.
%       - sigma_est1l, the lower bound of the estimated conductivity.
%       - sigma_est1u, the upper bound of the estimated conductivity

%% define some electromagnetic constants (NOTE: this estimate is for 1.5T)
mu0=1.2566370614*10^-6; 
eps0=8.85418782*10^-12;
f=63.87*10^6; 
omega = 2*pi*f; 
%% Determine conductivity:
% the numerical derivative below is noise sensitive and often doesn't work.
% dX = abs(mean(unique(diff(X_PT,1,1)))); 
% dY = abs(mean(unique(diff(Y_PT,1,2)))); 
% dZ = abs(mean(unique(diff(Z_PT,1,3)))); 
% L = del2(PT,dX,dY,dZ); 
% sigma_est2 = L/(2*mu0*omega); 

s=size(PT);
XY_data = zeros(numel(PT),3);
XY_data(:,1) = X_PT(:);
XY_data(:,2) = Y_PT(:);
XY_data(:,3) = Z_PT(:);
                   
%first find the origin of the poly, i.e. (X0,Y0,Z0). This is done 1D for 
%the middle profiles of the transceive phase.
poly2 = 'a*(x-b)^2+c';
PT_midx = double(smooth(squeeze(PT(:,round(end/2),round(end/2)))));

PT_midy = double(smooth(squeeze(PT(round(end/2),:,round(end/2)))));

PT_midz = double(smooth(squeeze(PT(round(end/2),round(end/2),:))));
const = zeros(3,1);

figure;
Xax= squeeze(X_PT(:,round(end/2),round(end/2))); 
startPoints = [(max(PT_midx)-min(PT_midx))/(0.5*(max(Xax)-min(Xax))), Xax(round(s(1)/2)) mean(PT_midx)];
f1 = fit(Xax,PT_midx,poly2,'Start', startPoints);
subplot(231),
plot(Xax,PT_midx)
hold on
plot(f1)
x_off=f1.b;
I=f1.a; 
const(1) = f1.c;
xlabel('x_{pos} (m)')
ylabel('PT (rad)')

Yax= squeeze(Y_PT(round(end/2),:,round(end/2))).'; 
startPoints = [(max(PT_midy)-min(PT_midy))/(0.5*(max(Yax)-min(Yax))), Yax(round(s(2)/2)) mean(PT_midy)];
f1 = fit(Yax,PT_midy,poly2,'Start', startPoints);
subplot(232),
plot(squeeze(Y_PT(round(end/2),:,round(end/2))).',PT_midy)
hold on
plot(f1)
y_off=f1.b;
J=f1.a; 
const(2) = f1.c;
xlabel('y_{pos} (m)')
ylabel('PT (rad)')

Zax= squeeze(Z_PT(round(end/2),round(end/2),:)); 
startPoints = [(max(PT_midz)-min(PT_midz))/(0.5*(max(Zax)-min(Zax))), Zax(round(s(3)/2)) mean(PT_midz)];
f1 = fit(Zax,PT_midz,poly2,'Start', startPoints);
subplot(233),
plot(Zax,PT_midz)
hold on
plot(f1)
z_off=f1.b;
K=f1.a; 
const(3) = f1.c;
xlabel('z_{pos} (m)')
ylabel('PT (rad)')

%initial guess
C = [mean(const) I J K]; 

%the polynomial that will fit the transceive phase distribution.
phi_of_xyz_no_linear_terms = @(C,XY_f)  C(1)...
                       +C(2)*(XY_f(:,1)-x_off).^2 ...
                       +C(3)*(XY_f(:,2)-y_off).^2 ...
                       +C(4)*(XY_f(:,3)-z_off).^2;
              
[C_est,~,resid,~,~,~,J] = lsqcurvefit(phi_of_xyz_no_linear_terms, C, XY_data, 0.5*double(PT(:)),[],[]); 

ci = nlparci(C_est,resid,'jacobian',J); 

sigma_est1 = 2*(C_est(2)+C_est(3)+C_est(4))/(mu0*omega);
sigma_est1l = 2*(ci(2,1)+ci(3,1)+ci(4,1))/(mu0*omega);
sigma_est1u = 2*(ci(2,2)+ci(3,2)+ci(4,2))/(mu0*omega);

%sigma_unc = 2/(mu0*omega)*sqrt(deltas(2)^2+deltas(3)^2+deltas(4)^2);

fit_results_PT = reshape(phi_of_xyz_no_linear_terms(C_est,XY_data), s);
%plot the resulting fitted transceive phase distribution.
subplot(4,3,7), imagesc(fit_results_PT(:,:,round(end/2)))
c=caxis;
subplot(4,3,10), imagesc(0.5*PT(:,:,round(end/2)),c)
subplot(4,3,8), imagesc(squeeze(fit_results_PT(round(end/2),:,:)))
title('fitted phase')
c=caxis;
subplot(4,3,11), imagesc(squeeze(0.5*PT(round(end/2),:,:)),c)
title('transceive phase')
subplot(4,3,9), imagesc(squeeze(fit_results_PT(:,round(end/2),:)))
c=caxis;
subplot(4,3,12), imagesc(squeeze(0.5*PT(:,round(end/2),:)),c)
end
