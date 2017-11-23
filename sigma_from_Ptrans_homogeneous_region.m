function [sigma_est1, sigma_est1l ,sigma_est1u ] = sigma_from_Ptrans_homogeneous_region(PT , X_PT, Y_PT, Z_PT)

mu0=1.2566370614*10^-6; 
eps0=8.85418782*10^-12;
f=63.87*10^6; 
omega = 2*pi*f; 
%% Determine conductivity:
dX = abs(mean(unique(diff(X_PT,1,1)))); 
dY = abs(mean(unique(diff(Y_PT,1,2)))); 
dZ = abs(mean(unique(diff(Z_PT,1,3)))); 
% L = del2(PT,dX,dY,dZ); 
% sigma_est2 = L/(2*mu0*omega); 
%% lets go 3D
s=size(PT);
XY_data = zeros(numel(PT),3);
XY_data(:,1) = X_PT(:);
XY_data(:,2) = Y_PT(:);
XY_data(:,3) = Z_PT(:);
                   
% options = optimoptions('lsqcurvefit'); 
% options.TolFun=1e-30;
% options.TolX=1e-30;
% options.MaxIter=2000;
% options.MaxFunEvals= 10000; 
% options.FinDiffRelStep=sqrt(eps/100);

fitted_phase = zeros(size(PT)); 
% poly2 = @(C,X) C(1)*X.^2-C(2);
% 
% poly2 = @(C,X) C(1)*X.^2-C(2);
% [C,~] = lsqcurvefit(poly2, [(min(X_mid)-max(X_mid))/s(1), round(s(1)/2)],(1:s(1)).',X_mid,[],[], options);
poly2 = 'a*(x-b)^2+c';
PT_midx = smooth(squeeze(PT(:,round(end/2),round(end/2))));

PT_midy = smooth(squeeze(PT(round(end/2),:,round(end/2))));

PT_midz = smooth(squeeze(PT(round(end/2),round(end/2),:)));
const = zeros(3,1);
h=figure;,
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
%[~, I] = min(smooth(squeeze(PT(:,round(end/2),round(end/2))))); 

C = [mean(const) I J K]; 

phi_of_xyz_no_linear_terms = @(C,XY_f)  C(1)...
                       +C(2)*(XY_f(:,1)-x_off).^2 ...
                       +C(3)*(XY_f(:,2)-y_off).^2 ...
                       +C(4)*(XY_f(:,3)-z_off).^2;
  
% 
% options.TolFun=1e-30;
% options.TolX=1e-30;
% options.MaxIter=1000;
% options.MaxFunEvals= 5000;                  
[C_est,~,resid,~,~,~,J] = lsqcurvefit(phi_of_xyz_no_linear_terms, C, XY_data, 0.5*double(PT(:)),[],[]); 

ci = nlparci(C_est,resid,'jacobian',J); 
%deltas = 0.5*diff(ci,1,2);

sigma_est1 = 2*(C_est(2)+C_est(3)+C_est(4))/(mu0*omega);
sigma_est1l = 2*(ci(2,1)+ci(3,1)+ci(4,1))/(mu0*omega);
sigma_est1u = 2*(ci(2,2)+ci(3,2)+ci(4,2))/(mu0*omega);

%sigma_unc = 2/(mu0*omega)*sqrt(deltas(2)^2+deltas(3)^2+deltas(4)^2);


fit_results_PT = reshape(phi_of_xyz_no_linear_terms(C_est,XY_data), s);

subplot(4,3,7), imagesc(fit_results_PT(:,:,round(end/2)))
c=caxis;
title('fitted phase')
subplot(4,3,10), imagesc(0.5*PT(:,:,round(end/2)),c)
title('transceive phase')
subplot(4,3,8), imagesc(squeeze(fit_results_PT(round(end/2),:,:)))
c=caxis;
subplot(4,3,11), imagesc(squeeze(0.5*PT(round(end/2),:,:)),c)
subplot(4,3,9), imagesc(squeeze(fit_results_PT(:,round(end/2),:)))
subplot(4,3,12), imagesc(squeeze(0.5*PT(:,round(end/2),:)))
