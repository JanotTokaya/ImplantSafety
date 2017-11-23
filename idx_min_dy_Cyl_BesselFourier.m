function [F_der] = idx_min_dy_Cyl_BesselFourier(Rcyl, PH, ka, M,normalize)
if nargin ==4
    normalize = 1;
end
% construct cylindrical Bessel-Fourier functions
F_der =[];
F_func=[];
for order =-M:M
    func = besselj(order, ka*Rcyl).*exp(1i*order*PH);
    func_der = 1i*exp(1i*(order+1)*PH)*0.5*ka.*(besselj(order-1, ka*Rcyl)-besselj(order+1, ka*Rcyl))-1i*order*exp(1i*(order+1)*PH).*besselj(order, ka*Rcyl)./(Rcyl);
    %func = func(:); func(index_zeros,:) = []; func = func/norm(func);
    F_der = [F_der func_der(:)];%[F func(:)/norm(func(:))];
    F_func = [F_func func(:)];
end

if normalize ==1
    for j =1:size(F_der,2)
        normcolumn = norm(F_func(:,j));
        F_der(:,j) = F_der(:,j)/normcolumn;
    end
    % previous 4 lines replaces this: F = normc(F); changed on 07-05-2013
    % because normc today does not accept complex data (!)
end
 
end
