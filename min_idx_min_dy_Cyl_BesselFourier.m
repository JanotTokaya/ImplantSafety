function [F] = min_idx_min_dy_Cyl_BesselFourier(Rcyl, PH, ka, M)
if nargin ==4
    normalize = 1;
end
% construct cylindrical Bessel-Fourier functions
F =[];
for order =-M:M
    func = besselj(order, ka*Rcyl).*exp(1i*order*PH);
    func_der = -1i*exp(1i*(order-1)*PH)*0.5*ka.*(besselj(order-1, ka*Rcyl)-besselj(order+1, ka*Rcyl))-1i*order*exp(1i*(order-1)*PH).*besselj(order, ka*Rcyl)./(Rcyl);
    %func = func(:); func(index_zeros,:) = []; func = func/norm(func);
    F = [F func_der(:)];%[F func(:)/norm(func(:))];
end
% if normalize ==1
%     for j =1:size(F,2)
%         normcolumn = norm(F(:,j));
%         if normcolumn == ~0;
%             F(:,j) = F(:,j)/normcolumn;
%         end
%     end
%     % previous 4 lines replaces this: F = normc(F); changed on 07-05-2013
%     % because normc today does not accept complex data (!)
 end
