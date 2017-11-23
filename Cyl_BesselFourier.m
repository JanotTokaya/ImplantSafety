function [F] = Cyl_BesselFourier(Rcyl, PH, ka, M, normalize)
if nargin ==4
    normalize = 1;
end
% construct cylindrical Bessel-Fourier functions
F =[];
for order =-M:M
    func = besselj(order, ka*Rcyl).*exp(1i*order*PH);
    %func = func(:); func(index_zeros,:) = []; func = func/norm(func);
    F = [F func(:)];%[F func(:)/norm(func(:))];
end
if normalize ==1
    for j =1:size(F,2)
        normcolumn = norm(F(:,j));
        F(:,j) = F(:,j)/normcolumn;
    end
    % previous 4 lines replaces this: F = normc(F); changed on 07-05-2013
    % because normc today does not accept complex data (!)
end
