function [F] = Sph_BesselFourier(R, TH, PH, ka, M,normalize)
% make spherical functions basis for solutions of spherical Helmholtz equation
if nargin ==5
    normalize = 1;
end
Rv=R(:);
kaR=ka*Rv;
THv=TH(:);
PHv=PH(:);
F=zeros(length(Rv),(M+1)^2);
ind = 0;
for n=0:M
    Pn=legendre(n,cos(THv));
    for m=-n:n
        Pnm=Pn(abs(m)+1,:).';
        jn=besselj(n+1/2,kaR).*sqrt(pi./(2*kaR));
        ind = ind+1;
        func = jn.*Pnm.*exp(1i*m*PHv);
        F(:,ind)=func;
    end
end
if normalize ==1
    for j =1:size(F,2)
        normcolumn = norm(F(:,j));
        if normcolumn == ~0;
            F(:,j) = F(:,j)/normcolumn;
        end
    end
    % previous 4 lines replaces this: F = normc(F); changed on 07-05-2013
    % because normc today does not accept complex data (!);
end
