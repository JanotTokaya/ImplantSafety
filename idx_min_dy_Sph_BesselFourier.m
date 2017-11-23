function [F_attempt2] = idx_min_dy_Sph_BesselFourier(R, TH, PH, ka, M,normalize)
% make spherical functions basis for solutions of spherical Helmholtz equation
if nargin ==5
    normalize = 1;
end
Rv=R(:);
kaR=ka*Rv;
THv=TH(:);
PHv=PH(:);
[F_func, F] = deal(zeros(length(Rv),(M+1)^2));
ind = 0;
for n=0:M %the degree of the harmonic
    Pn=legendre(n,cos(THv),'norm');
    for m=-n:n %the order of the harmonic
        Pnm=Pn(abs(m)+1,:).';
        if m < n
            Pn_mp1 =Pn(abs(m+1)+1,:).';
        else 
            Pn_mp1 = zeros(size(THv));
        end
        %spherical bessel function of the first kind
        jn      = sqrt(pi./(2*kaR)).*besselj(n+(1/2),kaR);
        jn_p1   = sqrt(pi./(2*kaR)).*besselj(n+(3/2),kaR);
        jn_m1   = sqrt(pi./(2*kaR)).*besselj(n-(1/2),kaR);

        %derivative of the spherical bessel function of the first kind w.r.t
        %r
        DjnDr=n*jn./Rv-ka*jn_p1;
        % this alternative expression
        % DjnDr_alt=(jn_m1-(n+1)./kaR.*jn)*ka;
        % gives the same result. It is an indentity of the equation above.
        
        %spherical harmonics
        Ynm=Pnm.*exp(1i*m*PHv); 
        Yn_mp1 = Pn_mp1.*exp(1i*(m+1).*PHv);
        
        %derivative of the spherical harmonic w.r.t theta
        DYnmDtheta=m*cot(THv).*Ynm+sqrt((n-m)*(n+m+1))*exp(-1i*PHv).*Yn_mp1;
        %derivative of the spherical harmonic w.r.t phi
        DYnmDphi=1i*m*Ynm;

        %the original function F
        func = jn.*Ynm;
        ind = ind+1;
        F_func(:,ind)=func;
        %the derivative (idx-dy) of F
        F_attempt2(:,ind) =1i*exp(1i*PHv).*sin(THv).*DjnDr.*Ynm...
            -exp(1i*PHv)./(Rv.*sin(THv)).*jn.*DYnmDtheta+...
            1i*exp(1i*PHv).*cos(THv)./Rv.*jn.*DYnmDphi; 
        F_attempt2(:,ind) =1i*exp(-1i*THv).*sin(PHv).*DjnDr.*Ynm...
            +exp(-1i*THv)./(Rv.*sin(PHv)).*jn.*DYnmDtheta+...
            1i*exp(-1i*THv).*cos(PHv)./Rv.*jn.*DYnmDphi; 
    end
end

if normalize ==1
    for j =1:size(F,2)
        normcolumn = norm(F_func(:,j));
        if normcolumn == ~0;
            F(:,j) = F(:,j)/normcolumn;
            F_attempt2(:,j) = F_attempt2(:,j)/normcolumn; 
        end
    end
    % previous 4 lines replaces this: F = normc(F); changed on 07-05-2013
    % because normc today does not accept complex data (!);
end
