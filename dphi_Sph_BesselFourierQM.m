function [F_der] = dphi_Sph_BesselFourierQM(R, TH, PH, ka, M)
% make derivatives of spherical functions
Rv=R(:);    %radial coordinate
THv=TH(:);  %radial coordinate
PHv=PH(:);  %radial coordinate
kaR=ka*Rv;  %ka=sqrt(eps_r*mu0*eps0*(2*pi*omega)^2-1i*sigma*(2*pi*omega)*mu0);

[F_func, F_der] = deal(zeros(length(Rv),(M+1)^2));
ind = 0;
for n=0:M %the degree of the harmonic
    Pn=legendre(n,cos(THv),'norm');
    for m=-n:n %the order of the harmonic
        %spherical bessel function of the first kind jn, jn+1 and jn-1
        jn      = sqrt(pi./(2*kaR)).*besselj(n+(1/2),kaR);
        Pnm=Pn(abs(m)+1,:).';

        Ynm = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(4*pi*factorial(n+m)))*Pnm.*exp(1i*m*PHv);

        %derivative of the spherical harmonic w.r.t phi
        DYnmDphi=1i*m*Ynm;

        %the original function F
        func = jn.*Ynm;
        ind = ind+1;
        F_func(:,ind)=func;
        %the derivative (idx-dy) of F in spherical coordinates
        F_der(:,ind) =jn.*DYnmDphi; 
    end
end
end
