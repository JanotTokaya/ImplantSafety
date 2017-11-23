function [F_der] = dx_Sph_BesselFourierQM(R, TH, PH, ka, M)
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
        %spherical harmonics
        Pnm=Pn(abs(m)+1,:).';
        Ynm = sqrt((2*n+1)*factorial(n-m)/(4*pi*factorial(n+m)))*Pnm.*exp(1i*m*PHv);         
        if m < n
            Pn_mp1 =Pn(abs(m+1)+1,:).';
            Yn_mp1 = sqrt((2*n+1)*factorial(n-(m+1))/(4*pi*factorial(n+(m+1))))*Pn_mp1.*exp(1i*(m+1)*PHv);            
        else 
            Yn_mp1 = zeros(size(THv));
        end
        
        %derivative of the spherical harmonic w.r.t theta
        DYnmDtheta=m*cot(THv).*Ynm+sqrt((n-m)*(n+m+1))*exp(-1i*PHv).*Yn_mp1;
        %both factors below are the same
        %sqrt(gamma(-m+n+1)*gamma(m+n+2)/(gamma(n-m)*gamma(n+m+1)))
        %sqrt((n-m)*(n+m+1))

        %the original function F
        func = jn.*Ynm;
        ind = ind+1;
        F_func(:,ind)=func;
        %the derivative (idx-dy) of F in spherical coordinates
        F_der(:,ind) =jn.*DYnmDtheta;
    end
end
end
