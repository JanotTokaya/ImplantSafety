function [F] = Sph_BesselFourierQM(R, TH, PH, ka, M)
% make spherical functions basis for solutions of spherical Helmholtz equation
Rv=R(:);
kaR=ka*Rv;
THv=TH(:);
PHv=PH(:);
F=zeros(length(Rv),(M+1)^2);
ind = 0;
for n=0:M
    Pn=legendre(n,cos(THv),'norm');
    for m=-n:n
        Pnm=Pn(abs(m)+1,:).';
        jn=besselj(n+1/2,kaR).*sqrt(pi./(2*kaR));
        ind = ind+1;
        %Ynm = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(4*pi*factorial(n+m)))*Pnm.*exp(1i*m*PHv); 
        Ynm = (-1)^m*sqrt((2*n+1)*factorial(n-m)/(4*pi*factorial(n+m)))*Pnm.*exp(1i*m*PHv); 
        %This normalization makes the harmonics orthonormal ans is taken
        %from quantum mechanics.
        func = jn.*Ynm;
        F(:,ind)=func;
    end
end

end
