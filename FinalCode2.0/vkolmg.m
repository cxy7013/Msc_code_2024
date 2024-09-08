function phz=vkolmg(D,dz,N,CN,wvl)
    % D phase screen size
    % N sampling number
    % dz propagation distance
    L0=1;% outer scale
    l0=0.002;% inner scale
    %D=1;
    delta=D/N;
    x=(-N/2:N/2-1)*delta;
    y=x;
    [X,Y]=meshgrid(x,y);
    del_f=1/(N*delta);
    fx=(-N/2:N/2-1)*del_f;
    [kx,ky]=meshgrid(2*pi*fx);
    k=2*pi/wvl; %wavenumber
    [~,ka]=cart2pol(kx,ky);%convert to polar
    km=5.92/l0;
    k0=2*pi/L0;
    
    %%
    PSD_phi=0.033*CN*exp(-(ka/km).^2)./(ka.^2+k0^2).^(11/6);%%spectrum
    cn=2*pi*k.^2*dz.*PSD_phi*(2*pi*del_f)^2;
    cn=fftshift(cn);
    phz_hi=N^2.*ifft2((randn(N)+1i.*randn(N)).*sqrt(cn));
    phz_hi=real(phz_hi);
    %phz=real(phz_hi);
    % low frequency compensation
    phz_lo=zeros(size(phz_hi));
    for p=1:3
        del_fp=1/(3^p*D);
        fx1=(-1:1)*del_fp;
        [kx1,ky1]=meshgrid(2*pi*fx1);
        [~,k1]=cart2pol(kx1,ky1);
        PSD_phi1=0.033*CN*exp(-(k1/km).^2)./(k1.^2+k0^2).^(11/6);
        PSD_phi1(2,2)=0;
        %random draws of Fourier coefficient
        cn1=2*pi*k.^2*dz.*PSD_phi1*(2*pi*del_fp).^2;
        cn1=fftshift(cn1);
        cn1=(randn(3)+1i*randn(3)).*sqrt(cn1);
        SH=zeros(N);
        for ii=1:9
            SH=SH+cn1(ii)*exp(1i*(kx1(ii)*X+ky1(ii)*Y));
        end
        phz_lo=phz_lo+SH;
    end
    phz=real(phz_hi+phz_lo);
end