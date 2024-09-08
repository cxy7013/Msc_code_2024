function lens_phase=lens_focus(D,N,wavl,f)
delta=D/N;
n=1;
k=n*2*pi/wavl;
x=(-N/2:N/2-1)*delta;
y=x;
[X,Y]=meshgrid(x,y);
[~,r]=cart2pol(X,Y);
lens_phase=exp(-1i*k/2/f*r.^2);
end
