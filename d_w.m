function d=d_w(N_MS,w)
d_w=zeros(N_MS,1);
for i=1:1:N_MS
    d_w(i)=exp(1j*(i-1)*w);
end
d=1./sqrt(N_MS)*d_w;
