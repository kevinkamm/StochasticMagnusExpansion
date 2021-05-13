T=10;
t0=0;
N0=11;
N1=1001;
Nlist = findN(N0,N1);
for i=1:1:length(Nlist)
    N=Nlist(i);
    ind_t=[1:1:N];
    ind_t(2:1:end)=ind_t(1:1:end-1).*floor((N1-1)/(N-1))+1;

    tref=linspace(t0,T,N1);
    dtref=diff(tref);
    t=linspace(t0,T,N);
    dt=diff(t);
    fprintf('diff of t %g\n',sum(abs(squeeze(tref(ind_t))-t),'all')) 
end

function N=findN(N0,N1)
    k=1;
    N=-1;
    for i=N0:1:N1
        if mod((N1-1),(i-1))==0
            N(k)=i;
            k=k+1;
        end
    end
end