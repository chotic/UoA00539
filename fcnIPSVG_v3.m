function o=fcnIPSVG_v3(X,uplim,opt)
%X is the time series, uplim is max k (donwsampling average parameter)
% opt=1 prints and plots comparisons of psvg and ipsvg for each k
%o indicates average either based on convergence or failing that
%a simple grand average over k
epsilon=0.1;
if nargin<2
    uplim=8;
end
if nargin<3
    opt=0;
end
N=length(X);
v=zeros(1,uplim);
for k=1:uplim
    if opt
        fprintf('%d ',k);
    end
    for m=1:k
        Xmk=zeros(1,floor((N-m)/k)+1);
        for p=0:floor((N-m)/k)
            Xmk(p+1)=X(m+p*k); % mqke our x1, x3, x5, ... x2, x4, x6, ... etc
        end
        vv=fcnPSVG_v3(Xmk);
        if opt
            fprintf('%4.3f ',vv);
        end
        v(k)=v(k)+vv/k; %averages the k values for a given k
    end
    if opt
        fprintf(': %4.3f\n',v(m));
    end
    if k>1 && abs(v(k-1)-v(k))<epsilon;
        o=v(k);
        break;
    end
    if k>1 && k==uplim
        if abs(v(k-1)-v(k))>=epsilon
            o=mean(v);
            break;
        end
    end
end
if opt
    plot(v);
end


        
            
            
            
