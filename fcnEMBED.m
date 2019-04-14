function [fdo, pko, fno, do]=fcnEMBED(iv,uf,tt,prt)
% iv time series uf calc fnns 
% e.g. embed(iv,0,0);
stt=2;
fin=12;
fd=zeros(1,fin);
pk=fd;
fn=fd;
if tt
    tic
end
for d=stt:fin
    [fd(d),pk(d),fn(d)]=fcnCD_PK_v2(iv,d,0,0,6,0,uf);
    if prt
        fprintf('%d %4.4f %4.4f\n',d,fd(d),fn(d));
    end
    if (uf==1 && fn(d)<=0.005)
        break;
    end
    if  fd(d)<1.01*fd(d-1)
        break;
    end
end
if fd(d)<fd(d-1)
        d=d-1;
end
fdo=fd(d);
pko=pk(d);
fno=fn(d);
do=d;
if prt
    fprintf('%d %4.4f\n',d,fdo);
end
if tt
%     fprintf('%4.2f %d\n',fdo,do);
    toc
end
