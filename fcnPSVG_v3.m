function o=fcnPSVG_v3(y,opt)
% y is the time series o is the negative of the log-log slope
% opt=1 plots the log-log slope and point values
if nargin<2
    opt=0;
end        
ln=length(y);
P=zeros(1,ln);
for a=1:ln
    for b=a+2:ln
         fl=1;
         for c=a+1:b-1
            if y(c)>=y(b)+(y(a)-y(b))*(b-c)/(b-a)
                fl=0;
                break;
            end
         end
         if fl %add one to the graph edges of length b-a
             P(b-a)=P(b-a)+1;
         end
    end
end
kk=1:ln;
%loglog(kk,P); %do a log-log plot of the distribution
if opt==1
    plot(log(kk),log(P),'.');
end
for x=2:ln %choose an interval in P with no zeros
    if P(x)==0
        break;
    end
end
p=polyfit(log(2:x-1),log(P(2:x-1)),1); %find the power law slope
if opt==1
    hold on
    x=[log(kk(1)):0.1:log(kk(length(kk)))];
    plot(x,p(1)*x+p(2),'r');
    text(2,9,['Slope=' num2str(p(1))]);
    hold off
    fprintf('Slope=%4.4f\n',p(1));
end
o=-p(1);