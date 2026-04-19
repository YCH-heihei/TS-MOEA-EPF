function Output = g6(PopDec,n)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   m=size(PopDec,2);
   denominator=-10*exp(-2/5)-exp(-1)+10+exp(1);
%    Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
   Output=repmat(-10*exp(-2/5*(1/m*sum(PopDec.^2,2)).^(1/2))+10+exp(1),1,n-m);
   for i=m+1:n
      Output(:,i-m)=Output(:,i-m)-exp(1/m*sum(cos(11*pi*PopDec+j(i-m)).^(3),2));
   end
   Output=Output./denominator;
end

