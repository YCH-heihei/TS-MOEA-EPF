function Output = g1(PopDec,n)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   m=size(PopDec,2);
   Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
   for i=m+1:n
       Output(:,i-m)=1/2*(1/m*sum(sin(1.5*pi*PopDec+j(i-m)),2)+1);
   end
end

