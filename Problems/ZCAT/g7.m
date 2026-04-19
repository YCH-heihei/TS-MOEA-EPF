function Output = g7(PopDec,n)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   m=size(PopDec,2);
   denominator=1+exp(1)-exp(-1);
%    Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
   Output=repmat(1/m*sum(PopDec,2)-exp(-1),1,n-m);
   for i=m+1:n
      Output(:,i-m)=Output(:,i-m)+exp(sin(sum(7*pi*1/m*sum(PopDec,2)-pi/2+j(i-m))));
   end
   Output=Output./denominator;
end

