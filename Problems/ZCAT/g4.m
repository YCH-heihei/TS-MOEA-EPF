function Output = g4(PopDec,n)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
   m=size(PopDec,2);
   Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
   for i=m+1:n
%        try
            Output(:,i-m)=1/2*(1/m*sum(PopDec.*cos(4*pi/m*sum(PopDec,2)+j(i-m)),2)+1);
%        catch
%            temp=0;
%        end
   end
end

