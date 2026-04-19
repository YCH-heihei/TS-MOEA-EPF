function Output = g9(PopDec,n)
%  UNTITLED 此处显示有关此函数的摘要
%  此处显示详细说明
   m=size(PopDec,2);
%    Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
   Output=repmat(1/m*sum(PopDec,2),1,n-m);
   for i=m+1:n
      Output(:,i-m)=1/2*(Output(:,i-m)-1/m*sum(abs(sin(2.5*pi*PopDec-1/2*pi+j(i-m))),2)+1);
   end
end

