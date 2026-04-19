function Output = g8(PopDec,n)
%  UNTITLED 此处显示有关此函数的摘要
%  此处显示详细说明
   m=size(PopDec,2);
   Output=zeros(size(PopDec,1),n-m);
   j=2*pi*((m+1:n)-(m+1))/n;
%    Output=repmat(sum(PopDec,2)-e(-1),1,n-m);
   for i=m+1:n
      Output(:,i-m)=1/m*sum(abs(sin(2.5*pi*(PopDec-0.5)+j(i-m))),2);
   end
end

