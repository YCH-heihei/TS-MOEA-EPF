function Output = Z1(PopDec,K)
   Output=10/size(PopDec,2)*sum(PopDec.^2,2);
   assert((all(Output<10)& all(Output>0)),'函数写错了');
end

