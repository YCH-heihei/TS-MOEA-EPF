function Output = Z2(PopDec,K)
   Output=10*max(abs(PopDec),[],2);
   assert(find(any(Output>10)||any(Output<10)),'函数写错了');
end

