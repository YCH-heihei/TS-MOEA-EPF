function Output = Z5(PopDec,K)
   Output=-0.7*Z3(PopDec,K)+10/size(PopDec,2)*sum(abs(PopDec).^0.002,2);
    assert(all(all(Output<10)&all(Output>0)),'函数写错了');
end

