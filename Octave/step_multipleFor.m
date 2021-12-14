function yout = step_multipleFor (t,t1,a1) 
  %t1=[10 20 30 40 100 120]; a1=[1 2 3 4 5 6 7];
    
    y=ones(size(t));
    yout=ones(size(t))*a1(length(a1));
    tm=t(2)-t(1);
    t2(1)=1;
    t2(2:length(t1)+1)=t1./tm;
    for i=1:length(t1)
      yout(t2(i):t2(i+1))=y(t2(i):t2(i+1)).*a1(i);
    end
end