function yout = step_multiple (t,t1,a1,offset) 
         
%    t1=[20 40 60]; a1=[0.8 0.6 1];
   y=ones(size(t));
   yout=y;
   tm=t(2)-t(1);
   t2=t1./tm;
   yout(1:t2(1))=y(1:t2(1)).*offset;
   yout(t2(1):t2(2))=y(t2(1):t2(2)).*a1(1);
   yout(t2(2):t2(3))=y(t2(2):t2(3)).*a1(2);
   yout(t2(3):end)=y(t2(3):end).*a1(3);
end