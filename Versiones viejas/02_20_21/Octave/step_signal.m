% Help
function [y]= step_signal(time, a, offset, delay);
    %step_signal(t,
%     time=0:0.1:1;
%     a=0.1; offset=0.3;delay=-0.3
    y=zeros(size(time));
    y=((time+delay)>=0);
    y= y*a + offset;
end