function [threshold] = FWHM(im,lev,limit_of_left_gauss)
 
    [counts,centers]= hist(im,255);
    positiv_indices = find(centers>=0);%find the positive values of x
    sub_counts =counts(positiv_indices);
    sub_centers =centers(positiv_indices);
%     figure();
%     plot(sub_centers,sub_counts);
%     title("Histogram of the plot area")
%     xlabel('Temperature (c\cdot)')
%     ylabel('Frequency')
    first_gauss_indices = find(sub_centers<limit_of_left_gauss);
    max_count =  max(sub_counts(first_gauss_indices));
    
    y = sub_counts /max_count;%normalization of the all counters in ratio to the max counter of  the first gaussian
    N = length(y);
    level =lev;
    x= sub_centers;
 
    if y(1) < level % find index of center (max or min) of pulse
        [garbage,centerindex]= max(sub_counts(first_gauss_indices));
%         disp('centerindex=')
        disp(centerindex)
%         Pol = +1;
%         disp('Pulse Polarity = Positive')
    else
        [garbage,centerindex]=min(y);
%         Pol = -1;
%         disp('Pulse Polarity = Negative')
    end
    i = 2;
 
    while sign(y(i)-level) == sign(y(i-1)-level) %get up until the level of max
        i = i+1;
    end
%      disp('end of get up  half at while and  new i=')
     disp(sign(y(i)-level))
     disp(sign(y(i-1)-level))
     disp(i-1)
 
    %first crossing is between v(i-1) & v(i)
    interp = (level-y(i-1)) / (y(i)-y(i-1));
%     disp('interp=')
%     disp(interp)
       tlead = x(i-1) + interp*(x(i)-x(i-1));
%     disp('F(x1)=')
%     disp(sub_counts(i-1) + interp*(sub_counts(i)-sub_counts(i-1)))
%     disp('tlead (x1)=')
%     disp(tlead)
    i = centerindex+1; %jump to a point after the peak   
 
     %start search for next crossing at center
    while ((sign(y(i)-level) == sign(y(i-1)-level)) & (i <= N-1))
        i = i+1;
    end
%         disp('stop search for next crossing at center')
    if i ~= N
%         disp('i!=N')
%         Ptype = 1;  
%         disp('Pulse is Impulse or Rectangular with 2 edges')
        interp = (level-y(i-1)) / (y(i)-y(i-1));
        ttrail = x(i-1) + interp*(x(i)-x(i-1));
%         disp('F(x2)=')
%         disp(sub_counts(i-1) + 0.5*interp*(sub_counts(i)-sub_counts(i-1)))
         disp('ttrail (x2)=')
         disp(ttrail)
%         width = ttrail - tlead;
    else
%         disp('enter else')
%         Ptype = 2; 
%         disp('Step-Like Pulse, no second edge')
        ttrail = NaN;
        width = NaN;
    end
%     new_thermal = zeros(size(im));
%      new_thermal(im<tlead)=0; 
%    new_thermal(find(im<ttrail))=0;
%     final_im = new_thermal; 
    threshold = ttrail;
 end

