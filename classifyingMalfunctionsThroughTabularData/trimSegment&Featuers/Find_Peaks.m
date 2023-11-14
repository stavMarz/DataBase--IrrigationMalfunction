function [UB_canopy,UB_background] = Find_Peaks(im)

[counts,binLocations]=hist(im(:),45);
filterd_idx= find(counts>5);
binLocations =  binLocations(filterd_idx);
counts =counts(filterd_idx);
[Y_peak,X_peak] = findpeaks(counts,binLocations,'SortStr','descend','npeaks',2);%Sort the peaks from tallest to shortest.
UB_canopy = min(X_peak)
UB_background = max(X_peak)

figure();
findpeaks(counts,binLocations,'SortStr','descend','npeaks',2)
text(X_peak+.5,Y_peak+5000,num2str((1:numel(Y_peak))'))
text(X_peak-1,Y_peak-15000,num2str(X_peak','%.2f'))
xlabel('Pixel Value')
ylabel('Frequency')
axis tight
grid on

end

