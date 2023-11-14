function  Show_all(I1,I2,I3,title1,title2,title3,link3together,main_title)
figure
ax(1)= subplottight(1,3,1), imshow(I1, 'border', 'tight'),title(num2str(title1))
ax(2)= subplottight(1,3,2), imshow(I2, 'border', 'tight'),title({['\fontsize{16}\color{red}',num2str(main_title)];['\fontsize{10}\color{black}',num2str(title2)]})
ax(3)=subplottight(1,3,3), imshow(I3, 'border', 'tight'),title(num2str(title3))
if(link3together)
    linkaxes([ax(1) ax(2) ax(3)],'xy')
else
   linkaxes([ax(1) ax(2)],'xy');
end

end

