function c1 = pcircle(x,y,radius,color)
c1 = rectangle('Position',[x-radius,y-radius,2*radius,2*radius],'Curvature',[1 1],'FaceColor',color,'EdgeColor',color,...
    'LineWidth',3);

end