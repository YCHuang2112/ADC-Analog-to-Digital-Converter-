% x = [1; 1; 1; 1; 4; 4; 4; 4];
% y = [1; 2; 2; 1; 1; 2; 2; 1];
% z = [1; 1; 3; 3; 3; 3; 1; 1];
% c = [0; 0; 0; 0; 1; 1; 1; 1];
% % figure
% view(3);
% patch(x,y,z,c,'FaceColor','interp')
% colorbar
set(gcf, 'color', [1 1 1]);
v = [8 1 1; 8 1 3; 8 3 3; 8 3 1;
     1 1 1; 1 1 3; 1 3 3; 1 3 1;
    ];
% % % f = [1 2 3 4; 5 6 7 8; 1 2 5 6; 2 3 6 7; 3 4 7 8; 4 1 8 5;];
f = [1 2 3 4; 5 6 7 8; 1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8;];
col = [0;0;0;0; 8;8;8;8; ];
view(3);
patch('Faces',f,'Vertices',v,'FaceVertexCData',col,'FaceColor','interp');
colorbar;
x_axis_tick = 1:1:8; y_axis_tick = [1, 5, 10]; z_axis_tick = [1, 2, 3, 4];
xticks(x_axis_tick); yticks(y_axis_tick); zticks(z_axis_tick);
xlabel('???','FontSize',32); ylabel('***','FontSize',32);  zlabel('<*O*>','FontSize',32);  
xlim([0 8]);ylim([0 10]);zlim([0 4]);
process_code = 'Example';
title ({['cube (',process_code,')'];['Yu Heng Cubic (Bad guy)']},'FontSize',14);
set(gca, 'FontSize',13);
set(gcf, 'Position',  [488   342   580   420]);
legend('boxoff');
