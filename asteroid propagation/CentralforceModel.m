%%
function dot_x  =CentralforceModel(x,muSun)

vel = x(1:3);
acc = -muSun/r^3*pos;


dot_x = [vel,acc];

return
end