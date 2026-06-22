%--------------------------------------------------------------------------
function ansd = vett(r1,r2)  
%
%subfunction that evaluates vector product
ansd(1)=(r1(2)*r2(3)-r1(3)*r2(2));
ansd(2)=(r1(3)*r2(1)-r1(1)*r2(3));
ansd(3)=(r1(1)*r2(2)-r1(2)*r2(1));

