% Convert a DCM to Euler angles
% function [roll, pitch, heading] = Dcm2Euler(dc)
function [roll, pitch, heading] = dcm2euler(dc)
pitch = atan(-dc(3,1)/sqrt(dc(3,2)^2 + dc(3,3)^2));
if dc(3,1) <= -0.999
    roll = NaN;
    heading = atan2((dc(2,3)-dc(1,2)),(dc(1,3)+dc(2,2)));
elseif dc(3,1) >= 0.999
    roll = NaN;
    heading = pi + atan2((dc(2,3)+dc(1,2)),(dc(1,3)-dc(2,2)));
else
    roll = atan2(dc(3,2), dc(3,3));
    heading = atan2(dc(2,1), dc(1,1));
end
