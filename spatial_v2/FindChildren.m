function [child, sister] = FindChildren(robot)
num = length(robot.parent);

for n=1:num
    robot.sister(n) = 0;
    robot.child(n)  = 0;
end

for n=1:num
    parent = robot.parent(n);  % My parent
    if parent ~= 0
        if robot.child(parent) == 0
            robot.child(parent) = n;  % I am the first child.
        else
            eldest_sister = robot.child(parent);  % I have a sister.
            SetYoungestSister(robot,eldest_sister,n);
        end
    end
end

child = robot.child;
sister = robot.sister;

end