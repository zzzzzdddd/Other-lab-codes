for i = 2:40
    fname = char(strcat("leave_one_out",string(i),"n.m"))
    copyfile('leave_one_out1n.m', fname)
end