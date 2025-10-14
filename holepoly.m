function domain = holepoly

t = 0.05:0.5:2*pi;
x1 = cos(t);
y1 = sin(t);
x2 = 0.5*cos(t);
y2 = 0.5*sin(t);
domain = polyshape({x1,x2},{y1,y2});

end