function area = get_area(x, y)
% polygnoal area
% https://zhuanlan.zhihu.com/p/110025234
area = 0;

for k = 1:length(x)
    if k == 1
        x0 = x(end); y0 = y(end);
        x1 = x(1); y1 = y(1);
    else
        x0 = x(k-1); y0 = y(k-1);
        x1 = x(k); y1 = y(k);
    end
    area = area + x0*y1 - x1*y0;
end

area = 0.5 * abs(area);
end