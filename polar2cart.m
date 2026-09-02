function [gridCart] = polar2cart(gridPolar)
    r = gridPolar.range;
    theta = gridPolar.azi;
    gridCart.x = r .* cos(theta);
    gridCart.y = r .* sin(theta);
end