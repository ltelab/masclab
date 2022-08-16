% approxime the perimeter of an ellipse of semi grand-axe a and semi
% petit-axe b using Ramanuja formula
% if a and b are same size vectors, the function computes the corresponding vector of
% perimeters
function p = ellipse_perim(a,b)

    p = pi.*(3.*(a+b) - sqrt((3.*a+b).*(a+3.*b)));

end
