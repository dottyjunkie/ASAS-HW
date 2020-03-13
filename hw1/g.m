function [output] = g(i,M)
    if i<=-M && i<=M
        output=cos(pi/2*i/M)^2
    else
        output=0
    end
end

