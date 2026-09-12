function dist = dist_nogeo(XI,XJ)
        X = abs(XI(:,1) - XJ(:,1));
        Y = abs(XI(:,2) - XJ(:,2));
        dist = sqrt(X.^2 + Y.^2);
end
