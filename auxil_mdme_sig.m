function Mtd = auxil_mdme_sig(M0, T1, T2, Td, TE, TR, alpha, theta, b1scale)

    % M0 = 1;
    cos_theta = cos(theta * b1scale);
    cos_alpha = cos(alpha * b1scale);
    
    Td = Td(:);
    TE = TE(:);
    
    Mtd = M0 .*...
        (1 - (1 - cos_theta) .* exp( -Td'./T1 ) - cos_theta * exp( -TR/T1) ) ./...
        (1 - cos_theta * cos_alpha * exp( -TR/T1 ));
    
    Mtd = exp(-TE/T2) * Mtd;

end