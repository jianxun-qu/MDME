function sysig = auxil_mdme_syflair(ti, te, tr, m0, t1, t2)

    
    sysig = m0 .* (1 - 2 * exp(-ti./t1) + exp(-tr./t1)) .* exp(-te./t2);


end