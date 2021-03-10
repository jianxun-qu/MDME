function sig = auxil_mdme_ir(M0, T1, T2, TE, TR, TI)

    sig = M0 * (1-2*exp(-TI/T1)+exp(-TR/T1)) * exp(-TE/T2);

end