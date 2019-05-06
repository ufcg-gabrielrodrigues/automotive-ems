function fit_sym = fitToSym(fit)

coeff_names = coeffnames(fit);
coeff_values = coeffvalues(fit);

fit_sym = str2sym(formula(fit));

for coeff_index = 1:numcoeffs(fit)
    fit_sym = subs(fit_sym, coeff_names{coeff_index}, num2str(coeff_values(coeff_index)));
end

end

