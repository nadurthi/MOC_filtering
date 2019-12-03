function f = FirstSignificant(a)
f = fix(abs(a) .* 10 .^ (-floor(log10(abs(a)))));
end