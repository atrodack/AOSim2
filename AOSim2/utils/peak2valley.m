function delta = peak2valley(V)

% delta = peak2valley(V)

delta = max(V(:)) - min(V(:));
