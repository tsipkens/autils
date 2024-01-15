
% LOGSPACE10  A simple bridging function to avoid log10 use in logpsace function.
%  AUTHOR: Timothy Sipkens, 2024-01-15

function v = logspace10(a, b, n)

v = logspace(log10(a), log10(b), n);

end
