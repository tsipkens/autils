
% UPDATE  A simple utility to update the mass-mobility parameters
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function prop = update(prop, f1, v1, f2, v2)

% Remove the relevant fields, to be replaced.
prop = rmfield(prop, {'zet', 'Dm', 'm0', 'rho0', 'm100', 'rho100'});

% Add new values.
prop.(f1) = v1;
prop.(f2) = v2;

% Fill out remaining values.
prop = massmob.gen(prop);

end

