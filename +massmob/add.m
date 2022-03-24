
% UPDATE  Add/update the mass-mobility parameters in an existing prop structure. 
%  
%  For example, this utility can add the mass-mobility information to a
%  structure containing flow rates, temperature, presssure, and other
%  parameters for evaluating the tranfer functions of classifiers. 
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function prop = update(prop, f1, v1, f2, v2)

% Remove the relevant fields, to be replaced.
prop = rmfield(prop, {'zet', 'Dm', 'm0', 'rho0', 'm100', 'rho100'});

% Add new values.
prop.(f1) = v1;
prop.(f2) = v2;

% Fill out remaining values.
prop = massmob.init(prop);

end

