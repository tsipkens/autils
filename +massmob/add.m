
% ADD  Add/update the mass-mobility parameters in an existing prop structure. 
%  
%  For example, this utility can add the mass-mobility information to a
%  structure containing flow rates, temperature, presssure, and other
%  parameters for evaluating the tranfer functions of classifiers. 
%  
%  PROP = massmob.add(PROP, F1, V1, F2, V2) uses the name-value pairs in
%  the form of (F1, V1) and (F2, V2) to add mass-mobility parameters to the
%  input PROP.
%  
%  PROP = massmob.add(PROP, TYPE) uses the particle type (e.g., 'soot' or
%  'salt') to add the mass-mobility parameters to the input PROP.
%  
%  NOTE: Overwrites existing mass-mobility information that may already be
%  in the structure. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function prop = add(prop, f1, v1, f2, v2)

% Remove the relevant fields, to be replaced.
fields = {'zet', 'Dm', 'm0', 'k', 'rho0', 'm100', 'rho100'};
for ii=1:length(fields)
    if isfield(prop, fields{ii})
        prop = rmfield(prop, fields{ii});
    end
end

if nargin == 2  % then second argument is string of particle type
    prop1 = massmob.init(f1);  % get properties for that string
    prop.zet = prop1.zet;  % copy into prop
    prop.rho100 = prop1.rho100;

else  % otherwise name-value pairs
    % Add new values.
    prop.(f1) = v1;
    prop.(f2) = v2;
end

% Fill out remaining values.
prop = massmob.init(prop);

end

