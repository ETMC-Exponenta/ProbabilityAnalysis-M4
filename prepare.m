path_array = { 'Ver'; 'NewLect'; 'My_Tools'; 'WorkProb' };

addpath('.');

[row, ~] = size( path_array );
for i = 1:row                    
	addpath( ['./' path_array{i,:} ] );
end