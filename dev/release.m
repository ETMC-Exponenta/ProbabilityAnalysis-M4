gendoc;
v = input(sprintf('Specify version: (%s) ', dev.gvp), 's');
if isempty(v)
    v = dev.gvp;
end
dev.release(v);