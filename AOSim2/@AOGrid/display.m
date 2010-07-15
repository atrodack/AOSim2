function display(G)
    fprintf('%s %s: %dx%d ',class(G),G.name,G.size);
    if(G.isX)
        fprintf('domain:SPACE ');
    else
        fprintf('domain:KAPPA ');
    end
    
    if(G.isCentered)
        fprintf('axis:CENTER ');
    else
        fprintf('axis:CORNER ');
    end
    
    fprintf('\n');
end
