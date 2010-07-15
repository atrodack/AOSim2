%% AOFIFO.  Support class for AOSim2.
% 20090425 JLCodona.

classdef AOFIFO < handle
    %AOFIFO A general delay line class for AOSim2.
    %  20090425 JLCodona.
    
    properties
        pipeline = {};
        iread = 1;
        iwrite = 1;
    
        name = '<name>';
    end
    
    methods
        function FIFO = AOFIFO(Length,name)
            FIFO.pipeline = cell(1,Length);
            if(nargin>1)
                FIFO.name = name;
            end
        end
    
        function FIFO = setDelay(FIFO,delay)
            if(delay>length(FIFO.pipeline))
                warning('AOFIFO:LONGER','Making the AOFIFO %s longer.',FIFO.name);
                FIFO.pipeline{delay} = [];
            end
            
            FIFO.pipeline{FIFO.iwrite} = data;
            FIFO.iread = FIFO.wrap(FIFO.iwrite-delay);
        end
        
        function FIFO = push(FIFO,data)
            FIFO.pipeline{FIFO.wrap(FIFO.iwrite)} = data;
            FIFO.iwrite = FIFO.iwrite + 1;
            FIFO.check
        end

        function data = peek(FIFO) % Read but don't increment.
            data = FIFO.pipeline{FIFO.wrap(FIFO.iread)};
        end
        
        function data = pop(FIFO)
            if(poll(FIFO))
                data = FIFO.pipeline{FIFO.wrap(FIFO.iread)};
                increment(FIFO);
            else
                    warning('AOFIFO:NOT_READY','Tried to read EMPTY AOFIFO %s.',FIFO.name);
            end
        end

        function  the_lag = lag(FIFO)
            the_lag = FIFO.iwrite - FIFO.iread;
        end

        function okay = check(FIFO)
            okay = FIFO.iread <= FIFO.iwrite ...
                && FIFO.iwrite - FIFO.iread < length(FIFO.pipeline);
        
            if(~okay)
                if(FIFO.iread > FIFO.iwrite) 
                    warning('AOFIFO:FUBAR','AOFIFO %s is EMPTY.',FIFO.name);
                elseif(FIFO.iwrite-FIFO.iread>length(FIFO.pipeline))
                    warning('AOFIFO:FUBAR','AOFIFO %s OVERFLOW.',FIFO.name);
                end
                
            end
        end
        
        function data_available = poll(FIFO)
            data_available = FIFO.iread < FIFO.iwrite;
        end
        
        function FIFO = increment(FIFO)
            FIFO.iread = FIFO.iread + 1;
            FIFO.check;
        end
        
        function FIFO = decrement(FIFO)
            FIFO.iread = FIFO.iread - 1;
            FIFO.check;
        end
        
        function FIFO = reset(FIFO)
            FIFO.iread = 1;
            FIFO.iwrite = 1;
        end
        
        function FIFO = flush(FIFO)
            FIFO.pipeline = cell(size(FIFO.pipeline));
            FIFO.iread = 1;
            FIFO.iwrite = 1;
        end
        
        function n = wrap(FIFO,n)
            L = length(FIFO.pipeline);
            n = mod(n-1,L)+1;
        end
    end
end

