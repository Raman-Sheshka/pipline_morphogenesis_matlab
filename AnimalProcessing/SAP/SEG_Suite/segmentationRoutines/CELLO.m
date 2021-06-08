%% CELLO class definition
% GOYA Yûki
% v0.1 (last update 2012-11-13)
%
classdef CELLO < handle   
    % properties definition
    properties (SetAccess=private)
        tag='';
        category=0;
        track=[];
    end
    % methods definition
    methods
        % constructor
        function C = CELLO(celltag,cellcategory,celltrack)
            C.tag=celltag;
            C.category=cellcategory;
            C.track=celltrack;
        end
        % tag getter
        function t=getTag(C)
           t=C.tag; 
        end
        % category getter
        function ctgr=getCategory(C)
           ctgr=C.category; 
        end
        % track getter (all)
        function trck=getFullTrack(C)
           trck=C.track; 
        end
        % track getter (select)
        function trf=getTrackAt(C,frame)
            if frame>0 && frame<=numel(C.track)
                trf=C.track(frame);
            else
                trf=-1;
                disp('ERROR: OUT OF BOUND')
            end
        end
        % category setter
        function setCategory(C,newcategory)
            C.category=newcategory; 
        end
   end % methods
end

% end of file