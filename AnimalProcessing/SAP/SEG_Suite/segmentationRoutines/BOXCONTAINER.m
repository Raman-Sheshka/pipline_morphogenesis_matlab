classdef BOXCONTAINER < handle
    %% boxcontainer class definition
    % GOYA Yûki
    % 2012-03-15
    
    %% propreties
    properties %(GetAccess=private)
        Boxes       % list of boxes
        Slot        % next slot to fill
        Capacity    % container capacity
    end
    properties (Dependent)
       UsedBoxes    % number of boxes used
    end
    %% methods
    methods
        % constructor
        function obj=BOXCONTAINER(capa)
            obj.Capacity=capa;
            obj.Slot=1;
            obj.Boxes=zeros(capa,4);
        end
        % access to boxes
        function box=getBox(obj,i)
            if i>0 && i<=obj.Capacity
                box=obj.Boxes(i,:); 
            else
                box=[];
            end
        end
        % access to capacity
        function cap=getCapacity(obj)
            cap=obj.Capacity;
        end
        % expand capacity
        function expandCapacity(obj,newcapa)
            if obj.Capacity < newcapa
               % newBoxes=zeros(newcapa,4);
                if sum(obj.Boxes(obj.Slot))==0 % empty box
                    % expand boxes
                    obj.Boxes=[obj.Boxes ; zeros(newcapa-obj.Capacity,4)];
                    %set new capacity
                    obj.Capacity=newcapa;    
                else % already full
                    %resort the content
                    obj.Boxes=[obj.Boxes(obj.Slot:end,:) ; obj.Boxes(1:obj.Slot-1,:) ; zeros(newcapa-obj.Capacity,4)];                    
                    % update Slot
                    obj.Slot=obj.Capacity+1;
                    % set new capacity
                    obj.Capacity=newcapa;
                end
            else
                disp('New capacity must be higher than current capacity!');
            end
        end
        % add a box
        function addBox(obj,box)
            obj.Boxes(obj.Slot,:)=box;
            obj.Slot=obj.Slot+1;
            % loop Slot if out of bound
            if obj.Slot>obj.Capacity, obj.Slot=1; end
        end
        % erase a box
        function eraseBox(obj,i)
            if i>0 && i<=obj.Capacity
                obj.Boxes(i,:)=[];
                obj.Boxes=[obj.Boxes ; zeros(1,4)];
                % correct Slot
                obj.Slot=obj.Slot-1;
                % loop if out of bound
                if obj.Slot==0, obj.Slot=obj.Capacity; end
            end
        end
        % access to UsedBoxes
        function UsedBoxes=get.UsedBoxes(obj)
            UsedBoxes=numel(find(obj.Boxes(:,1)~=0));
        end
        % is empty?
        function b=isempty(obj)
            b=obj.Boxes(1,1)==0;
        end
    end
end