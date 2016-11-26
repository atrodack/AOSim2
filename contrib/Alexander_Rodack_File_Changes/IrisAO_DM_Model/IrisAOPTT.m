function [ DM ] = IrisAOPTT( DM, segList, pistonList, tipList, tiltList, bump)
%**************************************************************************
%*********************************ANTIQUATED*******************************
%**************************************************************************
%****USE PTT METHOD IN AOAPERTURE INSTEAD---MUST MAP YOURSELF IF YOU DO****
%**************************************************************************
% [ DM ] = IrisAOPTT( DM, segList, pistonList, tipList, tiltList )
%   Applys Piston, Tip, and Tilt to an IrisAO model. This will do it in
%   hardware order, not in the Model's order. Input segList can be 1:37,
%   and this will map to the corresponding segment in DM.segList.
%
% The file PTT111 Numbered Array in the folder is a labeled picture of the
% hardware order of segments. Use this to assist you learn how the segments
% are layed out (rings outward from the center).
%
% pistonList is in units of meters---> on order of microns
% tipList and tiltList are in units of radians ----> on order of milliradians

%% The Ugly Nargin Section....
%all the fprintf's are commented because they take too much time to print
%on screen. They are left here for debugging purposes or so they can be
%turned on if desired

if nargin == 1
    segList = 1:37;
    pistonList = zeros(37,1);
    tipList = zeros(37,1);
    tiltList = zeros(37,1);
    bump = false;
%     fprintf('Flattening the Mirror\n');
elseif nargin == 3
    tipList = zeros(37,1);
    tiltList = zeros(37,1);
    bump = false;
%     fprintf('Applying Piston Only\n');
elseif nargin == 4
    tiltList = zeros(37,1);
    bump = false;
%     fprintf('Applying Piston and Tip\n');
elseif nargin == 5
%     maxval1 = max(abs(pistonList));
%     maxval2 = max(abs(tipList));
%     maxval3 = max(abs(tiltList));
%     if maxval1 <= 1e-10
%         if maxval2 == 0
%             if maxval3 == 0
%                 fprintf('Flattening the Mirror\n')
%             else
%                 fprintf('Applying Piston and Tilt\n');
%             end
%         else
%             if maxval3 == 0
%                 fprintf('Applying Piston and Tip\n');
%             else
%                 fprintf('Applying Piston, Tip, and Tilt\n');
%             end
%         end
%         fprintf('Applying Piston\n');
%     end
    
    bump = false;
elseif nargin == 6
%     maxval1 = max(abs(pistonList));
%     maxval2 = max(abs(tipList));
%     maxval3 = max(abs(tiltList));
%     if maxval1 <= 1e-10
%         if maxval2 == 0
%             if maxval3 == 0
%                 fprintf('Flattening the Mirror\n')
%             else
%                 fprintf('Applying Tilt\n');
%             end
%         else
%             if maxval3 == 0
%                 fprintf('Applying Tip\n');
%             else
%                 fprintf('Applying Tip and Tilt\n');
%             end
%         end
%     else
%         if maxval2 == 0
%             if maxval3 == 0
%                 fprintf('Applying Piston\n')
%             else
%                 fprintf('Applying Piston and Tilt\n');
%             end
%         else
%             if maxval3 == 0
%                 fprintf('Applying Piston and Tip\n');
%             else
%                 fprintf('Applying Piston, Tip, and Tilt\n');
%             end
%         end
%     end
%     if bump == true
%         fprintf('*****Bumping Segments*****\n');
%     else
%         fprintf('*****Overwriting Piston, Tip, and Tilt*****\n');
%     end
else
    error('Number of input arguments is incorrect');
end

if ~isa(DM,'AOAperture')
    error('DM must be AOSim2 class AOAperture');
end

%% Load in the Mapping Vector
load('IrisAO_SegMap.mat');
%% Do the Mapping and Apply PTT
if bump == false
    for ii = 1:37
        mapped_segment = IrisAO_SegMap(ii);
        DM.segList{ii}.piston = pistonList(mapped_segment);
        DM.segList{ii}.tiptilt = [tipList(mapped_segment),tiltList(mapped_segment)];
    end
else
    for ii = 1:37
        mapped_segment = find(IrisAO_MAP == ii);
        DM.segList{ii}.piston = DM.segList{ii}.piston + pistonList(mapped_segment);
        DM.segList{ii}.tiptilt = DM.segList{ii}.tiptilt + [tipList(mapped_segment),tiltList(mapped_segment)];
    end
end

%% Go Ahead and Render it
DM.touch;
DM.render;
end

