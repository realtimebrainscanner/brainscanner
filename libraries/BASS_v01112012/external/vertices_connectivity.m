function VertConn = vertices_connectivity(FV,VERBOSE);
%VERTICES_CONNECTIVITY - Generate the connections between vertices
% function VertConn = vertices_connectivity(FV,VERBOSE);
% function VertConn = vertices_connectivity(FV);
% FV is the standard matlab structure for faces and vertices,
%  where FV.faces is m x 3, and FV.vertices is n x 3.
%
% VertConn returned is vector of cells, one-to-one with each row of FV.vertices.
% VertConn{i} returns a row vector of the vertex numbers (rows in FV.vertices) that
%  are connected by faces to the ith vertex row in FV.vertices.
%
% Thus if we want to 'swell' the region around a vertex, VertConn{i} gives us the 
%  vertex numbers of those vertices that are adjacent.
%
% See also FACES_CONNECTIVITY

%<autobegin> ---------------------- 27-Jun-2005 10:46:02 -----------------------
% ------ Automatically Generated Comments Block Using AUTO_COMMENTS_PRE7 -------
%
% CATEGORY: Visualization
%
% At Check-in: $Author: silvin $  $Revision: 192 $  $Date: 2006-10-15 13:34:01 -0700 (Sun, 15 Oct 2006) $
%
% This software is part of BrainStorm Toolbox Version 27-June-2005  
% 
% Principal Investigators and Developers:
% ** Richard M. Leahy, PhD, Signal & Image Processing Institute,
%    University of Southern California, Los Angeles, CA
% ** John C. Mosher, PhD, Biophysics Group,
%    Los Alamos National Laboratory, Los Alamos, NM
% ** Sylvain Baillet, PhD, Cognitive Neuroscience & Brain Imaging Laboratory,
%    CNRS, Hopital de la Salpetriere, Paris, France
% 
% See BrainStorm website at http://neuroimage.usc.edu for further information.
% 
% Copyright (c) 2005 BrainStorm by the University of Southern California
% This software distributed  under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPL
% license can be found at http://www.gnu.org/copyleft/gpl.html .
% 
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%<autoend> ------------------------ 27-Jun-2005 10:46:02 -----------------------

% ----------------------------- Script History ---------------------------------
% 1994 by John C. Mosher, Ph.D.
% 19-Nov-1999 Based on May 1998 scripts.
% 19-May-2004 JCM Comments Cleaning
% ----------------------------- Script History ---------------------------------

if(~exist('VERBOSE','var')),
   VERBOSE = 1; % default non-silent running of waitbars
end

nv = size(FV.vertices,1);
[nf,ns] = size(FV.faces); % number of faces, number of sides per face

VertConn = cell(nv,1); % one empty cell per vertex

if(VERBOSE)
   % disp(sprintf('Processing %.0f faces',nf))
   hwait = waitbar(0,sprintf('Processing the Vertex Connectivity for %.0f faces',nf));
   drawnow %flush the display
end

for iv = 1:nf, % use each face's information
   if(VERBOSE)
      if(~rem(iv,max(round(nf/10),1))), % ten updates
         % fprintf(' %.0f',iv);d
         waitbar(iv/nf,hwait);
         drawnow %flush the display         
      end
   end
   for i = 1:ns, %each vertex of the face
      for j = 0:(ns-2), % each additional vertex
         VertConn{FV.faces(iv,i)}(end+1) = FV.faces(iv,rem(i+j,ns)+1);
      end
   end
end

if(VERBOSE),
   close(hwait);
   hwait = waitbar(0,sprintf('Now sorting the vertex results'));
   drawnow %flush the display
end

for i = 1:nv,
   if(VERBOSE)
      if(~rem(i,max(round(nf/10),1))), % ten updates
         waitbar(i/nv,hwait);
         drawnow %flush the display         
      end
   end
   VertConn{i} = unique(VertConn{i});
end
if VERBOSE 
    close(hwait);
end

return

