% Created for saving video of model evolution
% Author: Yu-Han Wang
% Date: 25 01 2025
fig = figure('Visible', 'off');
ouputFolderVideo = fullfile(parentFolder, 'post', outputFolder, 'localization_agu.avi');
[X,Y]  = meshgrid(xx,yy);
outputFolder 

% Delete outdated files
if isfile(ouputFolderVideo)
   delete(ouputFolderVideo);
   disp('video 2 deleted.')
end

% Setup VideoWriter for Figure 2
videoFile2 = VideoWriter(ouputFolderVideo, 'Motion JPEG AVI');
videoFile2.FrameRate = 10; 
open(videoFile2);

