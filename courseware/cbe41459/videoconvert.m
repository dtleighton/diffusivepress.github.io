% This script takes in an mp4 movie which has been saved from imovie (or
% other video program - including directly from the camera memory stick)
% and converts it to a sequence of jpegs for analysis.  It will take up to
% the first 99 frames of the clip.  The "basename" is the name of the movie
% file, which is the name it was "shared" under out of imovie.  The script
% will put the jpegs into a folder with that basename, and the basename is
% also the first part of each file name. This structure is compatible with
% the existing version of pivanal.
%
% The code is adapted from one written by chatGTP.

basename=input('Enter the base filename: ','s');
inputVideoFileName = [basename,'.mp4'];

% Create a VideoReader object to read the input video
videoReader = VideoReader(inputVideoFileName);

% We set the output folder to be the base filename
outputFolder = [basename,'/'];

% Create the output folder if it doesn't exist
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% We append the output folder to the matlab path
path(path,outputFolder)

% Loop through each frame in the video
frameCount = 1;
while hasFrame(videoReader)&&frameCount<100
    % Read the current frame
    currentFrame = readFrame(videoReader);
    
    % Define the output image file name
    outputImageFileName = sprintf(['%s',basename,' %02d.jpg'], outputFolder, frameCount);
   
    % Save the current frame as a JPEG image
    imwrite(currentFrame, outputImageFileName);
    
    % Move to the next frame
    frameCount = frameCount + 1;
end

fprintf('Conversion complete! %d frames saved as JPEG images.\n', frameCount - 1);
