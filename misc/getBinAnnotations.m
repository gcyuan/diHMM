function binAnnotation = getBinAnnotations
%   GETBINANNOTATIONs generates standard nucleosome-level annotations and colors
%
%   Author: Eugenio Marco


binAnnotation(1).name = 'Active Promoter';
binAnnotation(1).color = [255, 0, 0];

binAnnotation(2).name = 'Promoter Flanking';
binAnnotation(2).color = [255, 182, 193]; % Light pink

binAnnotation(3).name = 'Bivalent Promoter';
binAnnotation(3).color = [238, 130, 238];

binAnnotation(4).name = 'Poised Enhancer';
binAnnotation(4).color = [215, 137, 105];

binAnnotation(5).name = 'Strong Enhancer';
binAnnotation(5).color = [253, 199, 72];

binAnnotation(6).name = 'Weak Enhancer';
binAnnotation(6).color = [242, 232, 76];

binAnnotation(7).name = 'Transcribed Enhancer';
binAnnotation(7).color = [191, 255, 0]; % Lime green

binAnnotation(8).name = 'Transcriptional Elongation';
binAnnotation(8).color = [0, 181, 94];

binAnnotation(9).name = 'CTCF-Promoter';
binAnnotation(9).color = [123, 104, 238]; % mediumslateblue

binAnnotation(10).name = 'CTCF';
binAnnotation(10).color = [0, 187, 221];

binAnnotation(11).name = 'rRNA Rich';
binAnnotation(11).color = [0, 140, 140];

binAnnotation(12).name = 'Polycomb Repressed';
binAnnotation(12).color = [133, 125, 113];

binAnnotation(13).name = 'Heterochrom; Low Signal';
binAnnotation(13).color = [255, 255, 255];

binAnnotation(14).name = 'Repetitive/CNV';
binAnnotation(14).color = [238, 207, 161]; % navajowhite2

end