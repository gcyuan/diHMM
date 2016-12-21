function domainAnnotation = getDomainAnnotations
%   GETDOMAINANNOTATION generates standard nucleosome-level annotations
%   and colors
%
%   Author: Eugenio Marco

domainAnnotation(1).name = 'Broad Promoters';
domainAnnotation(1).color = [191, 0, 0];

domainAnnotation(2).name = 'Promoters/Exons';
domainAnnotation(2).color = [255, 153, 153];

domainAnnotation(3).name = 'Bivalent Promoters';
domainAnnotation(3).color = [238, 130, 238];

domainAnnotation(4).name = 'Poised Enhancers';
domainAnnotation(4).color = [215, 137, 105];

domainAnnotation(5).name = 'Super-Enhancers';
domainAnnotation(5).color = [253, 199, 72];

domainAnnotation(6).name = 'Upstream Enhancers';
domainAnnotation(6).color = [242, 232, 76];

domainAnnotation(7).name = 'Introns/Enhancers';
domainAnnotation(7).color = [242, 232, 145];

domainAnnotation(8).name = 'Transcribed';
domainAnnotation(8).color = [0, 181, 94];

domainAnnotation(9).name = 'Boundary';
domainAnnotation(9).color = [159, 249, 255];

domainAnnotation(10).name = 'Polycomb Repressed';
domainAnnotation(10).color = [133, 125, 113];

domainAnnotation(11).name = 'Heterochrom; Low Signal';
domainAnnotation(11).color = [255, 255, 255];

domainAnnotation(12).name = 'Satellite';
domainAnnotation(12).color = [238, 207, 161]; % navajowhite2

domainAnnotation(13).name = 'Low Coverage';
domainAnnotation(13).color = [130, 255, 138]; % navajowhite2

end