function projstruct = HimachalPradeshUTM()
% projstruct = HimachalPradeshUTM()
%return the projection structure for the UTM projection over the
%Himachal Pradesh region of India, where we have AVIRIS-NG data

projstruct = createUTMmstruct([32.753 77.26],'wgs84');

end