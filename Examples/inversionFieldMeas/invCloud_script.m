for c=1:4
    thisOne = darkCloud(c,:);
    refl = thisOne(1:7)';
    AP = defaultSMARTSinput('mlw','cosZ',cosd(thisOne(end)),'altit',5);
    S = SMARTS295Main(getSMARTShome,AP);
    [oS(c),stats(c),P(c)] = invertSnowCloudIntgRefl(S.spectralTbl.waveL,'nm',...
        [S.spectralTbl.HorzDirect S.spectralTbl.HorzDiffuse],refl,...
        {'radius','watereq'},'icecloud','cosZ',cosd(thisOne(end)),'R0',DarkLoam,...
        'waveu','nm','sensor','landsatoli','bands',[1:7]);
    disp(c)
end