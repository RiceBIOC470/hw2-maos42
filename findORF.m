function [MAXORF3, mstart, mend] = findORF(seq)
%Function to find the length of the longest open reading frame of a
%sequences called dnaseq

%fill in here.
rand_seq=upper(seq); %turns all to upper cases to avoid any problems

start=strfind(rand_seq,'ATG');
ends=[strfind(rand_seq,'TAG') strfind(rand_seq,'TGA') strfind(rand_seq,'TAA')];

nstart=length(start); %columnas
nends=length(ends); %rows

allORFs=bsxfun(@minus, ends, reshape(start(:),1,1,[])); %need to reshape this
rallORFs=reshape(allORFs,[nends,nstart]);

rallentre3=rallORFs./3;
rall=rallORFs./3;
readsall=(~mod(rall,1));
rall(rall<1 | ~readsall)=0;
deltaall=rall.*3; %now the difference is the real bp distance
ndeltaall=deltaall;
ndeltaall(ndeltaall<1)=NaN;
possibleORFs=min(ndeltaall(:,1:nstart));
possibleORFs(isnan(possibleORFs))=[];
nORFs=length(possibleORFs);

MAXORF=max(possibleORFs); %falta el +3, lo puse al final del resultado
MAXORF3=MAXORF+3;

if isempty(MAXORF)
    disp('No ORF found');
else
fprintf('\n The sequence for the longest ORF is:   %d   bp', MAXORF+3)

posstart=(MAXORF==possibleORFs);
posend=(MAXORF==(deltaall(:,posstart)));
mstart=start(posstart);
mend=ends(posend);
fprintf('\n starts at:  %d   bp', mstart);
fprintf('\n ends at:  %d   bp', mend);
end

end