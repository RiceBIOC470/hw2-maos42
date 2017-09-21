function [ probability ] = probabilityORF( N, L )

% Part 3. Write another function called probabilityORF that utilizes the functions from 
% Parts 1 and 2. It should take two inputs - a sequence length (N) and an length  of an ORF (N_ORF) and
% returns the probability that that a sequence of length N contains an ORF
% of at least length N_ORF

rand_seq=[];
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];

%de aquí puede empezar

a=1000; %las veces que se repita
count=0;
istore=0;

for loop=1:a %esto es el index

rand_seq=datasample(DNA,N);    
start=strfind(rand_seq,'ATG');
ends=[strfind(rand_seq,'TAG') strfind(rand_seq,'TGA') strfind(rand_seq,'TAA')];

nstart=length(start); %columnas
nends=length(ends); %rows

allORFs=bsxfun(@minus, ends, reshape(start(:),1,1,[])); %need to reshape this
rallORFs=reshape(allORFs,[nends,nstart]);

%rallentre3=rallORFs./3;
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

dORFs=((possibleORFs+3)>L);
gORFs=length(possibleORFs(dORFs));
istore(loop)=gORFs; %this was just for me to check each loop value.
count=sum(istore>=1);

end
probability=count/a;
end