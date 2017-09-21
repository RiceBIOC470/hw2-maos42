function [ probability ] = plotProbabilityORF( L )

% Part4. Write  a final function called plotProbabilityORF.m which takes
% N_ORF as an argument and makes a plot of the probabily of having an
% ORF at least this long as a function of the dnasequence length. Decide how the
% code should determine the lengths of dna sequence to test and implement
% your decision. 


rand_seq=[];
A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];
N=L*10; %the length of the ORFs should be around 10% of the sequence length
seql= 0:L:N*2;
seqn=1;
newprob=0;

%de aquí puede empezar
for N=seql
    
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
probability=count/a;
end

newprob(seqn)=probability;
seqn=seqn+1;

end

plot(seql, newprob);
xlabel('Sequence Length');
ylabel('Probability');
legend(num2str(L));

end

