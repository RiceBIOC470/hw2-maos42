function protein = dna2protein(seq)
% Input a dna sequence and it will return all the reading frames with the
% longest protein available for each.

rand_seq=upper(seq);

ds = datastore('codons.csv');
while(hasdata(ds))
import = read(ds);
 
end

cod=table2cell(import);
names=cod(:,1);
cell2mat(names);
triplets=cod(:,2);
cell2mat(triplets);
frequency=cod(:,3);
cell2mat(frequency);

codons(:,1)=names;
codons(:,2)=triplets;
codons(:,3)=frequency;

%for determining the ORF which will be read
N=length(rand_seq);
isd=N/3;
isit3=(~mod(isd,1));
iseven=N/2;
isit2=(~mod(iseven,1));

if isit3==1 && isit2==1
    RF1=rand_seq(1:N);
    RF2=rand_seq(2:N-2);
    RF3=rand_seq(3:N-1);
end
if  isit2==1 && isit3==0
    RF1=rand_seq(1:N-2);
    RF2=rand_seq(2:N-1);
    RF3=rand_seq(3:N);
end
if  isit2==0 && isit3==0
    RF1=rand_seq(1:N-1);
    RF2=rand_seq(2:N);
    RF3=rand_seq(3:N-2);
end
    
pRF1=replace(RF1,triplets, names );
pRF2=replace(RF2,triplets, names );
pRF3=replace(RF3,triplets, names );

pept1=1;
pept2=1;
pept3=1;

%for ORF1
start1=strfind(pRF1,'Met');
nstart1=length(start1);

emptys1=isempty(start1);
if emptys1==1
    start1=0;
end
end1=strfind(pRF1,'End');
nend1=length(end1);

emptye1=isempty(end1);
if emptye1==1
    end1=0;
end
if start1==0 
    pept1=0;
    maxpept1=0;

end
if end1==0
   pept1=0;
   maxpept1=0;

end

if pept1~=0
r1=bsxfun(@minus, end1, reshape(start1(:),1,1,[]));
rr1=reshape(r1,[nend1,nstart1]);
rr1(rr1<1)=NaN;;
pept1=min(rr1(:,1:nstart1));
maxpept1=max(pept1);
end

empty1=isempty(pept1);
if empty1==1;
    pept1=0;
    maxpept1=0;
end
if pept1<0;
    pept1=0;
    maxpept1=0;
end

% for ORF2
start2=min(strfind(pRF2,'Met'));
nstart2=length(start2);
emptys2=isempty(start2);
if emptys2==1
    start2=0;
end
end2=strfind(pRF2,'End');
nend2=length(end2);
emptye2=isempty(end2);
if emptye2==1
    end2=0;
end
if start2==0 
    pept2=0;
    maxpept2=0;
end
if end2==0
   pept2=0;
   maxpept2=0;
end

if pept2~=0
r2=bsxfun(@minus, end2, reshape(start2(:),1,1,[]));
rr2=reshape(r2,[nend2,nstart2]);
rr2(rr2<1)=NaN;;
pept2=min(rr2(:,1:nstart2));
maxpept2=max(pept2);
end

empty2=isempty(pept2);
if empty2==1;
    pept2=0;
    maxpept2=0;
end
if pept2<0;
    pept2=0;
    maxpept2=0;
end

%for ORF3
start3=min(strfind(pRF3,'Met'));
nstart3=length(start3);
emptys3=isempty(start3);
if emptys3==1
    start3=0;
end
end3=strfind(pRF3,'End');
nend3=length(end3);
emptye3=isempty(end3);
if emptye3==1
    end3=0;
end
if start3==0 
    pept3=0;
    maxpept3=0;

end
if end3==0
   pept3=0;
   maxpept3=0;

end

if pept3~=0
r3=bsxfun(@minus, end3, reshape(start3(:),1,1,[]));
rr3=reshape(r3,[nend3,nstart3]);
rr3(rr3<1)=NaN;;
pept3=min(rr3(:,1:nstart3));
maxpept3=max(pept3);
end

empty3=isempty(pept3);
if empty3==1;
    pept3=0;
end
if pept3<0;
    pept3=0;
end
peptide=max(max(maxpept1,maxpept2),maxpept3);

if maxpept1>maxpept2 && maxpept1>maxpept3
posstart=(peptide==pept1);
posend=(peptide==(rr1(:,posstart)));
mstart=start1(posstart);
mend=end1(posend);
mend=mend+2;
protein=pRF1(mstart:mend);
end

if maxpept2>maxpept1 && maxpept2>maxpept3
posstart=(peptide==pept2);
posend=(peptide==(rr2(:,posstart)));
mstart=start2(posstart);
mend=end2(posend);
mend=mend+2;
protein=pRF2(mstart:mend);
end

if maxpept3>maxpept2 && maxpept3>maxpept1
posstart=(peptide==pept3);
posend=(peptide==(rr3(:,posstart)));
mstart=start3(posstart);
mend=end3(posend);
mend=mend+2;
protein=pRF3(mstart:mend);
end

lprot=length(protein)/3;
lmp1=(maxpept1+3)/3;
if lmp1<3
   lmp1=0;
end

lmp2=(maxpept2+3)/3;
if lmp2<3
   lmp2=0;
end

lmp3=(maxpept3+3)/3;
if lmp3<3
   lmp3=0;
end

fprintf('\n in ORF1 your longest protein would be:  %d aa long', lmp1);
fprintf('\n in ORF2 your longest protein would be:  %d aa long', lmp2);
fprintf('\n in ORF3 your longest protein would be:  %d aa long', lmp3);

fprintf('\n \n And the longest possible peptide from you sequence is:  %d   aa long', lprot);

end