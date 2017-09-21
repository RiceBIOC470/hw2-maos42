function [ DNAseq ] = protein2dnaOptimized( pseq )

A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];
%rand_seq=datasample(DNA,N);
%pseq=('GlnGluTyrValAspGlyValIleGluGlySerGluGluSerLeuAspArgGlyAsnEndIlePheProHisHisThrSerProLeuThrSerCysIleGluGlyPheValArgValMetLeuSerTrpThrGlyProSerLeuArgValSerValCysLeuLeuSerAlaLysAlaLysPheSerArgIleTyrPheEndSerTyrArgSerTyrLeuValArgThrCysThrAlaLeuTyrMetThrSerAlaEndArgSerValAspLeuArgProAlaMetAspGlnAspPheSerGluTyrGlyGlnLysCysIleAlaLeuThrAlaAlaIleLysArgLeuHisTyrSerGlyProPheAspAsnAlaIleGlyArgAsnIleValLeuLysTrpValGluLysValGlnGlyEndAlaLeuLeuGlyLysLeuLeuPheSerAlaSerLeuLeuEndProArgHisArgSerPheLeuAlaValAsnArg')
%the sequence above was for checking the function to run properly.

ds = datastore('codons.csv');
while(hasdata(ds))
import = read(ds);
 % do stuff to t.
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

newnames= {'Gly','Glu','Asp','Val','Ala','Arg','Ser','Lys','Asn','Met','Ile','Thr','Trp','End','Cys','Tyr','Leu','Phe','Ser','Arg','Gln','His','Leu','Pro'};
newtriplets= {'GGC','GAG','GAC','GTG','GCC','AGA','AGC','AAG','AAC','ATG','ATC','ACC','TGG','TGA','TGC','TAC','TTG','TTC','TCC','CGC','CAG','CAC','CTG','CCC'};
DNAseq=replace(pseq,newnames,newtriplets);



end

