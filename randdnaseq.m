function rand_seq = randdnaseq(N)
% returns a random dna sequence of length N

A='A'; T='T'; G='G'; C='C';
DNA=[A,T,G,C];
rand_seq=datasample(DNA,N);
end

