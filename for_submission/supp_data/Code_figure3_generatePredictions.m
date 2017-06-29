close all;clear all;clc;

% start and end of lines (=compounds) in the data files you want to consider 
s=1;e=621;
% number of features for the MPM (mesothelioma) files
f=196;
%number of features for the DoceLessN (Prostate cancer) file
f1=38;



% load MPM files replicate 1
M01=readtable('MPM04_R1_V1.csv');
M03=readtable('MPM10_R1_V1.csv');
M04=readtable('MPM11_R1_V1.csv');
M05=readtable('MPM17_R1_V1.csv');
M06=readtable('MPM24_R1_V1.csv');
M07=readtable('MPM25_R1_V1.csv');
M08=readtable('MPM28_R1_V1.csv');
M09=readtable('MPM59_R1_V1.csv');
M10=readtable('MPM60_R1_V1.csv');

% load MPM files replicate 2
M012=readtable('MPM04_R2_V1.csv');
M032=readtable('MPM10_R2_V1.csv');
M042=readtable('MPM11_R2_V1.csv');
M052=readtable('MPM17_R2_V1.csv');
M062=readtable('MPM24_R2_V1.csv');
M072=readtable('MPM25_R2_V1.csv');
M082=readtable('MPM28_R2_V1.csv');
M092=readtable('MPM59_R2_V1.csv');
M102=readtable('MPM60_R2_V1.csv');

% load DoceLessN files, replicates 1&2
M11=readtable('Prestwick_R100_DoceLessN_R1_V1.csv');
M12=readtable('Prestwick_R100_DoceLessN_R2_V1.csv');



t=M012{s:e,1:f};



% average the two replicates and z-score normalize the features
F01=zscore((M01{s:e,1:f}+M012{s:e,1:f})/2);
F03=zscore((M03{s:e,1:f}+M032{s:e,1:f})/2);
F04=zscore((M04{s:e,1:f}+M042{s:e,1:f})/2);
F05=zscore((M05{s:e,1:f}+M052{s:e,1:f})/2);
F06=zscore((M06{s:e,1:f}+M062{s:e,1:f})/2);
F07=zscore((M07{s:e,1:f}+M072{s:e,1:f})/2);
F08=zscore((M08{s:e,1:f}+M082{s:e,1:f})/2);
F09=zscore((M09{s:e,1:f}+M092{s:e,1:f})/2);
F10=zscore((M10{s:e,1:f}+M102{s:e,1:f})/2);

F11=zscore((M11{s:e,9:46}+M11{s:e,9:46})/2);



% all the MPM data
FF=[F01,F03,F04,F05,F06,F07,F08,F09,F10];
% DoceLessN data
FFF=F11;



tic

% loop on the compounds
for test_ind=s:e

test_ind
% all the rows without the tested compound (test_ind)
train_ind=~ismember(1:e, test_ind);

% column 'TargetGeneNum' with a unique identifier for each target, with the rows corresponding to all the trained target (train_ind)
C01=M01{train_ind,201};

% parallelized for loop on the 9 MPM cell lines
parfor i=1:9   
  % train data
  Train=FF(train_ind,f*(i-1)+1:f*i);
  % test data
  Test=FF(test_ind,f*(i-1)+1:f*i);

  % train the model
  Mdl = fitensemble(Train,C01,'Bag',30,'Tree','type','classification');
  % get prediction and store it
  [K(test_ind,i),ClassProb(test_ind,i,:)]= predict(Mdl,Test);

end % end loop MPM cell lines

% same with the DoceLessN cell line
i=1;
  Train=FFF(train_ind,f1*(i-1)+1:f1*i);
  Test=FFF(test_ind,f1*(i-1)+1:f1*i);
  Mdl = fitensemble(Train,C01,'Bag',30,'Tree','type','classification');
  [K(test_ind,10+i),ClassProb(test_ind,10+i,:)]= predict(Mdl,Test);

end % end loop compounds



toc



% save ClassProb
save VoteEnsemble10.mat ClassProb;


