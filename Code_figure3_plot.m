% Code to plot the prediction results already computed in VoteEnsemble10.mat


close all;clear all;clc;
load VoteEnsemble10.mat;

% number of cell lines
n=10;
% number of compounds
o=614;
% number of targets
dc=114;

A=1:n;

% column with the gene ids
c=198;

% read the column 'TargetGeneNum' with a unique identifier for each target gene
M01=readtable('MPM04_R1_v1.csv');
label=M01{:,c};

% loop on the number of cell lines
for i=1:n
    i

% all combinations of i cell lines
B=combnk(A,i);

% loop on all combinations of i cell lines
for j=1:size(B,1)
    N=B(j,:);

    %%%%%% hard voting %%%%%%
    %S=K(:,N);
    %Kvote=mode(S,2);
    
    %%%%%% soft voting %%%%%%

    % ClassProb comes from VoteEnsemble10.mat
    temp=ClassProb(:,N,:);
    S=sum(temp,2);
    S1= reshape(S,o,dc);
    [m,K5vote]=sort(S1,2,'descend');
    top5=K5vote(:,1:5);

    % loop of the number of compounds (leave-one-out cross-validation)
    for k=1:o
        temp1(k)=any(top5(k,1:5)==label(k));
    end

    T5(j)=sum(temp1)/o;
    [m,Kvote]=max(S,[],3);
    [G,order] = confusionmat(label,Kvote);
    T1(j)=trace(G)/sum(G(:));

end % end loop combinations

    % keep values to plot
    C1(i,1)=mean(T1);
    C1(i,2)=std(T1);
    C5(i,1)=mean(T5);
    C5(i,2)=std(T5);

    %reinitialize
    T1=0;
    T5=0;

end % end loop on the cell lines

% plot and save
errorbar(1:n,C5(:,1)*100,C5(:,2)*100, 'linewidth',2, 'Color', [0.3 0.3 0.8]);hold on;
errorbar(1:n,C1(:,1)*100,C1(:,2)*100, 'linewidth',2, 'Color', [0.5 0.7 0.5]);hold on;
a=ones(1,n)*0.78;
l=ones(1,n)*0.28;
u=ones(1,n)*4.42;
plot(1:n,a,'r','LineWidth',2);hold on;
%plot(1:n,u,'r--','LineWidth',1);hold on;
%plot(1:n,l,'r--','LineWidth',1);hold on;
xlim([0.5 11.5]);
xlabel('Number of cell lines','FontSize', 16,'FontName','Helvetica');
ylabel('Percentage accuracy','FontSize', 16,'FontName','Helvetica');
%title('Soft voting on ensemble classifier prediction','FontSize', 20,'FontName','Helvetica');
%legend('Ensemble voting top 5','Ensemble voting top 1','Random chance', 'Location', 'southwestoutside');
print('emsembleVoting10.pdf','-r300','-dpdf');
saveas(gcf, 'ensembleVoting10.svg');
