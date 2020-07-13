[X0,genelist0]=sc_readmtxfile("0Mild.mtx","allGeneList.txt");
[X1,genelist1]=sc_readmtxfile("0Severe.mtx","allGeneList.txt");
[X0,g0]=sc_selectg(X0,genelist0,1,1);
[X1,g1]=sc_selectg(X1,genelist1,1,1);
[genelist,i,j]=intersect(g0,g1);
X0=X0(i,:);
X1=X1(j,:);
X0=sc_selectc(X0,100);
X1=sc_selectc(X1,100);
X=[X0 X1];

%%
rng default
s=sc_view2grpcells(X0,X1);
legend({'Mild','Severe'},"Location","best")
xlabel('t-SNE 1'); ylabel('t-SNE 2'); zlabel('t-SNE 3');

figure;
sc_view2grpcells(X0,X1,s);
legend({'Moderate (n=1125)','Severe (n=3735)'},"Location","best")
xlabel('t-SNE 1'); ylabel('t-SNE 2'); zlabel('t-SNE 3');
view(0,90);


tgene="CD68";
tgene="TLR2";

figure;
sc_stemscatter(s(:,1),s(:,2),log(1+X(genelist==tgene,:)));    
title(tgene)

%%
T=readtable('result_DR_genelist_mild_severe.txt');
FC=T.FC;
genelistx=string(T.genelist);
        pd = makedist('Gamma','a',0.5,'b',2);
        figure;
        qqplot(FC,pd);
        [~,i]=sort(FC);
        dt = datacursormode;
        dt.UpdateFcn = {@i_myupdatefcn1,genelistx(i)};
title('')
xlabel('Expected Manifold Distance FC')
ylabel('Observed Manifold Distance FC')


