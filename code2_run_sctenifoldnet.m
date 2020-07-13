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


rng default
T=sctenifoldnet_m(X0,X1,genelist,'nsubsmpl',20);
writetable(T,'res1.txt');

Tf=run_fgsea(T);
writetable(Tf,'res2.txt');

