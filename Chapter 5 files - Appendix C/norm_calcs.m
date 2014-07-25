load('stability test/rand1/1000terms/iteration1.mat')
%load('iteration1.mat')

Matrix1Finalp_it1 =Matrix1Finalp;
Matrix2Finalp_it1 =Matrix2Finalp;
Matrix3Finalp_it1 =Matrix3Finalp;

Matrix1Finaln_it1 =Matrix1Finaln;
Matrix2Finaln_it1 =Matrix2Finaln;
Matrix3Finaln_it1 =Matrix3Finaln;

load('stability test/rand1/1000terms/iteration2.mat')
%load('iteration2.mat')

Matrix1Finalp_it2 =Matrix1Finalp;
Matrix2Finalp_it2 =Matrix2Finalp;
Matrix3Finalp_it2 =Matrix3Finalp;

Matrix1Finaln_it2 =Matrix1Finaln;
Matrix2Finaln_it2 =Matrix2Finaln;
Matrix3Finaln_it2 =Matrix3Finaln;

normp1 =norm(Matrix1Finalp_it2-Matrix1Finalp_it1,'fro')
normp2 =norm(Matrix2Finalp_it2-Matrix2Finalp_it1,'fro')
normp3 =norm(Matrix3Finalp_it2-Matrix3Finalp_it1,'fro')

normn1 =norm(Matrix1Finaln_it2-Matrix1Finaln_it1,'fro')
normn2 =norm(Matrix2Finaln_it2-Matrix2Finaln_it1,'fro')
normn3 =norm(Matrix3Finaln_it2-Matrix3Finaln_it1,'fro')

columnsp1 = abs(sum(Matrix1Finalp_it2-Matrix1Finalp_it1,1))./abs(sum(Matrix1Finalp_it1,1));
rowsp1 = abs(sum(Matrix1Finalp_it2-Matrix1Finalp_it1,2))./abs(sum(Matrix1Finalp_it1,2));

columnsp2 = abs(sum(Matrix2Finalp_it2-Matrix2Finalp_it1,1))./abs(sum(Matrix2Finalp_it1,1));
rowsp2 = abs(sum(Matrix2Finalp_it2-Matrix2Finalp_it1,2))./abs(sum(Matrix2Finalp_it1,2));

columnsp3 = abs(sum(Matrix3Finalp_it2-Matrix3Finalp_it1,1))./abs(sum(Matrix3Finalp_it1,1));
rowsp3 = abs(sum(Matrix3Finalp_it2-Matrix3Finalp_it1,2))./abs(sum(Matrix3Finalp_it1,2));

columnsn1 = abs(sum(Matrix1Finaln_it2-Matrix1Finaln_it1,1))./abs(sum(Matrix1Finaln_it1,1));
rowsn1 = abs(sum(Matrix1Finaln_it2-Matrix1Finaln_it1,2))./abs(sum(Matrix1Finaln_it1,2));

columnsn2 = abs(sum(Matrix2Finaln_it2-Matrix2Finaln_it1,1))./abs(sum(Matrix2Finaln_it1,1));
rowsn2 = abs(sum(Matrix2Finaln_it2-Matrix2Finaln_it1,2))./abs(sum(Matrix2Finaln_it1,2));

columnsn3 = abs(sum(Matrix3Finaln_it2-Matrix3Finaln_it1,1))./abs(sum(Matrix3Finaln_it1,1));
rowsn3 = abs(sum(Matrix3Finaln_it2-Matrix3Finaln_it1,2))./abs(sum(Matrix3Finaln_it1,2));