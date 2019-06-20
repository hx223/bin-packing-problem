from scipy import sparse, misc
import numpy as np
from gurobipy import *

LOG = True

# Output is an instance of lil_matrix
def gene_instance_all(k):
	n = pow(3,k)
	m = 3 * int ((n+1) / 2)
	permu_num = misc.comb(m,3,exact=True)
	mtx = sparse.lil_matrix((permu_num,m),dtype=np.int8)
	row_id = 0
	for i in range(m):
		for j in range(i+1,m):
			for h in range(j+1,m):
				col_id = [i,j,h]
				mtx[row_id,col_id] = [1,1,1]
				row_id = row_id + 1
	return mtx

def lil2csr(M):
	A = M.transpose()
	A = A.tocsr()
	return A

def enumBinPacking(k,A):
	n = pow(3,k)
	m = 3 * int ((n+1) / 2)
	P = misc.comb(m,3,exact=True)

	model = Model("bpp")
	x={}
	for i in range(P):
		x[i] = model.addVar(obj=1, vtype="B", name="x[%d]"%i)
	model.update()


	for i in  range(m): 
		start = A.indptr[i] 
		end = A.indptr[i+1] 
		var = [x[j] for j in A.indices[start:end]] 
		coef = A.data[start:end] 
		model.addConstr(LinExpr(coef,var), ">=", 1, name="cnstr1[%d]"%i)

	if not LOG:
		model.Params.outputflag = 0
	model.optimize()
	[print(x[i].X) for i in x.keys()]

	return 0

k = 3
M = gene_instance_all(k)
A = lil2csr(M)
enumBinPacking(k,A)
