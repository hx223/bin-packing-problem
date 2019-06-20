import numpy as np
from scipy import sparse
# Input: k 
# Output: A 3 permuation matrix P of [n] where n=3^k, A indice matrix I of size m = 3(n+1)/2 and n+1, a size vector s of m items.
def gene_permu(K):
	P_init = [[1,2,3],[3,1,2],[2,3,1]]
	for k in range(1,K+1):
		n_div_3 = pow(3,k-1)
		if k is 1:
			P = P_init
		else:
			P_1 = P[0]
			P_2 = P[1]
			P_3 = P[2]

			P_1_new = P_1 + [x + y for x, y in zip(P_1, [n_div_3] * n_div_3)] + [x + y for x, y in zip(P_1, [2 * n_div_3] * n_div_3)]
			P_2_new = [x + y for x, y in zip(P_2, [2 * n_div_3] * n_div_3)] + P_2 + [x + y for x, y in zip(P_2, [ n_div_3] * n_div_3)] 
			P_3_new = [x + y for x, y in zip(P_3, [n_div_3] * n_div_3)] + [x + y for x, y in zip(P_3, [2 * n_div_3] * n_div_3)] + P_3

			P_new = [P_1_new,P_2_new,P_3_new]
			P = P_new
	return P

# Input P,
# Output a matrix of tuples after adding 0 at the end
def gene_tuple(P):

	tuple_P = []
	for p in P:
		p = p + [0]
		tuple_P = tuple_P + [p]
	
	tuple_P_new=[]
	for p in tuple_P:
		tuple_P_new_row = []
		for j in range(int(len(p)/2)):
			tuple_P_new_row = tuple_P_new_row + [(p[2*j],p[2*j+1])]
		tuple_P_new = tuple_P_new + [tuple_P_new_row]

	return tuple_P_new

# Input tuple paired matrx
# Output a dictionary with size 
def gene_size(tuple_P_new):
	d_piece_size = {}
	l = 0
	n = len(tuple_P_new[0]) * 2 - 1
	for row in tuple_P_new:
		l = l + 1
		j = 0

		row_key = "pi_"+ str(l)
		d_piece_size[row_key] = {}

		for t in row:
			j = j + 1
			size_t = 1.0/3.0 - (l + (2 * j / (n + 1.0) ))/48.0
			d_piece_size[row_key][t] = size_t

	return d_piece_size

# Input: P
# Output: patten matrix of size m times n+1, each row is an item, each column is a pattern
def gene_instance_x(P):
	tuple_P = []
	for p in P:
		p = p + [len(p)+1]
		tuple_P = tuple_P + [p]
	
	n = len(tuple_P[0]) - 1
	m = 3 * int((n + 1) / 2)

	l = -1 
	row=[]
	column=[]
	for p in tuple_P:
		l = l + 1
		for j in range(0,len(p),2):
			row_id  =  int(m / 3 * l + j/2) 
			col_id_1 = p[j] - 1
			col_id_2 = p[j+1] - 1
			row = row + [row_id] + [row_id]
			column = column + [col_id_1] + [col_id_2]
	
	data = np.ones(2*m, dtype=np.int8)
	matrix_P=sparse.csr_matrix((data, (row,column)), shape=(m , (n + 1) ), dtype=np.int8)  
	return matrix_P

# test
P = gene_permu(3)
T = gene_tuple(P)
#S = gene_size(T)
I = gene_instance_x(P)
start = I.indptr[5]
end = I.indptr[6]
print(I.indices[start:end])
print(I.data[start:end])