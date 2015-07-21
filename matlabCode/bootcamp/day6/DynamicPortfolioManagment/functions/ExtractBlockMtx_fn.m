function Block_t1_t2 = ExtractBlockMtx_fn(A,t1,t2,n1_,n2_)
%This function extracts the (t1,t2)-block out of the block-diagonal matrix A
%A is a matrix t_*n1_ x t_*n2_. The matrix has t_ blocks. Each block is
%n1_ x n2_

Block_t1_t2 = A(1+(t1-1)*n1_:t1*n1_,1+(t2-1)*n2_:t2*n2_);