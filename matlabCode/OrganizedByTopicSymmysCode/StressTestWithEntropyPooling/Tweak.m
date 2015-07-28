function db=Tweak(A,b,g)

[K,J]=size(A);

g_=[g; zeros(J,1)];

Aeq_=[zeros(1,K) ones(1,J)];
beq_=1;
lb_=[zeros(K,1); zeros(J,1)];
ub_=[Inf*ones(K,1); ones(J,1)];
A_=[-eye(K) A];
b_=b;
db_ = linprog(g_,A_,b_,Aeq_,beq_,lb_,ub_);
db = db_(1:K);

