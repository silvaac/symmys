function G = G_Hat_a(X)

G = (X(:,1)-X(:,end)).*X(:,2).*X(:,2);

