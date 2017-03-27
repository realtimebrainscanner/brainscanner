function QG=generate_basisfunctions(vert,face,s)
Nd=size(vert,1);

A     = spm_mesh_distmtx(struct('vertices',vert,'faces',face),0);
GL    = A - spdiags(sum(A,2),0,Nd,Nd);
GL    = GL*s/2;
Qi    = speye(Nd,Nd);
QG    = sparse(Nd,Nd);
for i = 1:8
    QG = QG + Qi;
    Qi = Qi*GL/i;
end
QG    = QG.*(QG > exp(-8));
QG    = QG*QG;
clear Qi A GL
fprintf(' - done\n')