function Seq = HaltonSeq(N, Base)
Seq = zeros(N,1);
NumBits = 1+ceil(log(N)/log(Base));
VetBase = Base.^(-(1:NumBits));
WorkVet = zeros(1,NumBits);
for i=1:N
j=1;
ok = 0;
while ok == 0
WorkVet(j) = WorkVet(j)+1;
if WorkVet(j) < Base
ok = 1;
else
WorkVet(j) = 0;
j = j+1;
end
end
Seq(i) = dot(WorkVet,VetBase);
end