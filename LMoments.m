function LM=LMoments(X,bias)
% Biased and unbiased calculation of L-Moments up to 4th order
% Written by George Pouliasis
beta0=mean(X);
X=sort(X,'descend');
n=length(X);
switch bias
    case 'biased'
    %{
    beta1=1/n*sum((1-(transpose(1:n)-0.35)/n).*X(transpose(1:n)));
    beta2=1/n*sum((1-(transpose(1:n)-0.35)/n).^2.*X(transpose(1:n)));
    beta3=1/n*sum((1-(transpose(1:n)-0.35)/n).^3.*X(transpose(1:n)));
    %}
    beta1=1/n*sum((1-(n-transpose(1:n)-0.65)./n).*X(n-transpose(1:n)+1));
    beta2=1/n*sum(((1-(n-transpose(1:n)-0.65)./n).^2).*X(n-transpose(1:n)+1));
    beta3=1/n*sum(((1-(n-transpose(1:n)-0.65)./n).^3).*X(n-transpose(1:n)+1));
    LM(1)=beta0;
    LM(2)=2*beta1-beta0;
    LM(3)=6*beta2-6*beta1+beta0;
    LM(4)=20*beta3-30*beta2+12*beta1-beta0;
    case 'unbiased'
    ubeta1=1/n*sum(((transpose(2:n)-1)/(n-1)).*X(n-transpose(2:n)+1));
    ubeta2=1/n*sum(((transpose(3:n)-1).*(transpose(3:n)-2)/(n-1)/(n-2)).*X(n-transpose(3:n)+1));
    ubeta3=1/n*sum(((transpose(3:n)-1).*(transpose(3:n)-2).*(transpose(3:n)-3)/(n-1)/(n-2)/(n-3)).*X(n-transpose(3:n)+1));
    LM(1)=beta0;
    LM(2)=2*ubeta1-beta0;
    LM(3)=6*ubeta2-6*ubeta1+beta0;
    LM(4)=20*ubeta3-30*ubeta2+12*ubeta1-beta0;
end