function result = teb(N, EbNo)

%N = 10000; # N = 32 pour visualisations temporelles
D = 10000000;
T = 1/D;
EbNodB = 20;

# 1)
bn = zeros(1, N);
for k=1:1:length(bn)
    bn(k) = round(rand());
end

# 2)
an = zeros(1, N);
for k=1:1:length(bn)
    an(k) = 2*bn(k)-1;
end

# 3-1)
F = 16; # Facteur de surechantillonage
st = zeros(1, N*F);
st(1) = F*an(1);
for k=1:1:length(an)-1
    st(k*F+1) = F*an(k+1);
end

#3-2)
L = 8;
alpha = 0.5;
Te = T/F;
t_filtre = [0 : T/F : L*T -T/F];

s_t = gen_filters2('srrc',t_filtre,T,F,L,alpha);

#3-3)
#QUESTION 3
ht = gen_filters2('srrc',t_filtre,T,F,L,0.5);
xt = conv(st, ht);

#Question 7
# P(xt) = T (car variance unitaire et filtre normalisé)
# Eb/No = (F/2)*(T/sigma_n^2) donc sigma_n^2 = (F*T/2)/(EbNo)
sigma_n = sqrt((F/2)/EbNo);
nt=sigma_n*randn(1,length(xt));
rt = xt + nt;

#Question 9
htr = fliplr(ht+L*T);

pt = conv(ht, htr);

t_filtre_pt = [0 : T/F : 2*(L*T + T/F)];

vect_result = zeros(1,17);

for k=1:1:8
vect_result(9-k) = pt(round((L*T-k*T+2*T/F)/(T/F)));
vect_result(9+k) = pt(round((L*T+k*T+2*T/F)/(T/F)));
end

vect_result(9) = pt(round((L*T+2*T/F)/(T/F)));

yt = conv(rt, htr);

yk = [];

for k=1:1:length(yt)-1
    if mod(k, F) == 0 && -L*F*T + k*T >= 0 && -L*F*T + k*T < (N*F)*T
        yk = [yk F*yt(k)];
    end
end

bnF = [];

for k=1:1:length(yk)
    val = yk(k)/16;
    if val > 1 - sigma_n
        bnF = [bnF 1];
    else
        bnF = [bnF 0];
    end
end

result = sum(xor(bn, bnF))/N;