% Когда приходит текущий фазовый вектор phi, подсчитывается оценка выигрыша при движении вдоль него от каждой оценки, которая хранится в памяти. Старые оценки также сохраняются.
% Отбрасываются слишком старые оценки в соответствии с заданной глубиной
% Оставшиеся сортируются по оценке суммарного выигрыша. Аутсайдеры забываются с тем, чтобы количество оставленных фаворитов не превышало заданный предельный объем.
% Текущий чемпион используется для синтеза управления.
% Происходит переход к следующему шагу.

clearvars
close all
prompt = {'\fontsize{15} Объект y_{k+1}=ay_k+bu_k+v_{k+1}; введите значения a и b (через пробел)',...
'\fontsize{15} Введите параметры алгоритма: \rho, глубина и объем памяти',...
'\fontsize{15} Введите начальные начения оценок a_0 и b_0',...
'\fontsize{15} Введите имя файла с результатами'};
titl = 'Параметры эксперимента';
dims = [1 180];
opts.Interpreter = 'tex';
definput = {'1.1 1.75','0.9 0 30', '1 -1','results.txt'};
answer = inputdlg(prompt, titl, dims, definput, opts);
ab = str2num(answer{1});
a = ab(1);
b = ab(2);
rdmn = str2num(answer{2});
rho = rdmn(1);
depth = rdmn(2); % глубина памяти
mn = rdmn(3); % объем памяти
tau = str2num(answer{3});
tau = tau';
tau = flipud(tau);
fname = answer{4};
rng('default');
v = rand(1)*2-1;
u = 0;
y = v;
Tau = tau;
a0 = tau(2);
b0 = tau(1);
D = 0;
K = 0;
phi = [u; y];
l50 = 150;
for k = 1:l50 
    u = -tau(2) * y/tau(1);
    phi = [u; y];
    v = rand(1)*2-1;
    y = a * y + b * u + v;
    Y(k) = y;
    U(k) = u;
    [newtau,newD,newK] = stripe(y, Tau, phi, 1.1/rho, rho, k, D, K);
    l = length(newK);
    [Tau, D, K] = squ(newtau, newD, newK, depth, mn, l, k);
    tau = Tau(1:2, length(K));
    de = tau - [b; a];
    dt(k) = de' * de;
end
figure(2)
plot(dt)
figure(3)
plot([b b],[a a],'g*')
figure(1)
plot(Y)
Y = abs(Y);
maxy = max(Y);
thr = 1.1/rho;
plot([1:l50], Y(1:l50), 'b-', [1,l50], [thr thr], 'r-', [1,l50], -[thr thr], 'r-')
title(['y_{k+1}=ay_k+bu_k+v_{k+1}, a=' num2str(a) ', b=' num2str(b) ', a_0=' num2str(a0) ', b_0=' num2str(b0)])
fid0 = fopen(fname, 'rt', 'n', 'UTF-8');

if fid0 > 0
    fclose(fid0);
else
    
fid = fopen(fname,'at', 'n', 'UTF-8');
fprintf(fid,'%s\n', 'ПРОТОКОЛ ЭКСПЕРИМЕНТОВ');
fprintf(fid,'\n%s', 'Обозначения:');
fprintf(fid,'\n%s', 'a и b --- коэффициенты уравнения объекта');
fprintf(fid,'\n%s', 'rho --- параметр Полоски');
fprintf(fid,'\n%s', 'depth --- глубина памяти');
fprintf(fid,'\n%s', 'mn --- максимальный объем памяти');
fprintf(fid,'\n%s', 'a0 и b0 --- начальные значения оценок');
fprintf(fid,'\n%s', 'max|y| --- амплитуда выхода');
fprintf(fid,'\n%s', 'K --- момент последней коррекции');
fprintf(fid,'\n%s', 'r2 --- квадрат расстояния от финальной оценки до');
fprintf(fid,'\n%s', ' вектора истинных параметров');
fprintf(fid,'\n\n%s', 'Figure 1 --- график выхода объекта');
fprintf(fid,'\n%s', 'Figure 2 --- квадрат расстояния от оценки до истины');
fprintf(fid,'\n%s', 'Figure 2 --- эволюция оценок');

    fprintf(fid,...
        '\n\n%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s',...
        '   a','   b','   rho','depth','mn','a0','b0','max|y|','  K','r2');
    fclose(fid);
end
fid = fopen(fname,'at+', 'n', 'UTF-8');
fprintf(fid,...
    '\n%6.2f\t%6.2f\t%6.2f\t%4.0f\t%4.0f\t%6.2f\t%6.2f\t%6.2f\t%6.2f\t%6.2f',...
    a, b, rho, depth, mn, a0, b0, maxy, K, dt(150));
fclose(fid);



function [newtau, newD, newK] = stripe(y, tau, phi, C, rho, k, D, K)
nc = size(tau);
nc = nc(2);
newtau = zeros(2, 2 * nc);
eta = y - tau' * phi;
newtau(1:2, 1:nc) = tau;
newD(1:nc) = D;
newK(1:nc) = K;
G2 = phi' * phi;
j = nc;
for i = 1:nc
    if abs(eta(i)) > C
        j = j + 1;
        figure(3)
        hold on
        newK(j) = k;
        delta = (eta(i) - C * rho * sign(eta(i)));
        newD(j) = D(i) + delta^2/G2;
        delta = delta/G2;
        newtau(1:2,j) = tau(1:2,i) + delta * phi;
        t1 = [-5,5];
        t2 = (y - C - t1 * phi(1))/phi(2);
        t3 = (y + C - t1*phi(1))/phi(2);
        t4 = (y - C/rho - t1*phi(1))/phi(2);
        t5 = (y + C/rho - t1*phi(1))/phi(2);
        plot(t1, t2, 'b-', t1, t3, 'b-', t1, t4, 'r-', t1, t5, 'r-')
        plot([tau(1,i), newtau(1,j)], [tau(2,i), newtau(2,j)], 'g-')
        axis([-5 5 -5 5]);
        tauk = ['\tau_{' num2str(k) '}'];
        text(newtau(1,j), newtau(2,j), tauk);
    end
end
newtau = newtau(1:2,1:j);
newD = newD(1:j);
newK = newK(1:j);
end

function [Taun,Dn,Kn] = squ(Tau,D,K,d,m,l,k)
ll = l;
for i = 1:l
    if K(i) < k - d
        Tau(1:2, i:ll-1) = Tau(1:2, i+1:ll);
        D(i:ll-1) = D(i+1:ll);
        K(i:ll-1) = K(i+1:ll);
        ll = ll-1;
    end
    if ll == i
        break
    end
end
[D,I] = sort(D);
bi = max(1, ll-m);
ll = max(ll, bi);
Dn = D(I(bi:ll));
Kn = K(I(bi:ll));
Taun = Tau(1:2, I(bi:ll));
end