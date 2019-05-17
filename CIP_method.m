clear
clc

%% äiéqì_ÇÃê›íË
X = 101;
Y = 101;
NT = 301;

f_p    = zeros(X,Y);
f_m    = zeros(X,Y);
g_p    = zeros(X,Y);
g_m    = zeros(X,Y);
fn_p   = zeros(X,Y);
fn_m   = zeros(X,Y);
gn_p   = zeros(X,Y);
gn_m   = zeros(X,Y);

P      = zeros(X,Y);
dx_P   = zeros(X,Y);
dy_P   = zeros(X,Y);

ZUx     = zeros(X,Y);
dx_ZUx  = zeros(X,Y);

ZUy     = zeros(X,Y);
dx_ZUy  = zeros(X,Y);

coeff10 = zeros(8);
coeff31 = zeros(32);

%% èâä˙ílÇÃê›íË

xc= (X-1) / 2;
yc= (Y-1) / 2;

dx = 5.e-2;
dy = 5.e-2;
dt = 5.e-5;

Ro = 1.21;
bm = 1.4235529e5;
c0 = sqrt(bm / Ro);
Z0 = sqrt(bm * Ro);
sigma = 0.2;
%% äiéqì_Ç≤Ç∆ÇÃï‚ä‘åvéZ
Ua = c0;
xi =-Ua * dt;
C  = c0 * dt / dx;
C2 = C  * C;
C3 = C2 * C;

coeff10(1)  = C * 0.5;
coeff10(2)  = (-C + 1.) * 0.5;

coeff31(1)  = (-2. * C3 + 3. * C2) * 0.5;
coeff31(2)  = (2. * C3 - 3. * C2 + 1.) * 0.5;
coeff31(3)  = xi * (C2 - C) * 0.5;
coeff31(4)  = xi * (C2 - 2. * C + 1.) * 0.5;
coeff31(5)  = 6. * (-C3 + C2) / xi * 0.5;
coeff31(6)  = 6. * (C3 - C2) / xi * 0.5;
coeff31(7)  = (3. * C2 - 2. * C) * 0.5;
coeff31(8)  = (3. * C2 - 4. * C + 1.) * 0.5;


Ua = -c0;
xi =-Ua * dt;
C  = c0 * dt / dx;
C2 = C  * C;
C3 = C2 * C;

coeff10(3)  = C * 0.5;
coeff10(4)  = (-C + 1.) * 0.5;

coeff31(9)  = (-2. * C3 + 3. * C2) * 0.5;
coeff31(10)  = (2. * C3 - 3. * C2 + 1.) * 0.5;
coeff31(11)  = xi * (C2 - C) * 0.5;
coeff31(12)  = xi * (C2 - 2. * C + 1.) * 0.5;
coeff31(13)  = 6. * (-C3 + C2) / xi * 0.5;
coeff31(14)  = 6. * (C3 - C2) / xi * 0.5;
coeff31(15)  = (3. * C2 - 2. * C) * 0.5;
coeff31(16)  = (3. * C2 - 4. * C + 1.) * 0.5;


Ua = c0;
xi =-Ua * dt;
C  = c0 * dt / dx;
C2 = C  * C;
C3 = C2 * C;

coeff10(5)  = C * 0.5;
coeff10(6)  = (-C + 1.) * 0.5;

coeff31(17)  = (-2. * C3 + 3. * C2) * 0.5;
coeff31(18)  = (2. * C3 - 3. * C2 + 1.) * 0.5;
coeff31(19)  = xi * (C2 - C) * 0.5;
coeff31(20)  = xi * (C2 - 2. * C + 1.) * 0.5;
coeff31(21)  = 6. * (-C3 + C2) / xi * 0.5;
coeff31(22)  = 6. * (C3 - C2) / xi * 0.5;
coeff31(23)  = (3. * C2 - 2. * C) * 0.5;
coeff31(24)  = (3. * C2 - 4. * C + 1.) * 0.5;


Ua = -c0;
xi =-Ua * dt;
C  = c0 * dt / dx;
C2 = C  * C;
C3 = C2 * C;

coeff10(7)  = C * 0.5;
coeff10(8)  = (-C + 1.) * 0.5;

coeff31(25)  = (-2. * C3 + 3. * C2) * 0.5;
coeff31(26)  = (2. * C3 - 3. * C2 + 1.) * 0.5;
coeff31(27)  = xi * (C2 - C) * 0.5;
coeff31(28)  = xi * (C2 - 2. * C + 1.) * 0.5;
coeff31(29)  = 6. * (-C3 + C2) / xi * 0.5;
coeff31(30)  = 6. * (C3 - C2) / xi * 0.5;
coeff31(31)  = (3. * C2 - 2. * C) * 0.5;
coeff31(32)  = (3. * C2 - 4. * C + 1.) * 0.5;

%%
for i = range(1, X-1)
    x = dx * i;
    for j = range(1, Y-1)
        y = dy * j;
        TX = x - xc * dx;
        TY = y - yc * dy;
        P(i,j)     = ...
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma^2));
        dx_P(i,j)  = -TX * ...
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma^2)) / sigma**2;
        dy_P(i,j)  = -TY * ...
                    exp(((-TX * TX) + (-TY * TY)) / (2. * sigma^2)) / sigma**2;
    end
end

%%
function LINEAR(coeff0, coeff1, f0, f1)
    return    coeff0 * f0 + coeff1 * f1;
end

function CIP(coeff0, coeff1, coeff2, coeff3, f0, f1, g0, g1)
    return    coeff0 * f0 + coeff1 * f1 + coeff2 * g0 + coeff3 * g1;
end
start = clock()
