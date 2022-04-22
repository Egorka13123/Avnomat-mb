import streamlit as st
from iapws import IAPWS97
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math as M
from sympy import *
#coding:utf-8
st.write("# Задание 1")

st.title('Курсовая работа ')
st.write('Быковский Егор, Фпэ-01-19, Вариант 4')
st.write('Ne = 214 МВт')
st.write('p0 = 12.7 МПа')
st.write('t0 = 544 C')
st.write('ppp = 2-3 МПа')
st.write('tpp = 546 C')
st.write('pk = 4,4 кПа')
st.write('tpv = 234 C')
st.write('z = 9 шт.')

Ne = 214e6
p0 = 12.7e6
t0 = 544
T0 = t0+273.15

age = st.slider('Укажите максимальную границу Pпп', min_value=2.0, max_value=3.0, step=0.1)
age = age + 0.01
P_pp = list(np.arange(2, age, 0.1))
ppp = [p * 1e6 for p in P_pp]
p_pp_min = float(ppp[0])
p_pp_max = float(ppp[-1])

delta_p_0 = 0.05 * p0
delta_p_pp = 0.08 * p_pp_max
delta_p = 0.03 * p_pp_max

tpp = 546
Tpp = tpp+273.15
pk = 4.4e3
tpv = 234
Tpv = tpv+273.15


z = 9

from bokeh.plotting import figure


def Calculate_eta_G0_Gk(N_e, p_0, T_0, p_pp, T_pp, p_k, T_pv):
    point_0 = IAPWS97(P=p_0 * 10 ** (-6), T=T_0)
    s_0 = point_0.s
    h_0 = point_0.h
    v_0 = point_0.v
    p_0_ = p_0 - 0.05 * p_0

    point_p_0_ = IAPWS97(P=p_0_ * 10 ** (-6), h=h_0)
    t_0_ = point_p_0_.T - 273.15
    s_0_ = point_p_0_.s
    v_0_ = point_p_0_.v

    p_1t = p_pp + 0.1 * p_pp
    point_1t = IAPWS97(P=p_1t * 10 ** (-6), s=s_0)
    t_1t = point_1t.T - 273.15
    h_1t = point_1t.h
    v_1t = point_1t.v

    point_pp = IAPWS97(P=p_pp * 10 ** (-6), T=T_pp)
    h_pp = point_pp.h
    s_pp = point_pp.s
    v_pp = point_pp.v

    H_0 = h_0 - h_1t
    eta_oi = 0.85
    H_i_cvd = H_0 * eta_oi

    h_1 = h_0 - H_i_cvd
    point_1 = IAPWS97(P=p_1t * 10 ** (-6), h=h_1)
    s_1 = point_1.s
    T_1 = point_1.T
    v_1 = point_1.v

    p_pp_ = p_pp - 0.03 * p_pp
    point_pp_ = IAPWS97(P=p_pp_ * 10 ** (-6), h=h_pp)
    s_pp_ = point_pp_.s
    v_pp_ = point_pp_.v
    point_kt = IAPWS97(P=p_k * 10 ** (-6), s=s_pp)
    T_kt = point_kt.T
    h_kt = point_kt.h
    v_kt = point_kt.v
    s_kt = s_pp
    H_0_csdcnd = h_pp - h_kt
    eta_oi = 0.85
    H_i_csdcnd = H_0_csdcnd * eta_oi
    h_k = h_pp - H_i_csdcnd
    point_k = IAPWS97(P=p_k * 10 ** (-6), h=h_k)
    T_k = point_k.T
    s_k = point_k.s
    v_k = point_k.v
    point_k_v = IAPWS97(P=p_k * 10 ** (-6), x=0)
    h_k_v = point_k_v.h
    s_k_v = point_k_v.s
    eta_oiI = (h_1 - h_0) / (h_1t - h_0)
    p_pv = 1.4 * p_0
    point_pv = IAPWS97(P=p_pv * 10 ** (-6), T=T_pv)
    h_pv = point_pv.h
    s_pv = point_pv.s
    ksi_pp_oo = 1 - (1 - (T_k * (s_pp - s_k_v)) / ((h_0 - h_1t) + (h_pp - h_k_v))) / (
                1 - (T_k * (s_pp - s_pv)) / ((h_0 - h_1t) + (h_pp - h_pv)))
    # T_0_= IAPWS97(P = p_pv*10**(-6),x = 0).T
    T_0_ = 374.2 + 273.15
    T_ = (point_pv.T - point_k.T) / (T_0_ - point_k.T)
    if T_ <= 0.631818:
        ksi1 = -1.7131 * T_ ** 2 + 2.3617 * T_ - 0.0142
    elif 0.631818 < T_ <= 0.718182:
        ksi1 = -2.5821 * T_ ** 2 + 3.689 * T_ - 0.4825
    elif 0.718182 < T_ <= 0.827273:
        ksi1 = -1.9864 * T_ ** 2 + 3.138 * T_ - 0.3626
    elif 0.827273 < T_ <= 0.936364:
        ksi1 = -2.0619 * T_ ** 2 + 3.3818 * T_ - 0.4814

    if T_ <= 0.631818:
        ksi2 = -1.7131 * T_ ** 2 + 2.3617 * T_ - 0.0142
    elif 0.631818 < T_ <= 0.718182:
        ksi2 = -2.5821 * T_ ** 2 + 3.689 * T_ - 0.4825
    elif 0.718182 < T_ <= 0.827273:
        ksi2 = -1.9864 * T_ ** 2 + 3.138 * T_ - 0.3626
    elif 0.827273 < T_ <= 0.936364:
        ksi2 = -2.0619 * T_ ** 2 + 3.3818 * T_ - 0.4814

    ksi = (ksi1 + ksi2) / 2
    ksi_r_pp = ksi * ksi_pp_oo
    eta_ir = (H_i_cvd + H_i_csdcnd) / (H_i_cvd + (h_pp - h_k_v)) * 1 / (1 - ksi_r_pp)
    H_i = eta_ir * ((h_0 - h_pv) + (h_pp - h_1))
    eta_m = 0.994
    eta_eg = 0.99
    G_0 = N_e / (H_i * eta_m * eta_eg * (10 ** 3))
    G_k = N_e / ((h_k - h_k_v) * eta_m * eta_eg * (10 ** 3)) * (1 / eta_ir - 1)

    return eta_ir, G_0, G_k


eta, G0, Gk = [], [], []

for p in ppp:
        eta_ = Calculate_eta_G0_Gk(N_e=Ne, p_0=p0, T_0=T0, p_pp=p, T_pp=Tpp, p_k=pk, T_pv=Tpv)
        eta.append(eta_[0])
        G0.append(eta_[1])
        Gk.append(eta_[2])


max: float = eta[0]
pos = 0
for i in range(len(eta)):
        if eta[i] > max: max = eta[i]; pos = i

ppp_f = [float(x) * 10 ** (-6) for x in ppp]
eta_f = [float(x) * 100 for x in eta]

ppp_eta = pd.DataFrame({"ppp, МПа": (ppp_f),
                        "eta, %": (eta_f),
                        "G_0, кг/с": (G0),
                        "G_k, кг/с": (Gk)
                        })
st.dataframe(ppp_eta)

ppp__eta = plt.figure()

plt.plot(ppp_f, eta_f)
plt.plot(ppp_f, eta_f, 'ro')
plt.title("Зависимость КПД от давления пром. перегрева")
plt.xlabel("P_пп, MПа")
plt.ylabel("КПД, %")
plt.grid()

st.pyplot(ppp__eta)

p_pp_max = ppp[pos]

fighs = plt.figure()

point_0 = IAPWS97(P=p0*1e-6, T=T0)
p_0_d = p0 - delta_p_0
point_0_d = IAPWS97(P=p_0_d*1e-6, h=point_0.h)
p_1t = p_pp_max + delta_p_pp
point_1t = IAPWS97(P=p_1t*10**(-6), s=point_0.s)
H_01 = point_0.h - point_1t.h
kpd_oi = 0.85
H_i_cvd = H_01 * kpd_oi
h_1 = point_0.h - H_i_cvd
point_1 = IAPWS97(P=p_1t*1e-6, h=h_1)
point_pp = IAPWS97(P=p_pp_max*1e-6, T=Tpp)
p_pp_d = p_pp_max - delta_p_pp
point_pp_d = IAPWS97(P=p_pp_d*1e-6, h=point_pp.h)
point_kt = IAPWS97(P=pk*1e-6, s=point_pp.s)
H_02 = point_pp.h - point_kt.h
kpd_oi = 0.85
H_i_csd_cnd = H_02 * kpd_oi
h_k = point_pp.h - H_i_csd_cnd
point_k = IAPWS97(P=pk*1e-6, h=h_k)

s_0 = [point_0.s-0.05,point_0.s,point_0.s+0.05]
h_0 = [IAPWS97(P = p0*1e-6,s = s_).h for s_ in s_0]
s_1 = [point_0.s-0.05,point_0.s,point_0.s+0.18]
h_1 = [IAPWS97(P=p_1t*1e-6, s = s_).h for s_ in s_1]
s_0_d = [point_0_d.s-0.05, point_0_d.s, point_0_d.s+0.05]
h_0_d = h_0
s_pp = [point_pp.s-0.05,point_pp.s,point_pp.s+0.05]
h_pp = [IAPWS97(P=p_pp_max*1e-6, s=s_).h for s_ in s_pp]
s_k = [point_pp.s-0.05,point_pp.s,point_pp.s+0.8]
h_k = [IAPWS97(P=pk*1e-6, s=s_).h for s_ in s_k]
s_pp_d = [point_pp_d.s-0.05,point_pp_d.s,point_pp_d.s+0.05]
h_pp_d = h_pp

plt.plot([point_0.s,point_0.s,point_0_d.s,point_1.s],[point_1t.h,point_0.h,point_0.h,point_1.h],'-or')
plt.plot([point_pp.s,point_pp.s,point_pp_d.s,point_k.s],[point_kt.h,point_pp.h,point_pp.h,point_k.h],'-or')
plt.plot(s_0,h_0)
plt.plot(s_1,h_1)
plt.plot(s_0_d,h_0_d)
plt.plot(s_pp,h_pp)
plt.plot(s_k,h_k)
plt.plot(s_pp_d,h_pp_d)

for x, y, ind in zip([point_pp.s, point_k.s], [point_pp.h, point_k.h], ['{пп}', '{к}']):
        plt.text(x-0.45, y+40, '$h_' + ind + ' = %.2f $'%y)
for x, y, ind in zip([point_kt.s, point_pp_d.s], [point_kt.h, point_pp_d.h], ['{кт}', '{ппд}']):
        plt.text(x+0.03, y+40, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip ([point_0.s, point_1.s], [point_0.h, point_1.h], ['{0}', '{1}']):
        plt.text(x-0.01, y+120, '$h_' + ind + ' = %.2f $'%y)

for x, y, ind in zip([point_1t.s, point_0_d.s], [point_1t.h, point_0_d.h], ['{1т}', '{0д}']):
        plt.text(x+0.03, y-60, '$h_' + ind + ' = %.2f $'%y)

plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")
plt.grid(True)

st.pyplot(fighs)

st.write("# Задание 2")

def iso_bar(wsp_point, min_s=-0.1, max_s=0.11, step_s=0.011, color = 'r'):
    if not isinstance(wsp_point,list):
        iso_bar_0_s = np.arange(wsp_point.s+min_s,wsp_point.s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point.P, s = i).h for i in iso_bar_0_s]
    else:
        iso_bar_0_s = np.arange(wsp_point[0].s+min_s,wsp_point[1].s+max_s,step_s).tolist()
        iso_bar_0_h = [IAPWS97(P = wsp_point[1].P, s = i).h for i in iso_bar_0_s]
    plt.plot(iso_bar_0_s,iso_bar_0_h,color)

d = 1.1 #m
p_0 = 12.7 #МПа
T_0 = 544+273 #K
n = 50 #Гц
G_0 = 172.54 #кг/с
H_0 = 115 #кДж/кг
rho = 0.05
l_1 = 0.015 #м
alpha_1 = 12 #град
b_1 = 0.06 #м
Delta = 0.003 #м
b_2 = 0.03 #м
kappa_vs = 0

def callculate_optimum(d, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs):
    u = M.pi*d*n
    point_0 = IAPWS97(P = p_0, T = T_0)
    H_0s = H_0*(1-rho)
    H_0r = H_0*rho
    h_1t = point_0.h - H_0s
    point_1t = IAPWS97(h = h_1t, s = point_0.s)
    c_1t = (2000*H_0s)**0.5
    M_1t = c_1t/point_1t.w
    mu_1 = 0.982 - 0.005*(b_1/l_1)
    F_1 = G_0*point_1t.v/mu_1/c_1t
    el_1 = F_1/M.pi/d/M.sin(M.radians(alpha_1))
    e_opt=5*el_1**0.5
    if e_opt > 0.85:
        e_opt = 0.85
    l_1 = el_1/e_opt
    fi_1 = 0.98 - 0.008*(b_1/l_1)
    c_1 = c_1t*fi_1
    alpha_1 = M.degrees(M.asin(mu_1/fi_1*M.sin(M.radians(alpha_1))))
    w_1 = (c_1**2+u**2-2*c_1*u*M.cos(M.radians(alpha_1)))**0.5
    betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1))/(M.cos(M.radians(alpha_1))-u/c_1)))
    Delta_Hs = c_1t**2/2*(1-fi_1**2)
    h_1 = h_1t + Delta_Hs*1e-3
    point_1 = IAPWS97(P = point_1t.P, h = h_1)
    h_2t = h_1 - H_0r
    point_2t = IAPWS97(h = h_2t, s = point_1.s)
    w_2t = (2*H_0r*1e3+w_1**2)**0.5
    l_2 = l_1 + Delta
    mu_2 = 0.965-0.01*(b_2/l_2)
    M_2t = w_2t/point_2t.w
    F_2 = G_0*point_2t.v/mu_2/w_2t
    betta_2 = M.degrees(M.asin(F_2/(e_opt*M.pi*d*l_2)))
    point_1w = IAPWS97(h = point_1.h+w_1**2/2*1e-3, s = point_1.s)
    psi = 0.96 - 0.014*(b_2/l_2)
    w_2 = psi*w_2t
    c_2 = (w_2**2+u**2-2*u*w_2*M.cos(M.radians(betta_2)))**0.5
    alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2))/(M.cos(M.radians(betta_2))-u/w_2)))
    if alpha_2<0:
        alpha_2 = 180 + alpha_2
    Delta_Hr = w_2t**2/2*(1-psi**2)
    h_2 = h_2t+Delta_Hr*1e-3
    point_2 = IAPWS97(P = point_2t.P, h = h_2)
    print(f'H_2 = {point_2t.h:.2f} кДж/кг')
    Delta_Hvs = c_2**2/2
    E_0 = H_0 - kappa_vs*Delta_Hvs
    etta_ol1 = (E_0*1e3 - Delta_Hs-Delta_Hr-(1-kappa_vs)*Delta_Hvs)/(E_0*1e3)
    etta_ol2 = (u*(c_1*M.cos(M.radians(alpha_1))+c_2*M.cos(M.radians(alpha_2))))/(E_0*1e3)
    return etta_ol2, alpha_2

face=plt.figure()

plt.ylabel('eta ')
plt.xlabel('U_cf м/с')

d = [i * 1e-2 for i in list(range(90, 111, 1))]
alpha1 = []
eta = []
ucf = []

for i in d:
    ucf_1 = M.pi * i * n / (2000 * H_0) ** 0.5
    ucf.append(ucf_1)
    eta_ol, alpha = callculate_optimum(i, p_0, T_0, n, G_0, H_0, rho, l_1, alpha_1, b_1, Delta, b_2, kappa_vs)
    print(i, eta_ol, alpha, ucf_1)
    eta.append(eta_ol)
    alpha1.append(alpha)

plt.plot(ucf, eta)
plt.show()


itog=pd.DataFrame({
"d , м": (list(np.arange(0.9, 1.1, 0.01))),
"eta_ol": (eta),
"alfa": (alpha1),
"G_k, кг/с": (ucf)
})


d = 1.1
u = M.pi*d*n
#st.write(f'u = {u:.2f} м/с')
point_0 = IAPWS97(P = p_0, T = T_0)
#st.write(f'h_0 = {point_0.h:.2f} кДж/кг')
#st.write(f's_0 = {point_0.s:.4f} кДж/(кг*К)')
H_0s = H_0*(1-rho)
H_0r = H_0*rho
h_1t = point_0.h - H_0s
#st.write(f'h_1т = {h_1t:.2f} кДж/кг')
point_1t = IAPWS97(h = h_1t, s = point_0.s)
c_1t = (2000*H_0s)**0.5
#st.write(f'c_1т = {c_1t:.2f} м/с')
M_1t = c_1t/point_1t.w
#st.write(f'M_1т = {M_1t:.2f}')
mu_1 = 0.982 - 0.005*(b_1/l_1)
F_1 = G_0*point_1t.v/mu_1/c_1t
#st.write(f'F_1 = {F_1:.4f} м^2')
el_1 = F_1/M.pi/d/M.sin(M.radians(alpha_1))
#st.write(f'el_1 = {el_1:.4f} м')
e_opt=6*el_1**0.5
if e_opt > 0.85:
    e_opt = 0.85
l_1 = el_1/e_opt

fighs=plt.figure()

plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")

def plot_hs_nozzle_t(x_lim, y_lim):
    plt.plot([point_0.s, point_1t.s],[point_0.h, point_1t.h],'ro-')
    iso_bar(point_0,-0.02,0.02,0.001,'c')
    iso_bar(point_1t,-0.02,0.02,0.001,'y')
    plt.xlim(x_lim)
    plt.ylim(y_lim)

plot_hs_nozzle_t([6.5, 6.7],[3300,3500])


#st.write(f'l_1 = {l_1:.4f} м')
if alpha_1 <= 10:
    NozzleBlade = 'C-90-09A'
    t1_ = 0.78
    b1_mod = 6.06
    f1_mod = 3.45
    W1_mod = 0.471
    alpha_inst1 = alpha_1-12.5*(t1_-0.75)+20.2
elif  10 < alpha_1 <= 13:
    NozzleBlade = 'C-90-12A'
    t1_ = 0.78
    b1_mod = 5.25
    f1_mod = 4.09
    W1_mod = 0.575
    alpha_inst1 = alpha_1-10*(t1_-0.75)+21.2
elif  13 < alpha_1 <= 16:
    NozzleBlade = 'C-90-15A'
    t1_ = 0.78
    b1_mod = 5.15
    f1_mod = 3.3
    W1_mod = 0.45
    alpha_inst1 = alpha_1-16*(t1_-0.75)+23.1
else:
    NozzleBlade = 'C-90-18A'
    t1_ = 0.75
    b1_mod = 4.71
    f1_mod = 2.72
    W1_mod = 0.333
    alpha_inst1 = alpha_1-17.7*(t1_-0.75)+24.2
#st.write('Тип профиля:',  NozzleBlade)
#st.write(f'Оптимальный относительный шаг t1_ = {t1_}')
z1 = (M.pi*d)/(b_1*t1_)
z1 = int(z1)
#if z1 % 2 == 0:
    #st.write(f'z1 = {z1}')
#else:
    #z1 = z1+1
   # st.write(f'z1 = {z1}')
t1_ = (M.pi*d)/(b_1*z1)
Ksi_1_ = (0.021042*b_1/l_1 + 0.023345)*100
k_11 = 7.18977510*M_1t**5 - 26.94497258*M_1t**4 + 39.35681781*M_1t**3 - 26.09044664*M_1t**2 + 6.75424811*M_1t + 0.69896998
k_12 = 0.00014166*90**2 - 0.03022881*90 + 2.61549380
k_13 = 13.25474043*t1_**2 - 20.75439502*t1_ + 9.12762245
Ksi_1 = Ksi_1_*k_11*k_12*k_13

fi_1 = M.sqrt(1-Ksi_1/100)
#st.write(f'mu_1 = {mu_1}')
#st.write(f'fi_1 = {fi_1}')

alpha_1 = 12
c_1 = c_1t*fi_1
#st.write(f'c_1 = {c_1:.2f} м/с')
alpha_1 = M.degrees(M.asin(mu_1/fi_1*M.sin(M.radians(alpha_1))))
#st.write(f'alpha_1 = {alpha_1:.2f} град.')
w_1 = (c_1**2+u**2-2*c_1*u*M.cos(M.radians(alpha_1)))**0.5
#st.write(f'w_1 = {w_1}')

c_1u = c_1*M.cos(M.radians(alpha_1))
c_1a = c_1*M.sin(M.radians(alpha_1))
w_1u = c_1u - u
#st.write(c_1u,w_1u)
w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]

face2=plt.figure()



ax = plt.axes()
plt.title("Треугольник скоростей")
ax.arrow(*c_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
plt.text(-2*c_1u/3, -3*c_1a/4, '$c_1$', fontsize=20)
plt.text(-2*w_1u/3, -3*c_1a/4, '$w_1$', fontsize=20)
plt.show()


betta_1 = M.degrees(M.atan(M.sin(M.radians(alpha_1))/(M.cos(M.radians(alpha_1))-u/c_1)))
#st.write(f'betta_1 = {betta_1:.2f}')
Delta_Hs = c_1t**2/2*(1-fi_1**2)
#st.write(f'Delta_Hs = {Delta_Hs:.2f} Дж/кг')

h_1 = h_1t + Delta_Hs*1e-3
#st.write(f'h_1 = {h_1:.2f} кДж/кг')
point_1 = IAPWS97(P = point_1t.P, h = h_1)
h_2t = h_1 - H_0r
#st.write(f'h_2t = {h_2t:.2f} кДж/кг')
point_2t = IAPWS97(h = h_2t, s = point_1.s)
w_2t = (2*H_0r*1e3+w_1**2)**0.5
#st.write(f'w_2t = {w_2t:.2f} м/с')
l_2 = l_1 + Delta
mu_2 = 0.965-0.01*(b_2/l_2)
#st.write(f'mu_2 = {mu_2:.2f}')
M_2t = w_2t/point_2t.w
#st.write(f'M_2t = {M_2t:.2f}')
F_2 = G_0*point_2t.v/mu_2/w_2t
#st.write(f'F_2 = {F_2:.2f}')
betta_2 = M.degrees(M.asin(F_2/(e_opt*M.pi*d*l_2)))
#st.write(f'betta_2 = {betta_2:.2f}')
point_1w = IAPWS97(h = point_1.h+w_1**2/2*1e-3, s = point_1.s)
fig4 = plt.figure(figsize=(10, 10))

face33=plt.figure()

plt.title("h - s диаграмма")
plt.xlabel("s, кДж/(кг*С)")
plt.ylabel("h, кДж/кг")

psi = 0.93
w_2 = psi * w_2t
c_2 = (w_2 ** 2 + u ** 2 - 2 * u * w_2 * M.cos(M.radians(betta_2))) ** 0.5
Delta_Hr = w_2t ** 2 / 2 * (1 - psi ** 2)
h_2 = h_2t + Delta_Hr * 1e-3
point_2 = IAPWS97(P=point_2t.P, h=h_2)
Delta_Hvs = c_2 ** 2 / 2
h_3 = h_2 + Delta_Hvs * 1e-3
point_3 = IAPWS97(P=point_2t.P, h=h_3)
point_2_ = IAPWS97(P=point_2t.P, h=point_0.h-H_0)

def plot_hs_stage_t(x_lim,y_lim):
    plot_hs_nozzle_t(x_lim,y_lim)
    plt.plot([point_0.s,point_1.s],[point_0.h,point_1.h],'bo-')
    plt.plot([point_1.s,point_2t.s],[point_1.h,point_2t.h], 'ro-')
    plt.plot([point_1.s,point_1.s],[point_1w.h, point_1.h],'ro-')

    plt.plot([point_1.s, point_2.s], [point_1.h, point_2.h], 'bo-')
    plt.plot([point_2.s, point_3.s], [point_2.h, point_3.h], 'bo-')

    iso_bar(point_2t,-0.02,0.02,0.001,'y')
    iso_bar(point_1w,-0.005,0.005,0.001,'c')
plot_hs_stage_t([6.575,6.649],[3300,3550])
plt.show()

if betta_2 <= 15:
    RotorBlade = 'P-23-14A'
    t2_ = 0.63
    b2_mod = 2.59
    f2_mod = 2.44
    W2_mod = 0.39
    beta_inst2 = betta_2-12.5*(t2_-0.75)+20.2

elif  15 < betta_2 <= 19:
    RotorBlade = 'P-26-17A'
    t2_ = 0.65
    b2_mod = 2.57
    f2_mod = 2.07
    W2_mod = 0.225
    beta_inst2 = betta_2-19.3*(t2_-0.6)+60


elif  19 < betta_2 <= 23:
    RotorBlade = 'P-30-21A'
    t2_ = 0.63
    b2_mod = 2.56
    f2_mod = 1.85
    W2_mod = 0.234
    beta_inst2 = betta_2-12.8*(t2_-0.65)+58


elif 23 < betta_2 <= 27:
    RotorBlade = 'P-35-25A'
    t2_ = 0.6
    b2_mod = 2.54
    f2_mod = 1.62
    W2_mod = 0.168
    beta_inst2 = betta_2-16.6*(t2_-0.65)+54.3

elif 27 < betta_2 <= 31:
    RotorBlade = 'P-46-29A'
    t2_ = 0.51
    b2_mod = 2.56
    f2_mod = 1.22
    W2_mod = 0.112
    beta_inst2 = betta_2-50.5*(t2_-0.6)+47.1


else:
    RotorBlade = 'P-50-33A'
    t2_ = 0.49
    b2_mod = 2.56
    f2_mod = 1.02
    W2_mod = 0.079
    beta_inst2 = betta_2-20.8*(t2_-0.6)+43.7

#st.write('Тип профиля:',  RotorBlade)
#st.write(f'Оптимальный относительный шаг t2_ = {t2_}')

z2 = int((M.pi*d)/(b_2*t2_))
#st.write(f'z2 = {z2}')
t2_ = (M.pi*d)/(b_2*z2)
Ksi_2_ = 4.364*b_2/l_2 + 4.22
k_21 = -13.79438991*M_2t**4 + 36.69102267*M_2t**3 - 32.78234341*M_2t**2 + 10.61998662*M_2t + 0.28528786
k_22 = 0.00331504*betta_1**2 - 0.21323910*betta_1 + 4.43127194
k_23 = 60.72813684*t2_**2 - 76.38053189*t2_ + 24.97876023
Ksi_2 = Ksi_2_*k_21*k_22*k_23

psi = M.sqrt(1-Ksi_2/100)
#st.write(f'psi = {psi:.2f}')


psi = 0.91

w_2 = psi*w_2t
#st.write(f'w_2 = {w_2:.2f} м/с')
c_2 = (w_2**2+u**2-2*u*w_2*M.cos(M.radians(betta_2)))**0.5
#st.write(f'c_2 = {c_2:.2f} м/с')
alpha_2 = M.degrees(M.atan(M.sin(M.radians(betta_2))/(M.cos(M.radians(betta_2))-u/w_2)))
#st.write(f'alpha_2 = {alpha_2:.2f}')
Delta_Hr = w_2t**2/2*(1-psi**2)
#st.write(f'Delta_Hr = {Delta_Hr:.2f} Дж/кг')
h_2 = h_2t+Delta_Hr*1e-3
point_2 = IAPWS97(P = point_2t.P, h = h_2)
Delta_Hvs = c_2**2/2
#st.write(f'Delta_Hvs = {Delta_Hvs:.2f} Дж/кг')
E_0 = H_0 - kappa_vs*Delta_Hvs
etta_ol1 = (E_0*1e3 - Delta_Hs-Delta_Hr-(1-kappa_vs)*Delta_Hvs)/(E_0*1e3)
#st.write(f'1. etta_ol = {etta_ol1}')
etta_ol2 = (u*(c_1*M.cos(M.radians(alpha_1))+c_2*M.cos(M.radians(alpha_2))))/(E_0*1e3)
#st.write(f'2. etta_ol = {etta_ol2}')

c_1u = c_1*M.cos(M.radians(alpha_1))
c_1a = c_1*M.sin(M.radians(alpha_1))
w_1u = c_1u - u
w_2a = w_2*M.sin(M.radians(betta_2))
w_2u = w_2*M.cos(M.radians(betta_2))
c_2u=w_2u-u
#st.write(c_1u,w_1u)
w_1_tr = [0, 0, -w_1u, -c_1a]
c_1_tr = [0, 0, -c_1u, -c_1a]
u_1_tr = [-w_1u, -c_1a, -u, 0]

w_2_tr = [0, 0, w_2u, -w_2a]
c_2_tr = [0, 0, c_2u, -w_2a]
u_2_tr = [w_2u,-w_2a, -u, 0]

face3=plt.figure()

ax = plt.axes()
ax.arrow(*c_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_1_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
ax.arrow(*c_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='r', ec='r')
ax.arrow(*w_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='b', ec='b')
ax.arrow(*u_2_tr, head_width=5, length_includes_head = True,head_length=20, fc='g', ec='g')
plt.text(-2*c_1u/3, -3*c_1a/4, '$c_1$', fontsize=20)
plt.text(-2*w_1u/3, -3*c_1a/4, '$w_1$', fontsize=20)
plt.text(2*c_2u/3, -3*w_2a/4, '$c_2$', fontsize=20)
plt.text(2*w_2u/3, -3*w_2a/4, '$w_2$', fontsize=20)
plt.show()
plt.title("Треугольник скоростей")


delta_a = 0.0025
z_per_up = 2
mu_a = 0.5
mu_r = 0.75
d_per = d + l_1
delta_r = d_per*0.001
delta_ekv = 1/M.sqrt(1/(mu_a*delta_a)**2+z_per_up/(mu_r*delta_r)**2)
#st.write("""Эквивалентный зазор в уплотнении по бандажу (периферийном)
           #streamlit run C:/Users/bykov/PycharmProjects/pythonProject4/main.py [ARGUMENTS] delta_ekv = %.3f мм""" % (delta_ekv*1000))
xi_u_b=M.pi*d_per*delta_ekv*etta_ol1/F_1*M.sqrt(rho+1.8*l_2/d)
#st.write("""Относительные потери от утечек через бандажные уплотнения
          # xi_u_b = %.3f""" % xi_u_b)
Delta_Hub = xi_u_b*E_0
#st.write("""Абсолютные потери от утечек через периферийное уплотнение ступени
           # Delta_Hub = %.3f кДж/кг""" % Delta_Hub)

k_tr=0.0007
Kappa_VS = 0
u = M.pi*d*n
c_f = M.sqrt(2000*H_0)
ucf = u/c_f
xi_tr = k_tr*d**2/F_1*ucf**3
#st.write("""Определяем u/c_ф для ступени
          #  U/c_ф = %.3f""" % ucf)
#st.write("""Относительные потери от трения диска
           # xi_tr =  = %.5f""" % xi_tr)
Delta_Htr = xi_tr*E_0
#st.write("""Абсолютные потери от трения диска
            #Delta_Htr = %.3f кДж/кг""" % Delta_Htr)

k_v = 0.065
m = 1
xi_v = k_v/M.sin(M.radians(alpha_1))*(1-e_opt)/e_opt*ucf**3*m
#st.write('Относительные вентиляционные потери',xi_v)
i_p = 4
B_2 = b_2*M.sin(M.radians(beta_inst2))
xi_segm = 0.25*B_2*l_2/F_1*ucf*etta_ol1*i_p
#st.write('Относительные сегментные потери',xi_segm)
xi_parc = xi_v+xi_segm
Delta_H_parc = E_0*xi_parc

H_i = E_0 - Delta_Hr*1e-3 - Delta_Hs*1e-3 - (1-Kappa_VS)*Delta_Hvs*1e-3 - Delta_Hub - Delta_Htr - Delta_H_parc
#st.write("""Использованный теплоперепад ступени
          # H_i = %.3f кДж/кг$""" % H_i)
eta_oi = H_i/E_0
st.write("""Внутренний относительный КПД ступени  
        eta_oi  = %.3f $""" % eta_oi)
N_i = G_0*H_i
st.write("""Внутреняя мощность ступени  
            N_i = = %.2f кВт""" % N_i)
itog
st.pyplot(face)
st.pyplot(face33)
st.pyplot(face3)

# Задание 3




st.write("# Задание 3")

st.write(" *Исходные данные:* ")



drs = 1.1
P0 = 12.7
h0 = 3459
G0 = 171.2003
etaoi = 0.743
Z = 9
Pz =(p_pp_max + 0.1*p_pp_max)*10**(-6)


deltaD = st.slider('deltaD', min_value=0.20, max_value=0.26, step=0.01)  # m
n = 50  # Гц
rho_s = st.slider('rho_s', min_value=0.03, max_value=0.07, step=0.01)
alfa = st.slider('alfa', min_value=10, max_value=16, step=1)  # град
fi = st.slider('fi', min_value=0.93, max_value=0.96, step=0.01)
mu1 = st.slider('𝜇', min_value=0.95, max_value=0.97, step=0.01)
delta = 0.003
tetta = 20

st.write("""# """)
st.write(" *Дано:* ")
st.write(""" P0 = """ + str(P0) + """ МПа""")
st.write(""" h0 = """ + str(h0) + """ кДж/кг""")
st.write(""" dрс = """ + str(drs) + """ м """)
st.write(""" Z = """ + str(Z) + """ шт """)
st.write(""" G0 = """ + str(G0) + """ кг/с """)
st.write(""" n = """ + str(n) + """ Гц """)
st.write(""" eta_oi = """ + str(etaoi) + """ """)
st.write(""" Pz = """ + str(Pz) + """ МПа """)



D1 = drs - deltaD
sat_steam = IAPWS97(P=P0, h=h0)
s_0 = sat_steam.s
t_0 = sat_steam.T
error = 2
i = 1
while error > 0.5:
    rho = rho_s + 1.8 / (tetta + 1.8)
    X = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho))
    H01 = 12.3 * (D1 / X) ** 2 * (n / 50) ** 2
    h2t = h0 - H01
    steam2t = IAPWS97(h=h2t, s=s_0)
    v2t = steam2t.v
    l11 = G0 * v2t * X / (M.pi ** 2 * D1 ** 2 * n * M.sqrt(1 - rho) * M.sin(M.radians(alfa)) * mu1)
    tetta_old = tetta
    tetta = D1 / l11
    error = abs(tetta - tetta_old) / tetta_old * 100
    i += 1

l21 = l11 + delta
d_s = D1 - l21
steam_tz = IAPWS97(P=Pz, s=s_0)
h_zt = steam_tz.h
H0 = h0 - h_zt
Hi = H0 * etaoi
h_z = h0 - Hi
steam_z = IAPWS97(P=Pz, h=h_z)
v_2z = steam_z.v
x = Symbol('x')
с = solve(x ** 2 + x * d_s - (l21 * (d_s + l21) * v_2z / v2t))
for j in с:
    if j > 0:
        l2z = j
d2z = d_s + l2z
tetta1 = (l21 + d_s) / l21
tettaz = (l2z + d_s) / l2z
rho1 = rho_s + 1.8 / (1.8 + tetta1)
rhoz = rho_s + 1.8 / (1.8 + tettaz)
X1 = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho1))
Xz = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rhoz))

DeltaZ = 1
ite = 0
while DeltaZ > 0:
    matr = []
    Num = 0
    SumH = 0
    for _ in range(int(Z)):
        li = (l21 - l2z) / (1 - Z) * Num + l21
        di = (D1 - d2z) / (1 - Z) * Num + D1
        tetta_i = di / li
        rho_i = rho_s + 1.8 / (1.8 + tetta_i)
        X_i = (fi * M.cos(M.radians(alfa))) / (2 * M.sqrt(1 - rho_i))
        if Num < 1:
            H_i = 12.3 * (di / X_i) ** 2 * (n / 50) ** 2
        else:
            H_i = 12.3 * (di / X_i) ** 2 * (n / 50) ** 2 * 0.95
        Num = Num + 1
        H_d = 0
        SumH = SumH + H_i
        matr.append([Num, round(di, 3), round(li, 3), round(tetta_i, 2), round(rho_i, 3), round(X_i, 3), round(H_i, 2),
                     round(H_d, 2)])
    H_m = SumH / Z
    q_t = 4.8 * 10 ** (-4) * (1 - etaoi) * H0 * (Z - 1) / Z
    Z_new = round(H0 * (1 + q_t) / H_m)
    DeltaZ = abs(Z - Z_new)
    # print(ite, Z)
    Z = Z_new
    ite += 1
DeltaH = (H0 * (1 + q_t) - SumH) / Z
a = 0
for elem in matr:
    matr[a][7] = round(elem[6] + DeltaH, 2)
    a += 1

N_ = []
di_ = []
li_ = []
tettai_ = []
rhoi_ = []
Xi_ = []
Hi_ = []
Hdi_ = []
a = 0
for elem in matr:
    N_.append(matr[a][0])
    di_.append(matr[a][1])
    li_.append(matr[a][2])
    tettai_.append(matr[a][3])
    rhoi_.append(matr[a][4])
    Xi_.append(matr[a][5])
    Hi_.append(matr[a][6])
    Hdi_.append(matr[a][7])
    a += 1

di_ = [float(x) for x in di_]
li_ = [float(x) for x in li_]
tettai_ = [float(x) for x in tettai_]
rhoi_ = [float(x) for x in rhoi_]
Xi_ = [float(x) for x in Xi_]
Hi_ = [float(x) for x in Hi_]
Hdi_ = [float(x) for x in Hdi_]

## Таблица
table = pd.DataFrame({"№ ступени": (N_),
                      "di, м": (di_),
                      "li, м": (li_),
                      "θi ": (tettai_),
                      "ρi ": (rhoi_),
                      "Xi ": (Xi_),
                      "Hi, кДж/кг": (Hi_),
                      "Hi + Δ, кДж/кг": (Hdi_)
                      }
                     )

st.dataframe(table)

## Графики
z = []
for a in range(1, Z + 1):
    z.append(a)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, di_, '-ro')
plt.title('Рисунок 1. Распределение средних диаметров по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, li_, '-ro')
plt.title('Рисунок 2. Распределение высот лопаток по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, tettai_, '-ro')
plt.title('Рисунок 3. Распределение обратной веерности по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, rhoi_, '-ro')
plt.title('Рисунок 4. Распределение степени реактивности по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, Xi_, '-ro')
plt.title('Рисунок 5. Распределение U/Cф по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, Hi_, '-ro')
plt.title('Рисунок 6. Распределение теплоперепадов по проточной части')
st.pyplot(fig)

st.write("#")
fig = plt.figure(figsize=(10, 5))
ax = fig.gca()
ax.set_xticks(np.arange(1, 16, 1))
plt.grid(True)
plt.plot(z, Hdi_, '-ro')
plt.title('Рисунок 7. Распределение теплоперепадов с учетом невязки по проточной части')
st.pyplot(fig)

#####
import pythoncom
pythoncom.CoInitializeEx(0)
import ezdxf


d_отв_ц=100;

rcp=270;

r_к=410; #Корневой диаметр

r_r=160; #Диаметр ротора

r_к_рег=502; #Диаметр корня регулирующей ступени

l_конц_п=480; #Длина переднего концевого уплотнения

l_сред=1440;

b_1_рег=60;# Хорда сопловой лопатки регулирующей ступени

lпром=120;

l_средконец=660;

b2_рег=120;

l_конц_3=240;

l_конц_3пр=960;

l_конц_пр=360;

r_rпр=180;

Constx=rcp,r_к,r_к,rcp

Consty=l_конц_п+l_сред

c=b_1_рег+b2_рег

d=l_средконец+l_конц_3

x = [0,0,l_конц_п,l_конц_п,Consty,Consty,Consty+b_1_рег,Consty+b_1_рег,Consty+c,Consty+c,Consty+c+b_1_рег,Consty+c+b_1_рег,Consty+c+c,Consty+c+c,Consty+c+c+b_1_рег,Consty+c+c+b_1_рег,Consty + 3*c,Consty+3*c,Consty+3*c+b_1_рег,Consty+3*c+b_1_рег,Consty+4*c,Consty+4*c,Consty+4*c+b_1_рег,Consty+4*c+b_1_рег,Consty+5*c,Consty+5*c,Consty+5*c+b_1_рег,Consty+5*c+b_1_рег,Consty+5*c+b_1_рег+l_средконец,Consty+5*c+b_1_рег+l_средконец,Consty+5*c+b_1_рег+l_средконец+b2_рег,Consty+5*c+b_1_рег+l_средконец+b2_рег,Consty+6*c+d,Consty+6*c+d,Consty+6*c+d+b_1_рег,Consty+6*c+d+b_1_рег,Consty+7*c+d,Consty+7*c+d,Consty+7*c+d+b_1_рег,Consty+7*c+d+b_1_рег,Consty+8*c+d,Consty+8*c+d,Consty+8*c+d+b_1_рег,Consty+8*c+d+b_1_рег,Consty+9*c+d,Consty+9*c+d,Consty+9*c+d+b_1_рег,Consty+9*c+d+b_1_рег,Consty+10*c+d,Consty+10*c+d,Consty+10*c+d+b_1_рег,Consty+10*c+d+b_1_рег,Consty+11*c+d,Consty+11*c+d,Consty+11*c+d+b_1_рег,Consty+11*c+d+b_1_рег,Consty+11*c+d+b_1_рег+l_конц_3пр,Consty+11*c+d+b_1_рег+l_конц_3пр,l_конц_пр+Consty+11*c+d+b_1_рег+l_конц_3пр,l_конц_пр+Consty+11*c+d+b_1_рег+l_конц_3пр]
y = [0,r_r,r_r,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к_рег,r_к_рег,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_к,r_к,rcp,rcp,r_rпр,r_rпр,0]
def Convert(lst):
    return [ -i for i in lst ]
y2 = Convert(y)
fig = plt.figure()
plt.plot(x, y,'k-')
plt.plot(x,y2,'k-')


plt.xlim([-30,6210])
plt.ylim([-3050,3050])
st.pyplot(fig)




from pyautocad import Autocad, APoint
doc = ezdxf.new('R2010')
msp = doc.modelspace()
with st.expander("Открыть чертеж в формате DWG (требуется AutoCAD любой версии)"):
    st.write("""
        Пожалуйста нажмите на галочку "Шаг 1" и подождите открытия AutoCAD, после нажмите на галочку "Шаг 2"
    """)
    agree = st.checkbox('Шаг 1')

    if agree:
        st.write('')
        acad = Autocad(create_if_not_exists=True)
        agree2 = st.checkbox('Шаг 2')
        print(acad.doc.Name)
        if agree2:

            for i in range(0, len(x) - 1):
                acad.model.AddLine(APoint(x[i], y[i]), APoint(x[i + 1], y[i + 1]))
            for i in range(0, len(x) - 1):
                acad.model.AddLine(APoint(x[i], y2[i]), APoint(x[i + 1], y2[i + 1]))

for i in range(0,len(x)-1):
    msp.add_line((x[i],y[i]),(x[i+1],y[i+1] ))

for i in range(0,len(x)-1):
    msp.add_line((x[i],y2[i]),(x[i+1],y2[i+1] ))
msp.add_line((-30,0),(6210,0))
doc.saveas('Rotor.dxf')

with open("Rotor.dxf", "r+") as file:
    btn = st.download_button(
        label="Скачать чертеж в формате DXF",
        data=file,
        file_name="Rotor.dxf",
        mime="DXF/dxf"
    )


