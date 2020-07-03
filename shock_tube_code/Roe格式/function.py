
import numpy as np
gamma = 1.4
N = 1000


def euler_generate(points, dis_location, *param):
    initial_l = param[:3]
    initial_r = param[3:]
    counts = points-dis_location+1
    rho = np.zeros(N+1)
    u = np.zeros(N+1)
    p = np.zeros(N+1)
    rho[:dis_location+1] = initial_l[0]
    rho[dis_location+1:] = initial_r[0]
    u[:dis_location+1] = initial_l[1]
    u[dis_location+1:] = initial_r[1]
    p[:dis_location+1] = initial_l[2]
    p[dis_location+1:] = initial_r[2]
    u = np.array([rho, u, p])

    return u


def physics_to_conservation(u_0):
    u = np.zeros((3, N+1))
    u[0, :] = u_0[0, :]
    u[1, :] = u_0[0, :]*u_0[1, :]
    u[2, :] = u_0[2, :]/(gamma-1)+0.5*u_0[0, :]*u_0[1, :]*u_0[1, :]

    return u


def fluxDiscre(dx, u_conservation):

    # u_physics = consevation_to_physics(u_conservation)  # 守恒量转化为基本量
    # 计算u_l和u_r，一阶插值

    U_l = u_conservation[:, 1:-2]
    U_r = u_conservation[:, 2:-1]

    # 计算Roe平均

    rho_l = U_l[0, :]
    u_l = U_l[1, :]/U_l[0, :]
    e_l = U_l[2, :]/U_l[0, :]

    p_l = (gamma-1)*(e_l-0.5*u_l*u_l)*rho_l
    H_l = e_l+p_l/rho_l

    rho_r = U_r[0, :]
    u_r = U_r[1, :]/U_r[0, :]
    e_r = U_r[2, :]/U_r[0, :]
    p_r = (gamma-1)*(e_r-0.5*u_r*u_r)*rho_r
    H_r = e_r+p_r/rho_r

    rho_roe = pow((0.5*(np.sqrt(rho_l)+np.sqrt(rho_r))), 2)
    u_roe = (np.sqrt(rho_l)*u_l+np.sqrt(rho_r)*u_r) / \
        (np.sqrt(rho_l)+np.sqrt(rho_r))
    H_roe = (np.sqrt(rho_l)*H_l+np.sqrt(rho_r)*H_r) / \
        ((np.sqrt(rho_l)+np.sqrt(rho_r)))
    p_roe = (gamma-1)/gamma*(rho_roe*H_roe-0.5*rho_roe*u_roe*u_roe)
    c_roe = np.sqrt((gamma-1)*(H_roe-0.5*u_roe*u_roe))

    # 计算特征分解矩阵

    lamda = np.zeros((N+2, 3, 3))
    S = np.zeros((N+2, 3, 3))
    S_inv = np.zeros((N+2, 3, 3))
    A_abs = np.zeros((N+2, 3, 3))
    F = np.zeros((3, N+2))

    lamda[:, 0, 0] = u_roe
    lamda[:, 1, 1] = u_roe-c_roe
    lamda[:, 2, 2] = u_roe+c_roe

    S[:, 0, 0] = 1
    S[:, 0, 1] = 1
    S[:, 0, 2] = 1
    S[:, 1, 0] = u_roe
    S[:, 1, 1] = u_roe-c_roe
    S[:, 1, 2] = u_roe+c_roe
    S[:, 2, 0] = 0.5*u_roe*u_roe
    S[:, 2, 1] = 0.5*u_roe*u_roe+gamma/(gamma-1)*p_roe/rho_roe-u_roe*c_roe
    S[:, 2, 2] = 0.5*u_roe*u_roe + gamma/(gamma-1)*p_roe/rho_roe + u_roe*c_roe

    S_inv[:, 0, 0] = 1 - (gamma-1)*u_roe*u_roe/(2*c_roe*c_roe)
    S_inv[:, 0, 1] = (gamma-1)*u_roe/(c_roe*c_roe)
    S_inv[:, 0, 2] = -(gamma-1)/(c_roe*c_roe)
    S_inv[:, 1, 0] = (gamma-1)*u_roe*u_roe/(4*c_roe*c_roe) + 0.5*u_roe/c_roe
    S_inv[:, 1, 1] = -(gamma-1)*u_roe/(2*c_roe*c_roe) - 1/(2*c_roe)
    S_inv[:, 1, 2] = (gamma-1)/(2*c_roe*c_roe)
    S_inv[:, 2, 0] = (gamma-1)*u_roe*u_roe/(4*c_roe*c_roe) - 0.5*u_roe/c_roe
    S_inv[:, 2, 1] = -(gamma-1)*u_roe/(2*c_roe*c_roe) + 1/(2*c_roe)
    S_inv[:, 2, 2] = (gamma-1)/(2*c_roe*c_roe)

    for i in range(N+2):
        aa = np.dot(S[i, :, :], abs(lamda[i, :, :]))
        A_abs[i, :, :] = np.dot(aa, S_inv[i, :, :])

    # 计算半点处重构函数通量

    F_l = np.array([rho_l*u_l, rho_l*u_l*u_l + p_l, u_l*(rho_l*e_l+p_l)])
    F_r = np.array([rho_r*u_r, rho_r*u_r*u_r + p_r, u_r*(rho_r*e_r+p_r)])

    for i in range(N+2):
        F[:, i] = 0.5*(F_l[:, i]+F_r[:, i] -
                       np.dot(A_abs[i, :, :], (U_r[:, i]-U_l[:, i])))

    # 结果

    dF = (F[:, 1:]-F[:, :-1])/dx  # 差分

    return dF


def consevation_to_physics(u_conservation):
    u_prim = np.zeros((3, N+1))
    u_prim[0, :] = u_conservation[0, :]
    u_prim[1, :] = u_conservation[1, :]/u_conservation[0, :]
    u_prim[2, :] = (u_conservation[2, :]-0.5*u_conservation[1, :]
                    * u_conservation[1, :]/u_conservation[0, :])*(gamma-1)

    return u_prim
