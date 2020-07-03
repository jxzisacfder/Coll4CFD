
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

    u_physics = consevation_to_physics(u_conservation)  # 守恒量转化为基本量
    E = u_physics[2, :]/(gamma-1)+0.5*u_physics[0, :] * \
        u_physics[1, :]*u_physics[1, :]  # 能量
    F = np.array([u_physics[0, :]*u_physics[1, :], u_physics[0, :]*u_physics[1, :]
                  * u_physics[1, :]+u_physics[2, :], (E+u_physics[2, :])*u_physics[1, :]])  # flux
    dF = 0.5*(F[:, 2:]-F[:, :-2])/dx  # 差分

    return dF


def fluxDiscre1(dx, u_conservation):

    u_physics = consevation_to_physics1(u_conservation)  # 守恒量转化为基本量
    E = u_physics[2, :]/(gamma-1)+0.5*u_physics[0, :] * \
        u_physics[1, :]*u_physics[1, :]  # 能量
    F = np.array([u_physics[0, :]*u_physics[1, :], u_physics[0, :]*u_physics[1, :]
                  * u_physics[1, :]+u_physics[2, :], (E+u_physics[2, :])*u_physics[1, :]])  # flux
    dF = (F[:, 2:]-F[:, 1:-1])/dx  # 差分

    return dF


def fluxDiscre2(dx, u_conservation):

    u_physics = consevation_to_physics1(u_conservation)  # 守恒量转化为基本量
    E = u_physics[2, :]/(gamma-1)+0.5*u_physics[0, :] * \
        u_physics[1, :]*u_physics[1, :]  # 能量
    F = np.array([u_physics[0, :]*u_physics[1, :], u_physics[0, :]*u_physics[1, :]
                  * u_physics[1, :]+u_physics[2, :], (E+u_physics[2, :])*u_physics[1, :]])  # flux
    dF = (F[:, 1:-1]-F[:, :-2])/dx  # 差分

    return dF


def consevation_to_physics1(u_conservation):
    u_prim = np.zeros((3, N+1))
    u_prim[0, :] = u_conservation[0, :]
    u_prim[1, :] = u_conservation[1, :]/u_conservation[0, :]
    u_prim[2, :] = (u_conservation[2, :]-0.5*u_conservation[1, :]
                    * u_conservation[1, :]/u_conservation[0, :])*(gamma-1)

    return u_prim


def consevation_to_physics(u_conservation):
    u_prim = np.zeros((3, N+3))
    u_prim[0, :] = u_conservation[0, :]
    u_prim[1, :] = u_conservation[1, :]/u_conservation[0, :]
    u_prim[2, :] = (u_conservation[2, :]-0.5*u_conservation[1, :]
                    * u_conservation[1, :]/u_conservation[0, :])*(gamma-1)

    return u_prim
