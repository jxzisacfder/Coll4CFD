# L_W格式，对本问题失效
import numpy as np
import matplotlib.pyplot as plt
import function as fc
# constants

N = 1000  # 划分数
cfl = 0.1  # 科朗数
xmax = 1
xmin = -1
x = np.linspace(xmin, xmax, N+1)  # 划分x
dx = (xmax-xmin)/N
x_dis = 0.5  # 不连续点位置
n_dis = int(0.5*N)  # 不连续点下标
t_max = 0.2  # 时间
gamma = 1.4
# param0 = [0.445, 0.698, 3.528, 0.5, 0, 0.571]  # 前三个为左侧气体的密度速度压强，后三个为右侧气体对应项
param0 = [1, 0, 2.5, 0.125, 0.0, 0.25]
flag = 1
steps = 0
maxsteps = 1e5
current_time = 0

# initial

u_tt = u_0 = fc.euler_generate(N, n_dis, *param0)  # 初始化原始量
u_vector = fc.physics_to_conservation(u_tt)  # 把原始量转化为守恒量
# u_vector = np.column_stack((u_vector[:,0], u_vector, u_vector[:,-1]))
u_t = u_vector  # 不同时间的守恒量
dF = np.zeros((3, N+1))  # flux的差分
c = np.zeros((1, N+3))
# L-W discrete function


def fluxdiffL_W(dt, dx, u, F):
    dF = 0.5*(F[:, 2:]-F[:, :-2])/dx  # 基础差分部分

    # 计算Jacobi矩阵

    u_l = np.zeros((3, N+1))
    J_L = np.zeros((N+1, 3, 3))
    u_l = 0.5*(u[:, :-2]+u[:, 1:-1])  # 计算Jacobi矩阵的自变量
    J_L[:, 0, 0] = 0
    J_L[:, 0, 1] = 1
    J_L[:, 0, 2] = 0
    J_L[:, 1, 0] = 0.5*(gamma-3)*pow(u_l[1, :]/u_l[0, :], 2)
    J_L[:, 1, 1] = (3-gamma)*u_l[1, :]/u_l[0, :]
    J_L[:, 1, 2] = gamma-1
    J_L[:, 2, 0] = (-gamma)*u_l[1, :]*u_l[2, :]/pow(u_l[0, :],
                                                    2)+(gamma-1)*pow(u_l[1, :]/u_l[0, :], 3)
    J_L[:, 2, 1] = gamma*u_l[2, :]/u_l[0, :] - \
        (gamma-1)*1.5*pow(u_l[1, :]/u_l[0, :], 2)
    J_L[:, 2, 2] = gamma*u_l[1, :]/u_l[0, :]

    u_r = np.zeros((3, N+1))
    J_R = np.zeros((N+1, 3, 3))
    u_r = 0.5*(u[:, 1:-1]+u[:, 2:])
    J_R[:, 0, 0] = 0
    J_R[:, 0, 1] = 1
    J_R[:, 0, 2] = 0
    J_R[:, 1, 0] = 0.5*(gamma-3)*pow(u_r[1, :]/u_r[0, :], 2)
    J_R[:, 1, 1] = (3-gamma)*u_r[1, :]/u_r[0, :]
    J_R[:, 1, 2] = gamma-1
    # (-gamma)*u_r[1,:]*u_r[2,:]/pow(u_r[0,:],2)+(gamma-1)*pow(u_r[1,:]/u_r[0,:],3)
    J_R[:, 2, 0] = (gamma-1)*pow(u_r[1, :]/u_r[0, :], 3) - \
        gamma*u_r[1, :]*u_r[2, :]/pow(u_r[0, :], 2)
    J_R[:, 2, 1] = gamma*u_r[2, :]/u_r[0, :] - \
        (gamma-1)*1.5*pow(u_r[1, :]/u_r[0, :], 2)
    J_R[:, 2, 2] = gamma*u_r[1, :]/u_r[0, :]

    # 计算修正量

    dF2 = np.zeros((3, N+1))
    for i in range(1, N+2):
        j_l = J_L[i-1, :, :]
        j_r = J_R[i-1, :, :]
        mid_l = np.dot(j_l, (F[:, i]-F[:, i-1]))
        mid_r = np.dot(j_r, (F[:, i+1]-F[:, i]))
        dF2[:, i-1] = 0.5*dt/(dx*dx)*(mid_r-mid_l)

    dF = dF-dF2
    return dF


def fluxDiscre(dt, dx, u_conservation):
    u_physics = fc.consevation_to_physics(u_conservation)  # 守恒量转化为基本量

    F = np.array([u_physics[0, :]*u_physics[1, :], u_physics[0, :]*u_physics[1, :]*u_physics[1, :] +
                  u_physics[2, :], (u_conservation[2, :]+u_physics[2, :])*u_physics[1, :]])  # flux,注意中间有个u_conservation
    dF = fluxdiffL_W(dt, dx, u_conservation, F)

    return dF


def weird_division(n, d):
    return n / d if d else 0


# execute
'''
c = np.sqrt(gamma*u_vector[2, :]/u_vector[0, :])  # 声速
dt = cfl*dx/max(abs(u_vector[1, :])+c)  # 时间步长
dF = fluxDiscre(dt, dx, u_vector)
'''

while flag:
    # 初始化dt以及判断
    steps = steps+1
    u_0 = fc.consevation_to_physics(np.column_stack(
        (u_vector[:, 0], u_vector, u_vector[:, -1])))
    for i in range(N+2):
        c[0, i] = np.sqrt(weird_division(gamma*u_0[2, i], u_0[0, i]))  # 声速
    aa = max(abs(u_0[1, :])+c)
    dt = weird_division(cfl*dx, aa[0])  # 时间步长
    #dt = 0.1*dx
    if steps > maxsteps:
        print('warning: maxsteps reach')
        break

    if current_time+dt >= t_max:
        dt = t_max-current_time
        flag = 0

    current_time = current_time+dt
    '''
    # 先随时间推进，再计算边界，再计算dF,更新u_vector
    u_t = u_vector[:, 1:-1]-dt*dF  # 推进
    u_t = np.column_stack((u_t[:, 0], u_t, u_t[:, -1]))
    dF = fluxDiscre(dt, dx, u_t)  # 计算dF
    u_vector = u_t.copy()  # 更新
    '''
    u_t = np.column_stack((u_t[:, 0], u_t, u_t[:, -1]))
    dF = fluxDiscre(dt, dx, u_t)
    u_vector = u_t[:, 1:-1]-dt*dF+0.04 * \
        (u_t[:, 2:]-2*u_t[:, 1:-1]+u_t[:, :-2])
    u_t = u_vector.copy()

u_t = np.column_stack((u_t[:, 0], u_t, u_t[:, -1]))
uu = fc.consevation_to_physics(u_t)

plt.figure(1)
plt.style.use("ggplot")
plt.plot(x, uu[0, 1:-1], '-b')
plt.xlabel('x')
plt.ylabel('density')

plt.figure(2)
plt.style.use("ggplot")
plt.plot(x, uu[1, 1:-1], '-b')
plt.xlabel('x')
plt.ylabel('velocity')

plt.figure(3)
plt.style.use("ggplot")
plt.plot(x, uu[2, 1:-1], '-b')
plt.xlabel('x')
plt.ylabel('pressure')

plt.show()
