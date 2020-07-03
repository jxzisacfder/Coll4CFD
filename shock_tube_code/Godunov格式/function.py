import numpy as np

N = 100
gamma = 1.4


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


'''
def fluxDiscre(dx, dt, u_conservation):

    u_physics1 = consevation_to_physics1(u_conservation[0, :])  # 守恒量转化为基本量
    u_physics2 = consevation_to_physics1(u_conservation[1, :])
    u_physics3 = consevation_to_physics1(u_conservation[2, :])
    # E = u_physics[2, :]/(gamma-1)+0.5*u_physics[0, :] * \
    #   u_physics[1, :]*u_physics[1, :]  # 能量
    u11 = _draw(u_physics1, u_physics2, gamma, dx, dt, N=100)
    u22 = _draw(u_physics2, u_physics3, gamma, dx, dt, N=100)
    F1 = np.array([u22[0]*u22[1], u22[0]*u22[1]
                   * u22[1]+u22[2], (u22[2]+u22[2])*u22[1]])  # flux
    F2 = np.array([u11[0]*u11[1], u11[0]*u11[1]
                   * u11[1]+u11[2], (u11[2]+u11[2])*u11[1]])
    dF = np.zeros((3, 1))
    dF = (F1-F2)/dx  # 差分

    return dF
'''


def fluxDiscre(dx, dt, u_conservation):  # N+5

    FF = np.zeros((3, N+3))
    dFF = np.zeros((3, N+3))
    u_p = consevation_to_physics1(u_conservation)  # N+5
    # E = u_physics[2, :]/(gamma-1)+0.5*u_physics[0, :] * \
    #   u_physics[1, :]*u_physics[1, :]  # 能量

    u_pp = np.zeros((3, N+2))
    for i in range(0, N+2):
        u_pp[:, i] = _draw(u_p[:, i+1], u_p[:, i+2], gamma, dx, dt, 100)

    F1 = np.array([u_pp[0, :]*u_pp[1, :], u_pp[0, :]*u_pp[1, :]
                   * u_pp[1, :]+u_pp[2, :], (u_conservation[2, 1:-2]+u_pp[2, :])*u_pp[1, :]])  # flux
    dFF = (F1[:, 1:]-F1[:, :-1])/dx  # 差分

    return dFF


def consevation_to_physics(u_conservation):
    u_prim = np.zeros((3, N+1))
    u_prim[0, :] = u_conservation[0, :]
    u_prim[1, :] = u_conservation[1, :]/u_conservation[0, :]
    u_prim[2, :] = (u_conservation[2, :]-0.5*u_conservation[1, :]
                    * u_conservation[1, :]/u_conservation[0, :])*(gamma-1)

    return u_prim


def consevation_to_physics1(u_conservation):
    u_prim = np.zeros((3, N+5))
    u_prim[0, :] = u_conservation[0, :]
    u_prim[1, :] = u_conservation[1, :]/u_conservation[0, :]
    u_prim[2, :] = (u_conservation[2, :]-0.5*u_conservation[1, :]
                    * u_conservation[1, :]/u_conservation[0, :])*(gamma-1)

    return u_prim


'''
def consevation_to_physics1(u_conservation):
    u_prim = np.zeros((3, 1))
    u_prim[0, 0] = u_conservation[0, 0]
    u_prim[1, 0] = u_conservation[0, 1]/u_conservation[0, 0]
    u_prim[2, 0] = (u_conservation[0, 2]-0.5*u_conservation[0, 1]
                    * u_conservation[0, 1]/u_conservation[0, 0])*(gamma-1)

    return u_prim
'''


def _pressure_function(p, uu, gamma):
    '''
    Helper function, computing f_l(p,W_k), df_l(p,W_k)
    :param p: pressure
    :param uu: primitive state variables
    :param eos: equation of state
    :return:
    '''
    [rho_k, v_k, p_k] = uu
    if (p > p_k):  # shock

        A_k, B_k = 2 / ((gamma + 1) * rho_k), (gamma - 1) / (gamma + 1) * p_k

        f = (p - p_k) * np.sqrt(A_k / (p + B_k))

        df = np.sqrt(A_k / (B_k + p)) * (1 - (p - p_k) / (2 * (B_k + p)))

    else:  # rarefaction
        a_k = np.sqrt(gamma * p_k / rho_k)

        f = 2 * a_k / (gamma - 1) * ((p / p_k) **
                                     ((gamma - 1) / (2 * gamma)) - 1)

        df = 1 / (rho_k * a_k) * (p / p_k) ** (- (gamma + 1) / (2 * gamma))

    return f, df


def _pressure_function_init(uu_l, uu_r, gamma):
    '''
    guess the initial value of pressure at contact discontinuity
    :param uu_l: left state primitive variables
    :param uu_r: right state primitive variables
    :param gamma: equation of state

    '''
    [rho_l, v_l, p_l] = uu_l

    [rho_r, v_r, p_r] = uu_r

    return 0.5 * (p_l + p_r)


def _solve_contact_discontinuity(uu_l, uu_r, gamma):
    '''
    solve function f(p,W_l,W_r) = f_l(p,W_l) + f_r(p,w_r) + delta v
    :param p: pressure guess at contact discontinuity
    :param uu_l: left state primitive variables
    :param uu_r: right state primitive variables
    :param eos: equation of state
    :return: p,v around the contact discontinuity
    '''

    [rho_l, v_l, p_l] = uu_l
    [rho_r, v_r, p_r] = uu_r

    d_v = v_r - v_l

    MAX_ITE = 100
    TOLERANCE = 1.0e-12
    found = False

    p_old = _pressure_function_init(uu_l, uu_r, gamma)

    for i in range(MAX_ITE):
        f_l, df_l = _pressure_function(p_old, uu_l, gamma)
        f_r, df_r = _pressure_function(p_old, uu_r, gamma)

        p = p_old - (f_l + f_r + d_v) / (df_l + df_r)

        if(p < 0.0):
            p = TOLERANCE

        if (2 * abs(p - p_old) / (p + p_old) < TOLERANCE):
            found = True
            break
        p_old = p

    if not found:
        print('Divergence in Newton-Raphason iteration')

    v = 0.5 * (v_l + v_r + f_r - f_l)

    return p, v


def _sample_solution(p_m, v_m, uu_l, uu_r, gamma, s):
    '''
    :param p: pressure at contact discontinuity
    :param v: velocity at contact discontinuity
    :param uu_l: left primitive variables
    :param uu_r: right primitive variables
    :param gamma: equation of state
    :param s: sample point satisfies s = x/t
    :return: rho, v , p at (t,x)
    '''

    [rho_l, v_l, p_l] = uu_l
    a_l = np.sqrt(gamma * p_l / rho_l)

    [rho_r, v_r, p_r] = uu_r
    a_r = np.sqrt(gamma * p_r / rho_r)

    if s < v_m:
        # sampling point lies to the left of the contact discontinuity
        if p_m < p_l:

            # left rarefaction
            s_l = v_l - a_l
            a_ml = a_l * (p_m / p_l) ** ((gamma - 1) / (2 * gamma))
            s_ml = v_m - a_ml
            if s < s_l:  # left state
                rho, v, p = rho_l, v_l, p_l

            elif s < s_ml:  # left rarefaction wave
                rho = rho_l * (2 / (gamma + 1) + (gamma - 1) /
                               ((gamma + 1) * a_l) * (v_l - s)) ** (2 / (gamma - 1))
                v = 2 / (gamma + 1) * (a_l + (gamma - 1) / 2.0 * v_l + s)
                p = p_l * (2 / (gamma + 1) + (gamma - 1) / ((gamma + 1)
                                                            * a_l) * (v_l - s)) ** (2 * gamma / (gamma - 1))

            else:  # left contact discontinuity
                rho, v, p = rho_l * (p_m / p_l) ** (1 / gamma), v_m, p_m
        else:
            # left shock

            s_shock = v_l - a_l * \
                np.sqrt((gamma + 1) * p_m / (2 * gamma * p_l) +
                        (gamma - 1) / (2 * gamma))

            if s < s_shock:
                rho, v, p = rho_l, v_l, p_l
            else:
                rho = rho_l * (p_m / p_l + (gamma - 1) / (gamma + 1)) / \
                    ((gamma - 1) * p_m / ((gamma + 1) * p_l) + 1)
                v = v_m
                p = p_m
    else:
        # sampling point lies to the right of the contact discontinuity

        if p_m < p_r:

            # right rarefaction
            s_r = v_r + a_r
            a_mr = a_r * (p_m / p_r) ** ((gamma - 1) / (2 * gamma))
            s_mr = v_m + a_mr
            if s > s_r:  # left state
                rho, v, p = rho_r, v_r, p_r

            elif s > s_mr:  # left rarefaction wave
                rho = rho_r * (2 / (gamma + 1) - (gamma - 1) /
                               ((gamma + 1) * a_r) * (v_r - s)) ** (2 / (gamma - 1))
                v = 2 / (gamma + 1) * (-a_r + (gamma - 1) / 2.0 * v_r + s)
                p = p_r * (2 / (gamma + 1) - (gamma - 1) / ((gamma + 1)
                                                            * a_r) * (v_r - s)) ** (2 * gamma / (gamma - 1))

            else:  # left contact discontinuity
                rho, v, p = rho_r * (p_m / p_r) ** (1/gamma), v_m, p_m
        else:
            # right shock

            s_shock = v_r + a_r * \
                np.sqrt((gamma + 1) * p_m / (2 * gamma * p_r) +
                        (gamma - 1) / (2 * gamma))

            if s > s_shock:  # after shock
                rho, v, p = rho_r, v_r, p_r
            else:
                # preshock, contact discontinuity
                rho = rho_r * (p_m / p_r + (gamma - 1) / (gamma + 1)) / \
                    ((gamma - 1) * p_m / ((gamma + 1) * p_r) + 1)
                v = v_m
                p = p_m
    return rho, v, p


def _compute_domain(uu_l, uu_r, gamma, t):
    '''
    :param uu_l: left primitive variables
    :param uu_r: right primitive variables
    :param gamma: equation of state
    :param t: solution at time t
    return L, the computational domain [-L,L], which presents all feature of the solution
    '''

    [rho_l, v_l, p_l] = uu_l
    a_l = np.sqrt(gamma * p_l / rho_l)

    [rho_r, v_r, p_r] = uu_r
    a_r = np.sqrt(gamma * p_r / rho_r)

    return 2.0 * max(abs(v_l) + a_l, abs(v_r) + a_r) * t


def _draw(uu_l, uu_r, gamma, dx, dt, N):
    # draw the solution
    # _compute_domain(uu_l, uu_r, gamma, t)
    x = np.linspace(-dx, dx, num=N)

    uu = np.empty([1, 3], dtype=float)

    p_m, v_m = _solve_contact_discontinuity(uu_l, uu_r, gamma)

    uu = _sample_solution(p_m, v_m, uu_l, uu_r, gamma, 0)
    return uu
