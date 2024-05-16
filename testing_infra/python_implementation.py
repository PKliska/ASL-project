import csv
import sys

import numpy


def build_up_b(b, rho, dt, u, v, dx, dy):
    b[1:-1, 1:-1] = (rho * (1 / dt *
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) /
                             (2 * dx) + (v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) -
                            ((u[1:-1, 2:] - u[1:-1, 0:-2]) / (2 * dx)) ** 2 -
                            2 * ((u[2:, 1:-1] - u[0:-2, 1:-1]) / (2 * dy) *
                                 (v[1:-1, 2:] - v[1:-1, 0:-2]) / (2 * dx)) -
                            ((v[2:, 1:-1] - v[0:-2, 1:-1]) / (2 * dy)) ** 2))

    return b

def pressure_poisson(p, dx, dy, b, nit):
    pn = numpy.empty_like(p)
    pn = p.copy()

    for q in range(nit):
        pn = p.copy()
        p[1:-1, 1:-1] = (((pn[1:-1, 2:] + pn[1:-1, 0:-2]) * dy ** 2 +
                          (pn[2:, 1:-1] + pn[0:-2, 1:-1]) * dx ** 2) /
                         (2 * (dx ** 2 + dy ** 2)) -
                         dx ** 2 * dy ** 2 / (2 * (dx ** 2 + dy ** 2)) *
                         b[1:-1, 1:-1])

        p[:, -1] = p[:, -2]  # dp/dx = 0 at x = 2
        p[0, :] = p[1, :]  # dp/dy = 0 at y = 0
        p[:, 0] = p[:, 1]  # dp/dx = 0 at x = 0
        p[-1, :] = 0  # p = 0 at y = 2

    return p

def cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nx, ny, nit):
    un = numpy.empty_like(u)
    vn = numpy.empty_like(v)
    b = numpy.zeros((ny, nx))

    for n in range(nt):
        un = u.copy()
        vn = v.copy()

        b = build_up_b(b, rho, dt, u, v, dx, dy)
        p = pressure_poisson(p, dx, dy, b, nit)

        u[1:-1, 1:-1] = (un[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (un[1:-1, 1:-1] - un[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (un[1:-1, 1:-1] - un[0:-2, 1:-1]) -
                         dt / (2 * rho * dx) * (p[1:-1, 2:] - p[1:-1, 0:-2]) +
                         nu * (dt / dx ** 2 *
                               (un[1:-1, 2:] - 2 * un[1:-1, 1:-1] + un[1:-1, 0:-2]) +
                               dt / dy ** 2 *
                               (un[2:, 1:-1] - 2 * un[1:-1, 1:-1] + un[0:-2, 1:-1])))

        v[1:-1, 1:-1] = (vn[1:-1, 1:-1] -
                         un[1:-1, 1:-1] * dt / dx *
                         (vn[1:-1, 1:-1] - vn[1:-1, 0:-2]) -
                         vn[1:-1, 1:-1] * dt / dy *
                         (vn[1:-1, 1:-1] - vn[0:-2, 1:-1]) -
                         dt / (2 * rho * dy) * (p[2:, 1:-1] - p[0:-2, 1:-1]) +
                         nu * (dt / dx ** 2 *
                               (vn[1:-1, 2:] - 2 * vn[1:-1, 1:-1] + vn[1:-1, 0:-2]) +
                               dt / dy ** 2 *
                               (vn[2:, 1:-1] - 2 * vn[1:-1, 1:-1] + vn[0:-2, 1:-1])))

        u[0, :] = 0
        u[:, 0] = 0
        u[:, -1] = 0
        u[-1, :] = 1  # set velocity on cavity lid equal to 1
        v[0, :] = 0
        v[-1, :] = 0
        v[:, 0] = 0
        v[:, -1] = 0

    return u, v, p

if __name__ == '__main__':
    rho = 1
    nu = .1
    dt = .001
    nit = 50

    custom_dimension = sys.argv[2]
    if has_user_provided_custom_dimension := (custom_dimension is not None) and (custom_dimension != ""):
        ny = nx = int(custom_dimension)
    else:
        raise ValueError("Must provide an int value for matrix dimension (nx) as 2nd argument!")

    # custom_nt = sys.argv[3]
    # if has_user_provided_custom_nt := (custom_nt is not None) and (custom_nt != ""):
    #     nt = int(custom_nt)
    # else:
    nt = 100 # that's what we use as default arg in the C implementation

    # basically equivalent to this "dy = dx = 0.025" but we do it this way to keep it more bit aligned with the C implementations
    dy = dx = (0.025 * (nx-1)) / (nx - 1)

    u = numpy.zeros((ny, nx))
    v = numpy.zeros((ny, nx))
    p = numpy.zeros((ny, nx))
    b = numpy.zeros((ny, nx))

    u, v, p = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu, nx, ny, nit)

    custom_destination = sys.argv[1]
    if has_user_provided_destination := (custom_destination is not None) and (custom_destination != ""):
        output_destination = custom_destination
    else:
        output_destination = './.test_results/python_simulation_data.csv'

    # Output to CSV
    with open(output_destination, 'w') as csvfile:
        writer = csv.writer(csvfile, lineterminator='\n')

        # Write metadata (first line)l
        # writer.writerow([nx, ny, rho, nu])
        writer.writerow([f'{nx}', f'{ny}', f'{rho:.06f}', f'{nu:.06f}'] + ["" for _ in range(nx-3)])
        # writer.writerow([f'{elem:.06f}' for elem in [nx, ny, rho, nu]])

        def write_matrix(m):
            for row in m:
                # if ya wann know why ".17g" check out: https://stackoverflow.com/a/21162120/5338761
                writer.writerow([f'{elem:.6g}' for elem in row])

        # Write U matrix
        # writer.writerows(u)
        write_matrix(u)

        writer.writerow(["" for _ in range(nx+1)])

        # Write V matrix
        # writer.writerows(v)
        write_matrix(v)

        writer.writerow(["" for _ in range(nx+1)])

        # Write P matrix
        # writer.writerows(p)
        write_matrix(p)


    # print("Data saved to 'python_simulation_data.csv'")

