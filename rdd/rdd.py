import math
import sys


def rdd(radius, energy, S, r0,  alg):
    """
    energy: [MeV/nucleon]
    S:  [MeV*cm2/g]
    radius: [cm]
    r0
    alg:
    1 = Kraft-Scholz
    2 = Geiss
    returns dose in [MeV/g]
    """

    if alg == 1:
        rmax = rmax_scholz(energy)
    if alg == 2:
        rmax = rmax_geiss(energy)

    k = S / ((2 * math.pi) * (0.5 + math.log(rmax/r0)) * r0*r0)

    if radius <= r0:
        return k
    elif radius < rmax:
        return k * (r0 / radius) * (r0 / radius)
    else:
        return 0.0


def rmax_scholz(energy):
    """
    energy : [MeV/nucleon]
    returns maximum radius in [g/cm2]

    Eq. 2.5 in Scholz' habil. No good below 1 MeV/nucleon, according to Scholz.
    Originally in [um] here, in [cm]
    """

    return 5.0e-6 * math.pow(energy, 1.7)


def rmax_geiss(energy):
    """
    energy : [MeV/nucleon]
    returns maximum radius in [g/cm2]
    """

    return 4.0e-5 * energy * math.sqrt(energy)


def main(args=sys.argv[1:]):

    let = 260.0  # [MeV*cm2/g]
    r0 = 1e-6  # [g/cm2 (1e-6 = 10 nm track core size at rho = 1)]

    energy_p = 1.004   # [MeV/nucl]
    energy_he = 6.37   # [MeV/nucl]
    energy_c = 100.15  # [MeV/nucl]

    r = 1e-7  # [g/cm2]

    g = 1.602e-10  # [ 1 MeV/g in Gy ]

    while r < 1.0:

        d1 = g * rdd(r, energy_p, let, r0, 1)
        d2 = g * rdd(r, energy_he, let, r0, 1)
        d3 = g * rdd(r, energy_c, let, r0, 1)

        print("{:e} {:e} {:e} {:e}".format(r, d1, d2, d3))
        r *= 1.1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
