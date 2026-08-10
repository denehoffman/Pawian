# ruff: noqa: INP001
"""Shared laddu 0.21 channel and Zlm model for the three-way comparison."""

from __future__ import annotations

import math

import laddu as ld

POLARIZATION = 0.3519
POLARIZATION_ANGLE = 0.0
MASS_EDGES = tuple(round(1.0 + 0.1 * index, 1) for index in range(11))


def channel() -> ld.Channel:
    kaon = ld.particles.K_SHORT
    return ld.Channel(
        'gamma p -> K_S K_S p',
        edges=[
            ld.Edge(
                'beam',
                p4='beam',
                particle=ld.particles.PHOTON,
                output=True,
                initial_momentum=ld.InitialMomentum.energy(
                    8.0,
                    direction=[0.0, 0.0, 1.0],
                ),
            ),
            ld.Edge(
                'target',
                p4='target',
                particle=ld.particles.PROTON,
                output=True,
                initial_momentum=ld.InitialMomentum.momentum([0.0, 0.0, 0.0]),
            ),
            ld.Edge('X', mass_proposal=ld.MassProposal(2.0 * kaon.mass, high=2.0)),
            ld.Edge('recoil', p4='recoil', particle=ld.particles.PROTON, output=True),
            ld.Edge('ks1', p4='ks1', particle=kaon, output=True),
            ld.Edge('ks2', p4='ks2', particle=kaon, output=True),
        ],
        vertices=[
            ld.Vertex(
                'production',
                incoming=['beam', 'target'],
                outgoing=['X', 'recoil'],
                generation=ld.VertexProposal.t_exchange(
                    incoming='beam',
                    outgoing='X',
                    slope=3.0,
                    uniform_fraction=0.2,
                ),
            ),
            ld.Vertex(
                'decay',
                incoming=['X'],
                outgoing=['ks1', 'ks2'],
                generation=ld.VertexProposal.isotropic(),
            ),
        ],
    )


def angles(reaction: ld.Channel) -> tuple[ld.Expr, ld.Expr, ld.Expr]:
    """Match Examples/gammap/amptools/Zlm.cc's helicity-frame angles."""
    production = reaction.vertex('production')
    decay = reaction.vertex('decay')
    beam = production.vec3('beam')
    recoil = production.vec3('recoil')
    z_axis = -decay.vec3('recoil')
    y_hint = beam.cross(-recoil)
    costheta = decay.costheta('ks1', z_axis=z_axis, y_hint=y_hint)
    phi = decay.phi('ks1', z_axis=z_axis, y_hint=y_hint)
    epsilon = ld.Vec3(
        math.cos(POLARIZATION_ANGLE),
        math.sin(POLARIZATION_ANGLE),
        0.0,
    )
    production_normal = beam.unit().cross(-recoil.unit()).unit()
    big_phi = ld.atan2(
        production_normal.dot(epsilon),
        beam.unit().dot(epsilon.cross(production_normal)),
    )
    return costheta, phi, big_phi


def zlm_components(
    reaction: ld.Channel,
    ell: int,
    projection: int,
) -> tuple[ld.Expr, ld.Expr]:
    """Positive-reflectivity Zlm real/imaginary incoherent components."""
    costheta, phi, big_phi = angles(reaction)
    rotated = ld.spherical_harmonic(
        ell,
        projection,
        costheta=costheta,
        phi=phi,
    ) * ld.cis(-big_phi)
    return (
        math.sqrt(1.0 + POLARIZATION) * rotated.real(),
        math.sqrt(1.0 - POLARIZATION) * rotated.imag(),
    )


def truth_model(reaction: ld.Channel) -> ld.Model:
    """An f0(1500) S wave coherently interfering with an f2(1270) D wave."""
    s_amplitudes, d_amplitudes = truth_wave_amplitudes(reaction)
    intensity = sum(
        (s_amplitude + d_amplitude).norm_sqr()
        for s_amplitude, d_amplitude in zip(s_amplitudes, d_amplitudes, strict=True)
    )
    return ld.Model(intensity)


def truth_wave_amplitudes(
    reaction: ld.Channel,
) -> tuple[tuple[ld.Expr, ld.Expr], tuple[ld.Expr, ld.Expr]]:
    """Return the two incoherent S0+ and D2+ amplitudes used for generation."""
    s = reaction.s('X')
    kaon_mass = reaction.particle('ks1').mass
    f0 = ld.relativistic_breit_wigner(
        s,
        mass=1.500,
        width=0.109,
        mass1=kaon_mass,
        mass2=kaon_mass,
        l=0,
    )
    f2 = ld.relativistic_breit_wigner(
        s,
        mass=1.2754,
        width=0.1866,
        mass1=kaon_mass,
        mass2=kaon_mass,
        l=2,
    )
    d_coupling = 0.70 * ld.cis(0.65)
    s_real, s_imag = zlm_components(reaction, 0, 0)
    d_real, d_imag = zlm_components(reaction, 2, 2)
    return (
        (f0 * s_real, f0 * s_imag),
        (d_coupling * f2 * d_real, d_coupling * f2 * d_imag),
    )


def truth_wave_models(reaction: ld.Channel) -> tuple[ld.Model, ld.Model]:
    """Return diagonal S0+ and D2+ intensity models for truth projections."""
    return tuple(
        ld.Model(sum(amplitude.norm_sqr() for amplitude in amplitudes))
        for amplitudes in truth_wave_amplitudes(reaction)
    )


def bin_wave_models(reaction: ld.Channel) -> tuple[ld.Model, ld.Model]:
    """Return unit-coefficient S0+ and D2+ intensity models for fitted bins."""
    s_real, s_imag = zlm_components(reaction, 0, 0)
    d_real, d_imag = zlm_components(reaction, 2, 2)
    return (
        ld.Model(s_real * s_real + s_imag * s_imag),
        ld.Model(d_real * d_real + d_imag * d_imag),
    )


def bin_model(reaction: ld.Channel, prefix: str) -> ld.Model:
    """Mass-independent S0+ + D2+ angular model for one mass bin."""
    magnitude = ld.parameter(
        f'{prefix}_d_magnitude',
        initial=0.5,
        bounds=(0.0, 5.0),
    )
    phase = ld.parameter(
        f'{prefix}_d_phase',
        initial=0.0,
        bounds=(-math.pi, math.pi),
        periodic=True,
    )
    coefficient = ld.polar_complex(magnitude, phase)
    s_real, s_imag = zlm_components(reaction, 0, 0)
    d_real, d_imag = zlm_components(reaction, 2, 2)
    intensity = (s_real + coefficient * d_real).norm_sqr()
    intensity += (s_imag + coefficient * d_imag).norm_sqr()
    return ld.Model(intensity)
