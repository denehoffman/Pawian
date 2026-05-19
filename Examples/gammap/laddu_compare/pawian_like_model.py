#!/usr/bin/env python3
"""
A laddu model that keeps the photon SDME contraction explicit.

This is the comparison-side model, not Pawian production code. It is intended
to mirror the structure of gammapBaseLh closely enough that we can compare
normalizations and spin-convention choices before introducing a real
reflectivity-basis production formalism in Pawian.
"""

from __future__ import annotations

from fractions import Fraction

import laddu as ld


def reaction_variables() -> tuple[ld.Reaction, ld.Decay, ld.Mass, ld.Polarization]:
    beam = ld.Particle.stored('beam')
    target = ld.Particle.missing('target')
    kshort1 = ld.Particle.stored('kshort1')
    kshort2 = ld.Particle.stored('kshort2')
    kk = ld.Particle.composite('kk', (kshort1, kshort2))
    recoil = ld.Particle.stored('recoil')
    reaction = ld.Reaction.two_to_two(beam, target, kk, recoil)
    decay = reaction.decay('kk')
    polarization = reaction.polarization('pol_magnitude', 'pol_angle')
    return reaction, decay, decay.parent_mass(), polarization


def breit_wigner(tag: str, spin: int, mass_name: str, width_name: str, mass0: float) -> ld.Expression:
    _reaction, _decay, res_mass, _polarization = reaction_variables()
    return ld.BreitWigner(
        tag,
        mass=ld.parameter(mass_name, mass0),
        width=ld.parameter(width_name, initial=0.1, bounds=(0.001, 1.0)),
        l=spin,
        daughter_1_mass=ld.Mass(['kshort1']),
        daughter_2_mass=ld.Mass(['kshort2']),
        resonance_mass=res_mass,
    )


def x_decay_factor(tag: str, spin: int, projection: Fraction | float) -> ld.Expression:
    _reaction, decay, _res_mass, _polarization = reaction_variables()
    return decay.canonical_factor(
        tag,
        f'{tag}_D{spin}{projection}',
        spin=spin,
        projection=projection,
        orbital_l=spin,
        coupled_spin=0,
        daughter='kshort1',
        daughter_1_spin=0,
        daughter_2_spin=0,
        lambda_1=0,
        lambda_2=0,
    )


def wave_amplitude(tag: str, spin: int, resonance: ld.Expression, projection: Fraction | float) -> ld.Expression:
    coeff = ld.ComplexScalar(
        tag,
        f'{tag}_c',
        re=ld.parameter(f'{tag} re', initial=1.0),
        im=ld.parameter(f'{tag} im', initial=0.0),
    )
    return coeff * resonance * x_decay_factor(tag, spin, projection)


def helicity_amplitude(
    photon_helicity: int, target_spin: Fraction | float, recoil_spin: Fraction | float
) -> ld.Expression:
    projection = photon_helicity + target_spin - recoil_spin
    amplitudes = []

    f0 = breit_wigner('S0', 0, 'f0_mass', 'f0_width', 1.506)
    f2 = breit_wigner('D', 2, 'f2_mass', 'f2_width', 1.517)

    if projection == 0:
        amplitudes.append(wave_amplitude('S0', 0, f0, 0))
    if abs(float(projection)) <= 2:  # noqa: PLR2004
        amplitudes.append(wave_amplitude('D', 2, f2, projection))
    return ld.expr_sum(amplitudes) if amplitudes else ld.Zero()


def intensity_model() -> ld.Expression:
    """Return sum_spin rho(lambda,lambda') A_lambda A_lambda'^*."""
    _reaction, _decay, _res_mass, polarization = reaction_variables()
    target_spins = [Fraction(-1, 2), Fraction(1, 2)]
    recoil_spins = [Fraction(-1, 2), Fraction(1, 2)]
    photon_helicities = [-1, 1]

    terms = []
    for target_spin in target_spins:
        for recoil_spin in recoil_spins:
            amps = {helicity: helicity_amplitude(helicity, target_spin, recoil_spin) for helicity in photon_helicities}
            for helicity in photon_helicities:
                for helicity_prime in photon_helicities:
                    rho = ld.PhotonSDME(
                        f'rho_{helicity}_{helicity_prime}',
                        helicity=helicity,
                        helicity_prime=helicity_prime,
                        polarization=polarization,
                    )
                    terms.append((rho * amps[helicity] * amps[helicity_prime].conj()).real())
    return ld.expr_sum(terms)


def zlm_baseline_model() -> ld.Expression:
    """The standard laddu Zlm model for convention cross-checks."""
    _reaction, decay, res_mass, polarization = reaction_variables()
    angles = decay.angles('kshort1')
    z00p = ld.Zlm('S0+', 'Z00+', l=0, m=0, r='+', angles=angles, polarization=polarization)
    z22p = ld.Zlm('D2+', 'Z22+', l=2, m=2, r='+', angles=angles, polarization=polarization)
    bw_f0 = ld.BreitWigner(
        'S0+',
        mass=ld.parameter('f0_mass', 1.506),
        width=ld.parameter('f0_width', initial=0.112, bounds=(0.001, 1.0)),
        l=0,
        daughter_1_mass=ld.Mass(['kshort1']),
        daughter_2_mass=ld.Mass(['kshort2']),
        resonance_mass=res_mass,
    )
    bw_f2 = ld.BreitWigner(
        'D2+',
        mass=ld.parameter('f2_mass', 1.517),
        width=ld.parameter('f2_width', initial=0.086, bounds=(0.001, 1.0)),
        l=2,
        daughter_1_mass=ld.Mass(['kshort1']),
        daughter_2_mass=ld.Mass(['kshort2']),
        resonance_mass=res_mass,
    )
    s0p = ld.Scalar('S0+', value=ld.parameter('S0+ re', initial=1.0))
    d2p = ld.ComplexScalar(
        'D2+',
        re=ld.parameter('D2+ re', initial=0.5),
        im=ld.parameter('D2+ im', initial=0.5),
    )
    amp_re = s0p * bw_f0 * z00p.real() + d2p * bw_f2 * z22p.real()
    amp_im = s0p * bw_f0 * z00p.imag() + d2p * bw_f2 * z22p.imag()
    return amp_re.norm_sqr() + amp_im.norm_sqr()
