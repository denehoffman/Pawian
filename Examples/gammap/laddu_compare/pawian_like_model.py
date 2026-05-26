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
from common import DECAY, POLARIZATION, REACTION


def x_decay_factor(spin: int, projection: Fraction | float) -> ld.Expression:
    return DECAY.canonical_factor(
        spin=spin,
        projection=projection,
        orbital_l=spin,
        coupled_spin=0,
        daughter='kshort1',
        daughter_1_spin=0,
        daughter_2_spin=0,
        lambda_1=0,
        lambda_2=0,
        frame='Helicity',
    )


def polar_parameter(mag_name: str, phi_name: str) -> ld.Expression:
    return ld.PolarComplexScalar(
        r=ld.parameter(mag_name),
        theta=ld.parameter(phi_name),
    )


def production_ls_factor(
    total_spin: Fraction,
    total_projection: Fraction,
    x_spin: int,
    x_projection: Fraction | int,
    recoil_spin: Fraction,
    orbital_l: int,
    coupled_spin: Fraction,
) -> ld.Expression:
    daughter_helicity = Fraction(x_projection) - recoil_spin
    if abs(float(total_projection)) > float(total_spin):
        return ld.Zero()
    if abs(float(daughter_helicity)) > float(total_spin):
        return ld.Zero()
    if abs(float(daughter_helicity)) > float(coupled_spin):
        return ld.Zero()
    return REACTION.production().canonical_factor(
        spin=total_spin,
        projection=total_projection,
        orbital_l=orbital_l,
        coupled_spin=coupled_spin,
        produced_spin=x_spin,
        recoil_spin=Fraction(1, 2),
        lambda_produced=x_projection,
        lambda_recoil=recoil_spin,
    )


def pawian_ls_wave(
    resonance: ld.Expression,
    total_projection: Fraction,
    recoil_spin: Fraction,
    x_spin: int,
    x_projection: Fraction | int,
    production_terms: tuple[tuple[int, Fraction, str, str], ...],
    decay_mag_name: str,
    decay_phi_name: str,
) -> ld.Expression:
    production = []
    for orbital_l, coupled_spin, mag_name, phi_name in production_terms:
        production.append(
            polar_parameter(mag_name, phi_name)
            * production_ls_factor(
                Fraction(1, 2),
                total_projection,
                x_spin,
                x_projection,
                recoil_spin,
                orbital_l,
                coupled_spin,
            )
        )
    decay_coeff = polar_parameter(decay_mag_name, decay_phi_name)
    return ld.expr_sum(production) * resonance * decay_coeff * x_decay_factor(x_spin, x_projection)


def pawian_ls_helicity_amplitude(
    photon_helicity: int,
    target_spin: Fraction,
    recoil_spin: Fraction,
) -> ld.Expression:
    total_projection = Fraction(photon_helicity, 1) + target_spin
    amplitudes = []
    # Pawian's current gammap cfg leaves these decays at WoDynamics.
    f0 = ld.One()
    f2 = ld.One()

    amplitudes.append(
        pawian_ls_wave(
            f0,
            total_projection,
            recoil_spin,
            0,
            0,
            (
                (
                    1,
                    Fraction(1, 2),
                    'L1S1/2_J1/2P-1C0Tof0(1500)_protonMag',
                    'L1S1/2_J1/2P-1C0Tof0(1500)_protonPhi',
                ),
            ),
            'L0S0_f0(1500)ToKshort1_Kshort2Mag',
            'L0S0_f0(1500)ToKshort1_Kshort2Phi',
        )
    )

    amplitudes.extend(
        [
            pawian_ls_wave(
                f2,
                total_projection,
                recoil_spin,
                2,
                x_projection,
                (
                    (
                        1,
                        Fraction(3, 2),
                        "L1S3/2_J1/2P-1C0Tof2'(1525)_protonMag",
                        "L1S3/2_J1/2P-1C0Tof2'(1525)_protonPhi",
                    ),
                    (
                        3,
                        Fraction(5, 2),
                        "L3S5/2_J1/2P-1C0Tof2'(1525)_protonMag",
                        "L3S5/2_J1/2P-1C0Tof2'(1525)_protonPhi",
                    ),
                ),
                "L2S0_f2'(1525)ToKshort1_Kshort2Mag",
                "L2S0_f2'(1525)ToKshort1_Kshort2Phi",
            )
            for x_projection in range(-2, 3)
        ]
    )
    return ld.expr_sum(amplitudes)


def pawian_ls_model() -> ld.Expression:
    """LS-production comparison model using laddu's canonical-factor machinery."""
    target_spins = [Fraction(-1, 2), Fraction(1, 2)]
    recoil_spins = [Fraction(-1, 2), Fraction(1, 2)]
    photon_helicities = [-1, 1]
    terms = []
    for target_spin in target_spins:
        for recoil_spin in recoil_spins:
            amps = {
                helicity: pawian_ls_helicity_amplitude(helicity, target_spin, recoil_spin)
                for helicity in photon_helicities
            }
            for helicity in photon_helicities:
                for helicity_prime in photon_helicities:
                    rho = ld.PhotonSDME(
                        helicity=helicity,
                        helicity_prime=helicity_prime,
                        polarization=POLARIZATION,
                    )
                    terms.append((rho * amps[helicity] * amps[helicity_prime].conj()).real())
    return ld.expr_sum(terms)
