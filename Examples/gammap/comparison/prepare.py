#!/usr/bin/env python3
"""Write mass-bin AmpTools configs and PAWIAN configs/start parameters."""

from __future__ import annotations

import argparse
from pathlib import Path

from model import MASS_EDGES


def amptools_config(index: int, directory: Path) -> str:
    common = directory.parent
    data = (common / f'data-{index:02d}.root').resolve()
    mc = (common / f'mc-{index:02d}.root').resolve()
    normint = (directory / f'norm-{index:02d}.ni').resolve()
    return f"""fit zlm-bin-{index:02d}

reaction ksks Beam Proton KShort KShort
data ksks ROOTDataReader {data}
genmc ksks ROOTDataReader {mc}
accmc ksks ROOTDataReader {mc}
normintfile ksks {normint}

# Positive reflectivity: sqrt(1+P) Re[Ye^(-iPhi)] and
# sqrt(1-P) Im[Ye^(-iPhi)] are incoherent sums with shared coefficients.
sum ksks PosRe
sum ksks PosIm
amplitude ksks::PosRe::S0+ Zlm 0 0 1 1
amplitude ksks::PosIm::S0+ Zlm 0 0 -1 -1
amplitude ksks::PosRe::D2+ Zlm 2 2 1 1
amplitude ksks::PosIm::D2+ Zlm 2 2 -1 -1

constrain ksks::PosRe::S0+ ksks::PosIm::S0+
constrain ksks::PosRe::D2+ ksks::PosIm::D2+
initialize ksks::PosRe::S0+ cartesian 1.0 0.0 real
initialize ksks::PosIm::S0+ cartesian 1.0 0.0 real
initialize ksks::PosRe::D2+ cartesian 0.5 0.0
initialize ksks::PosIm::D2+ cartesian 0.5 0.0
"""


def prepare_amptools(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    for index in range(len(MASS_EDGES) - 1):
        (directory / f'bin-{index:02d}.cfg').write_text(
            amptools_config(index, directory),
            encoding='utf-8',
        )


def prepare_parameters(defaults: Path, output: Path) -> list[str]:
    s_stem = 'Reflectivity_f0(1500)_J0_M0_R+_K0'
    d_stem = 'Reflectivity_f2(1270)_J2_M2_R+_K0'
    free = {f'{d_stem}Mag', f'{d_stem}Phi'}
    fixed = []
    rows = []
    for raw in defaults.read_text(encoding='utf-8').splitlines():
        fields = raw.split()
        if not fields:
            continue
        name = fields[0]
        value = 0.0
        step = 0.0
        limits = ''
        if name == f'{s_stem}Mag':
            value = 1.0
        elif name == f'{d_stem}Mag':
            value, step = 0.5, 0.1
            limits = '\t0\t5'
        elif name == f'{d_stem}Phi':
            value, step = 0.35, 0.2
        elif 'ToKshort1_Kshort2Mag' in name:
            value = 1.0
        elif name.endswith('Scaling'):
            value, step = 1.0, 0.01
            limits = '\t0\t20'
        rows.append(f'{name}\t{value:.17g}\t{step:.17g}{limits}')
        if name not in free and not name.endswith('Scaling'):
            fixed.append(name)
    output.write_text('\n'.join(rows) + '\n', encoding='utf-8')
    return fixed


def pawian_config(
    index: int,
    directory: Path,
    parameters: Path,
    fixed: list[str],
) -> str:
    low, high = MASS_EDGES[index : index + 2]
    common = directory.parent
    fixed_lines = '\n'.join(f'mnParFix = {name}' for name in fixed)
    return f"""noOfThreads = 1
errLogMode = warning
datFile = {(common / 'data.txt').resolve()}
mcFile = {(common / 'mc.txt').resolve()}
unitInFile = GEV
orderInFile = E Px Py Pz
initial4Vec = 8.9382720813 0.0 0.0 8.0
projectile4Vec = 8.0 0.0 0.0 8.0
ratioMcToData = 1000
useDataEventWeight = true
useMCEventWeight = true
paramFile = {parameters.resolve()}
mode = pwa
verbose = false
productionFormalism = Reflectivity
usePolarization = true
cloneParticle = Kshort Kshort1
cloneParticle = Kshort Kshort2
finalStateParticle = proton
finalStateParticle = Kshort1
finalStateParticle = Kshort2
production = f0(1500) proton
production = f2(1270) proton
decay = Cano f0(1500) To Kshort1 Kshort2
decay = Cano f2(1270) To Kshort1 Kshort2
massRangeCuts = {low:.17g} {high:.17g} Kshort1 Kshort2
pdgTableFile = /Particle/pdtNew.table
cacheAmps = true
fixPrimaryChannelScaling = false
lmax = 3
{fixed_lines}
"""


def prepare_pawian(defaults: Path, directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    parameters = directory / 'start.dat'
    fixed = prepare_parameters(defaults, parameters)
    for index in range(len(MASS_EDGES) - 1):
        bin_directory = directory / f'bin-{index:02d}'
        bin_directory.mkdir(exist_ok=True)
        for stale_result in bin_directory.glob('finalResult_*.dat'):
            stale_result.unlink()
        (directory / f'bin-{index:02d}.cfg').write_text(
            pawian_config(index, directory, parameters, fixed),
            encoding='utf-8',
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest='command', required=True)
    amptools = subparsers.add_parser('amptools')
    amptools.add_argument('--directory', type=Path, required=True)
    pawian = subparsers.add_parser('pawian')
    pawian.add_argument('--defaults', type=Path, required=True)
    pawian.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'amptools':
        prepare_amptools(args.directory)
    else:
        prepare_pawian(args.defaults, args.directory)


if __name__ == '__main__':
    main()
