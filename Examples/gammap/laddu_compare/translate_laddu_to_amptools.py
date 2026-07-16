#!/usr/bin/env python3
"""Convert laddu Kshort samples to AmpTools trees."""

from __future__ import annotations

import argparse
from array import array
from pathlib import Path

import ROOT
from common import FINAL_STATE, open_dataset


def write_amptools_tree(input_path: Path, output_path: Path, *, tree: str = 'kin') -> None:
    dataset = open_dataset(input_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_file = ROOT.TFile(str(output_path), 'RECREATE')
    tree = ROOT.TTree(tree, tree)
    n_fs = len(FINAL_STATE)

    NumFinalState = array('i', [n_fs])
    E_Beam = array('f', [0.0])
    Px_Beam = array('f', [0.0])
    Py_Beam = array('f', [0.0])
    Pz_Beam = array('f', [0.0])
    Weight = array('f', [0.0])
    E_FinalState = array('f', [0.0] * n_fs)
    Px_FinalState = array('f', [0.0] * n_fs)
    Py_FinalState = array('f', [0.0] * n_fs)
    Pz_FinalState = array('f', [0.0] * n_fs)

    tree.Branch('NumFinalState', NumFinalState, 'NumFinalState/I')
    tree.Branch('E_Beam', E_Beam, 'E_Beam/F')
    tree.Branch('Px_Beam', Px_Beam, 'Px_Beam/F')
    tree.Branch('Py_Beam', Py_Beam, 'Py_Beam/F')
    tree.Branch('Pz_Beam', Pz_Beam, 'Pz_Beam/F')
    tree.Branch('Weight', Weight, 'Weight/F')
    tree.Branch('E_FinalState', E_FinalState, 'E_FinalState[NumFinalState]/F')
    tree.Branch('Px_FinalState', Px_FinalState, 'Px_FinalState[NumFinalState]/F')
    tree.Branch('Py_FinalState', Py_FinalState, 'Py_FinalState[NumFinalState]/F')
    tree.Branch('Pz_FinalState', Pz_FinalState, 'Pz_FinalState[NumFinalState]/F')

    for event in dataset.events_local:
        p4_beam = event.p4('beam')
        E_Beam[0] = p4_beam.e
        Px_Beam[0] = p4_beam.px
        Py_Beam[0] = p4_beam.py
        Pz_Beam[0] = p4_beam.pz
        Weight[0] = event.weight
        for i, name in enumerate(FINAL_STATE):
            p4 = event.p4(name)
            E_FinalState[i] = p4.e
            Px_FinalState[i] = p4.px
            Py_FinalState[i] = p4.py
            Pz_FinalState[i] = p4.pz
        tree.Fill()

    output_file.Write()
    output_file.Close()

    print(f'wrote {output_path} with final-state order: {", ".join(FINAL_STATE)}')


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path, help='laddu parquet input')
    parser.add_argument('output', type=Path, help='Pawian ASCII output')
    args = parser.parse_args()
    write_amptools_tree(args.input, args.output)


if __name__ == '__main__':
    main()
