# Contributing

## Developer tools

To facilitate development, a set of developer tools is defined in the
[requirements-dev.txt](./requirements-dev.txt) file. You can install these with

```bash
python3 -m pip install -r requirements-dev.txt # --user
pre-commit install
```

The flag `--user` is optional. Currently, this installs the following:

- [`pre-commit`](https://pre-commit.com/): runs a set of checks on staged
  files. These checks are defined as Git hooks in the file
  [.pre-commit-config.yaml](./.pre-commit-config.yaml).

Another nice tool is [EditorConfig](https://editorconfig.org/), which is
available as a plugin to almost all editors. If you have activated this plugin
in your editor, it will use the [`.editorconfig`](./.editorconfig) file to
standardise some simple formatting conventions. This improves line-based
comparisons.


## Coding conventions

We use [Clang-Format](https://clang.llvm.org/docs/ClangFormat.html) to apply
formatting conventions to the C++ code. The conventions are defined in the
[`.clang-format` file](./.clang-format) and are based on the [LLVM coding
standards](https://llvm.org/docs/CodingStandards.html).
