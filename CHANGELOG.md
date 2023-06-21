# Changelog

### [0.6.2](https://www.github.com/bihealth/clinvar-tsv/compare/v0.6.1...v0.6.2) (2023-06-21)


### Bug Fixes

* failure to determine location only goes to debug level ([#22](https://www.github.com/bihealth/clinvar-tsv/issues/22)) ([76bc510](https://www.github.com/bihealth/clinvar-tsv/commit/76bc5105e6f0b26146dcaca43465c1dd59d1aee2))
* higher verbosity in Snakemake rules ([#26](https://www.github.com/bihealth/clinvar-tsv/issues/26)) ([a3fc327](https://www.github.com/bihealth/clinvar-tsv/commit/a3fc3271a8984b39d52c8852f6fe77b05b1193c0))
* interpret --cores argument ([#24](https://www.github.com/bihealth/clinvar-tsv/issues/24)) ([3b09803](https://www.github.com/bihealth/clinvar-tsv/commit/3b098038c47fb33bffabf94510cc9b6fee3f7d43))
* map "low penetrance" to "uncertain significance" ([#25](https://www.github.com/bihealth/clinvar-tsv/issues/25)) ([b2708d7](https://www.github.com/bihealth/clinvar-tsv/commit/b2708d75ad37d4270c253bf6928056c7deba8d84))
* no verbose output by default ([#27](https://www.github.com/bihealth/clinvar-tsv/issues/27)) ([0ad10cb](https://www.github.com/bihealth/clinvar-tsv/commit/0ad10cb8122480d2f46fbb6d2fba1e063be6da3c))
* reduce tqdm progress display unless on TTY ([#21](https://www.github.com/bihealth/clinvar-tsv/issues/21)) ([770b0c8](https://www.github.com/bihealth/clinvar-tsv/commit/770b0c833a2707a88e5ee9c0f2a0eb1435defdc6))

### [0.6.1](https://www.github.com/bihealth/clinvar-tsv/compare/v0.6.0...v0.6.1) (2023-06-21)


### Bug Fixes

* missing/problematic clinvar version ([#19](https://www.github.com/bihealth/clinvar-tsv/issues/19)) ([b11a8a4](https://www.github.com/bihealth/clinvar-tsv/commit/b11a8a435d9269031589106cf8929169893db5ef))

## [0.6.0](https://www.github.com/bihealth/clinvar-tsv/compare/v0.5.0...v0.6.0) (2023-06-21)


### Features

* allow providing clinvar version ([#17](https://www.github.com/bihealth/clinvar-tsv/issues/17)) ([dd80f2d](https://www.github.com/bihealth/clinvar-tsv/commit/dd80f2d10fceab350c61fa8de61bbe6264ad2008))


### Documentation

* adding badges to README ([#15](https://www.github.com/bihealth/clinvar-tsv/issues/15)) ([6e7ac01](https://www.github.com/bihealth/clinvar-tsv/commit/6e7ac013be4c0c7c34df41b69594581a7b8116f9))

## [0.5.0](https://www.github.com/bihealth/clinvar-tsv/compare/v0.4.1...v0.5.0) (2023-05-03)


### Features

* export structural variants ([#13](https://www.github.com/bihealth/clinvar-tsv/issues/13)) ([db44d87](https://www.github.com/bihealth/clinvar-tsv/commit/db44d8739f6f619266f806611950f339b0842352))

## 0.4.1

- Also writing out ``set_type`` column (#10).

## 0.4.0

- Greatly refining record merging strategy (#6).
  Also, providing both a ClinVar-like and a paranoid merging scheme.
- Improving CI (#7)

## 0.3.0

- Various refinements of the code.
- Adding tests and CI.

## 0.2.2

- Fixing bug with quotes.

## 0.2.1

- Fixing bug in setting clinical significance flags.

## 0.2.0

- Complete refurbishing of XML parsing, using models based on python-attrs.
- Removing old tests.

## 0.1.1

- Fixing installation of Snakefile.

## 0.1.0

- First actual release, versioning done using versioneer.
- Everything is new!
